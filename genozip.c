
// ------------------------------------------------------------------
//   genozip.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>
#ifndef _WIN32
#include <execinfo.h>
#include <signal.h>
#endif
#include <getopt.h>

#include "genozip.h"
#include "text_help.h"
#include "version.h" // automatically incremented by the make when we create a new distribution
#include "txtfile.h"
#include "zfile.h"
#include "zip.h"
#include "piz.h"
#include "crypt.h"
#include "file.h"
#include "vblock.h"
#include "endianness.h"
#include "vcf.h"
#include "stream.h"
#include "url.h"
#include "strings.h"
#include "stats.h"
#include "arch.h"
#include "license.h"
#include "vcf.h"
#include "dict_id.h"
#include "reference.h"
#include "refhash.h"
#include "context.h"
#include "random_access.h"
#include "codec.h"

// globals - set it main() and never change
const char *global_cmd = NULL; 
ExeType exe_type;

// primary_command vs command: primary_command is what the user typed on the command line. command is what is 
// running now - for example, when ZIP is unzipping a reference, primary_command=ZIP and command=PIZ
CommandType command = NO_COMMAND, primary_command = NO_COMMAND; 
const char *command_line = NULL;

uint32_t global_max_threads = DEFAULT_MAX_THREADS; 
uint32_t global_max_memory_per_vb = 0; // ZIP only: used for reading text file data

static void print_call_stack (void) 
{
#ifndef _WIN32
#   define STACK_DEPTH 15
    void *array[STACK_DEPTH];
    size_t size = backtrace(array, STACK_DEPTH);
    
    fprintf (stderr, "Call stack:\n");
    backtrace_symbols_fd (array, size, STDERR_FILENO);
#endif
}

static bool im_in_main_exit = false, exit_completed = false;
void main_exit (bool show_stack, bool is_error) 
{
    im_in_main_exit = true;

    if (false /* show_stack */) print_call_stack(); //this is useless - doesn't print function names

    buf_test_overflows_all_vbs("exit_on_error");

    url_kill_curl();
    file_kill_external_compressors(); 

    // if we're in ZIP - remove failed genozip file (but don't remove partial failed text file in PIZ - it might be still useful to the user)
    if (primary_command == ZIP && z_file && z_file->name && !flag.reading_reference) {
        char save_name[strlen (z_file->name)+1];
        strcpy (save_name, z_file->name);

        // if we're not the main thread - cancel the main thread before closing z_file, so that the main 
        // thread doesn't attempt to access it (eg. z_file->data_type) and get a segmentation fault.
        if (!arch_am_i_io_thread()) 
            cancel_io_thread(); 

        file_close (&z_file, false); // also frees file->name

        // note: logic to avoid a race condition causing the file not to be removed - if another thread seg-faults
        // because it can't access a z_file filed after z_file is freed, and main_sigsegv_handler aborts
        file_remove (save_name, true);
    }

    exit_completed = true;

    if (is_error)
        abort();
    else
        exit (0);
} 

#ifndef _WIN32
static void main_sigsegv_handler (int sig) 
{
    // note: during exit_on_error, we close z_file which might cause compute threads access z_file fields to 
    // throw a segmentation fault. we don't show it in this case, as we have already displayed the error itself
    if (!im_in_main_exit) 
        fprintf (stderr, "\nError: Segmentation fault\n");

    // busy-wait for exit_on_error to complete before exiting cleanly
    else {
        //fprintf (stderr, "Segmentation fault might appear now - it is not an error - you can safely ignore it\n");
        while (!exit_completed) 
            usleep (10000); // 10 millisec
        exit (1);
    }

    //print_call_stack(); //this is useless - doesn't print function names
    abort();
}
#endif

static void main_print_help (bool explicit)
{
    static const char **texts[NUM_EXE_TYPES+1] = {help_genozip, help_genounzip, help_genols, help_genocat, help_genozip_developer}; // same order as ExeType
    static unsigned sizes[] = {sizeof(help_genozip), sizeof(help_genounzip), sizeof(help_genols), sizeof(help_genocat), sizeof(help_genozip_developer)};
    
    if (flag.force) exe_type = NUM_EXE_TYPES; // -h -f shows developer help

    str_print_text (texts[exe_type], sizes[exe_type] / sizeof(char*), 
                    flag.force ? "                          " : "                     ",  "\n", 0);
    str_print_text (help_footer, sizeof(help_footer) / sizeof(char*), "", "\n", 0);

// in Windows, we ask the user to click a key - this is so that if the user double clicks on the EXE
// from Windows Explorer - the terminal will open and he will see the help
#ifdef _WIN32
    if (!explicit) {
        printf ("Press any key to continue...\n");
        getc(stdin);
    }
#endif
}

static void main_print_version()
{
    printf ("version=%s\n", GENOZIP_CODE_VERSION);  
}

static void main_list_dir(); // forward declaration

static void main_genols (const char *z_filename, bool finalize, const char *subdir, bool recursive) 
{
    if (!finalize) {
        // no specific filename = show entire directory
        if (!z_filename) {
            main_list_dir("."); 
            return;
        }

        char *last_c = (char *)&z_filename[strlen(z_filename)-1];
        if (*last_c == '/' 
#ifdef _WIN32
            || *last_c == '\\'
#endif
            ) {
            *last_c = 0; // remove trailing '/' or '\'
            main_list_dir (z_filename);
            return;
        }
 
        // filename is a directory - show directory contents (but not recursively)
        if (!subdir && file_is_dir (z_filename)) {
            main_list_dir (z_filename);
            return;
        }
    }

    static bool first_file = true;
    static unsigned files_listed=0, files_ignored=0;
    static int64_t total_uncompressed_len=0, total_compressed_len=0;

    const unsigned FILENAME_WIDTH = 40;

    const char *head_format = "\n%-5.5s %11s %10s %10s %6s %s  %*s %s\n";
    const char *foot_format = "\nTotal:                  %10s %10s %5.*fX\n";
    const char *item_format = "%-5.5s %11s %10s %10s %5.*fX  %s  %s%s%*s %s\n";

    const char *head_format_bytes = "\n%-5.5s %11s %15s %15s %6s  %*s\n";
    const char *foot_format_bytes = "\nTotal:            %15s %15s %5.*fX\n";
    const char *item_format_bytes = "%-5.5s %11s %15s %15s %5.*fX  %s%s%*s\n";

    // we accumulate the string in str_buf and print in the end - so it doesn't get mixed up with 
    // warning messages regarding individual files
    static Buffer str_buf = EMPTY_BUFFER; 
    
    if (finalize) {
        if (files_listed > 1) {
            double ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            if (flag.bytes) 
                bufprintf (evb, &str_buf, foot_format_bytes, str_int_s (total_compressed_len).s, 
                           str_int_s (total_uncompressed_len).s, ratio < 100, ratio)
            else 
                bufprintf (evb, &str_buf, foot_format, str_size(total_compressed_len).s, 
                           str_size(total_uncompressed_len).s, ratio < 100, ratio);
        }
        
        ASSERTW (!files_ignored, "Ignored %u file%s that %s not have a " GENOZIP_EXT " extension", 
                 files_ignored, files_ignored==1 ? "" : "s", files_ignored==1 ? "does" : "do");
        
        goto finish;
    }

    if (first_file) {
        if (flag.bytes) 
            bufprintf (evb, &str_buf, head_format_bytes, "Type", "Records", "Compressed", "Original", "Factor", -(int)FILENAME_WIDTH, "Name")
        else
            bufprintf (evb, &str_buf, head_format, "Type", "Records", "Compressed", "Original", "Factor", " MD5 of original textual file    ", -(int)FILENAME_WIDTH, "Name", "Creation");
        
        first_file = false;
    }
    
    if (!file_has_ext (z_filename, GENOZIP_EXT) || access (z_filename, F_OK)!=0) {
        files_ignored++;
        goto finish;
    }

    bool is_subdir = subdir && (subdir[0] != '.' || subdir[1] != '\0');

    z_file = file_open (z_filename, READ, Z_FILE, 0); // open global z_file

    Md5Hash md5_hash_bound; 
    uint64_t txt_data_size, num_lines;
    char created[FILE_METADATA_LEN];
    if (!zfile_read_genozip_header (&md5_hash_bound, &txt_data_size, &num_lines, created))
        goto finish;

    double ratio = z_file->disk_size ? ((double)txt_data_size / (double)z_file->disk_size) : 0;
    
    // TODO: have an option to print ref_file_name and ref_file_md5

    DataType dt = (z_file->data_type == DT_SAM && (z_file->flags & SEC_GENOZIP_HEADER_FL_TXT_IS_BIN)) ? DT_BAM : z_file->data_type;

    if (flag.bytes) 
        bufprintf (evb, &str_buf, item_format_bytes, dt_name (dt), str_uint_commas (num_lines).s, 
                   str_int_s (z_file->disk_size).s, str_int_s (txt_data_size).s, ratio < 100, ratio, 
                   (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
                   is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH,
                   z_filename)
    else 
        bufprintf (evb, &str_buf, item_format, dt_name (dt), str_uint_commas (num_lines).s,
                   str_size (z_file->disk_size).s, str_size (txt_data_size).s, ratio < 100, ratio, 
                   md5_display (md5_hash_bound).s,
                   (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
                   is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH,
                   z_filename, created);

    total_compressed_len   += z_file->disk_size;
    total_uncompressed_len += txt_data_size;
    
    files_listed++;

    // if --unbind, OR if the user did genols on one file (not a directory), show bound components, if there are any
    if (flag.unbind || (!flag.multiple_files && !recursive && z_file->num_components >= 2)) {
        buf_add_string (evb, &str_buf, "Components:\n");
        const SectionListEntry *sl_ent = NULL;
        uint64_t num_lines_count=0;
        while (sections_get_next_section_of_type (&sl_ent, SEC_TXT_HEADER, false, false)) {
            zfile_read_section_header (evb, sl_ent->offset, sl_ent->vblock_i, SEC_TXT_HEADER);

            SectionHeaderTxtHeader *header = FIRSTENT (SectionHeaderTxtHeader, evb->compressed);
            
            num_lines_count += BGEN64 (header->num_lines);
            bufprintf (evb, &str_buf, item_format, "", str_uint_commas (BGEN64 (header->num_lines)).s, "", 
                       str_size (BGEN64 (header->txt_data_size)).s, 
                       0, 0.0, md5_display (header->md5_hash_single), "", "",
                       -(int)FILENAME_WIDTH, header->txt_filename, "");

            buf_free (&evb->compressed);
        }

        if (num_lines_count != num_lines)
            buf_add_string (evb, &str_buf, "\nNote: the difference between the file's Records and the total of its components' is the number of lines of the 1st component's header\n");
    }
    file_close (&z_file, false);

finish:
    if (!recursive) {
        buf_print (&str_buf, false);
        buf_free (&str_buf);
    }
}

// if this is a bound file, and we don't have --unbind or --force, we ask the user
static void main_ask_about_unbind (void)
{
    if (!isatty(0) || !isatty(2)) return; // if we stdin or stderr is redirected - we cannot ask the user an interactive question
    if (flag.test || flag.genocat_info_only) return; // other cases we don't ask test

    fprintf (stderr, "\n%s: %s contains %u bound files. You may either:\n"
                     "y) uncompress and unbind - retrieve the individual files (or use --unbind=<prefix> to add a prefix) ; or-\n"
                     "n) uncompress to one large file\n"
                     "Note: in the future, you may silence this question with --force\n\n", 
                     global_cmd, z_name, z_file->num_components);
        
    // read all chars available on stdin, so that if we're processing multiple files - and we ask this question
    // for a subsequent file later - we don't get left overs of this response
    char read_buf[1000];
    str_query_user ("Do you wish to unbind this file? ([y] or n) ", read_buf, sizeof(read_buf), str_verify_y_n, "Y");

    if (read_buf[0] == 'Y') flag.unbind = ""; // unbind with no prefix
    fprintf (stderr, "\n");
}

static void main_genounzip (const char *z_filename, const char *txt_filename, bool is_last_file)
{
    // save flag as it might be modified - so that next file has the same flags
    SAVE_FLAGS;

    txtfile_header_initialize();
    
    // get input FILE
    ASSINP0 (z_filename, "Error: z_filename is NULL");

    // we cannot work with a remote genozip file because the decompression process requires random access
    ASSINP (!url_is_url (z_filename), 
            "%s: genozip files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    ASSINP (!txt_filename || !url_is_url (txt_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, txt_filename);

    // skip this file if its size is 0
    RETURNW (file_get_size (z_filename),, "Cannot decompress file %s because its size is 0 - skipping it", z_filename);

    z_file = file_open (z_filename, READ, Z_FILE, DT_NONE);    

    // read the genozip header:
    // 1) verify the data type deduced from the file name, or set the data type if it wasn't deduced
    // 2) if an external reference is not specified, check if the file needs one, and if it does - set it from the header
    // 3) identify skip cases (DT_NONE returned) - empty file, unzip of a reference
    // 4) reset flag.unbind if file contains only one component
    if (!zfile_read_genozip_header (0,0,0,0)) goto done; 

    // if we're genocatting a BAM file, output it as a SAM unless user requested otherwise
    if (z_file->data_type == DT_SAM && (z_file->flags & SEC_GENOZIP_HEADER_FL_TXT_IS_BIN) && flag.out_dt==-1 && exe_type == EXE_GENOCAT)
        flag.out_dt = DT_SAM;

    // if this is a bound file, and we don't have --unbind or --force, we ask the user
    if (z_file->num_components >= 2 && !flag.unbind && !flag.force)
        main_ask_about_unbind();

    // case: reference not loaded yet bc --reference wasn't specified, and we got the ref name from zfile_read_genozip_header()   
    if (flag.reference == REF_EXTERNAL && !ref_is_reference_loaded()) {
        ASSINP (flag.reference != REF_EXTERNAL || !flag.show_ref_seq, "%s: Error: --show-ref-seq cannot be used on a file that requires a reference file: use genocat --show-ref-seq on the reference file itself instead", global_cmd);

        if (!flag.genocat_info_only) {
            SAVE_VALUE (z_file); // actually, read the reference first
            ref_load_external_reference (false, false /* parameter ignored for REF_EXTERNAL */);
            RESTORE_VALUE (z_file);
        }
    }

    flags_update_piz_one_file();
    
    // set txt_filename from genozip file name (inc. extensions if translating or --bgzip)
    if (!txt_filename && !flag.to_stdout && !flag.unbind) 
        txt_filename = txtfile_piz_get_filename (z_filename, "", true);

    // open output txt file (except if unbinding or outputting to stdout)
    if (txt_filename || flag.to_stdout) {
        ASSERT0 (!txt_file, "Error: txt_file is unexpectedly already open"); // note: in bound mode, we expect it to be open for 2nd+ file
        txt_file = file_open (txt_filename, WRITE, TXT_FILE, flag.out_dt);
    }
    else if (flag.unbind) {
        // do nothing - the component files will be opened by txtfile_genozip_to_txt_header()
    }
    else {
        ABORT0 ("Error: unrecognized configuration for the txt_file");
    }
    
    // a loop for decompressing all components in unbind mode. in non-unbind mode, it collapses to one a single iteration.
    uint32_t component_i=0;
    while (piz_one_file (component_i, is_last_file) && flag.unbind) component_i++;

    if (!flag.to_stdout && !flag.unbind) {
        // don't close stdout - we might still need it for the next file
        // don't close in unbind mode - piz_one_file() opens and closes each component
        file_close (&txt_file, false); 
    }

    file_close (&z_file, false);

    if (flag.replace && txt_filename && z_filename) file_remove (z_filename, true); 

done:
    RESTORE_FLAGS;
}

// run the test genounzip after genozip - for the most reliable testing that is nearly-perfectly indicative of actually 
// genounzipping, we create a new genounzip process
static void main_test_after_genozip (char *exec_name, char *z_filename, bool is_last_file)
{
    const char *password = crypt_get_password();

    // is we have a loaded reference and it is no longer needed, unload it now, to free the memory for the testing process
    if (is_last_file) ref_unload_reference (true);

    StreamP test = stream_create (0, 0, 0, 0, 0, 0, 0,
                                  "To use the --test option",
                                  exec_name, "--decompress", "--test", z_filename,
                                  flag.quiet       ? "--quiet"        : SKIP_ARG,
                                  password         ? "--password"     : SKIP_ARG,
                                  password         ? password         : SKIP_ARG,
                                  flag.show_memory ? "--show-memory"  : SKIP_ARG,
                                  flag.show_time   ? "--show-time"    : SKIP_ARG,
                                  flag.threads_str ? "--threads"      : SKIP_ARG,
                                  flag.threads_str ? flag.threads_str : SKIP_ARG,
                                  flag.reference == REF_EXTERNAL ? "--reference" : SKIP_ARG,
                                  flag.reference == REF_EXTERNAL ? ref_filename  : SKIP_ARG,
                                  NULL);

    // wait for child process to finish, so that the shell doesn't print its prompt until the test is done
    int exit_code = stream_wait_for_exit (test);

    TEMP_VALUE (primary_command, TEST_AFTER_ZIP); // make exit_on_error NOT delete the genozip file in this case, so its available for debugging
    ASSINP (!exit_code, "genozip test exited with status %d\n", exit_code);
    RESTORE_VALUE (primary_command); // recover in case of more non-concatenated files
}

static void main_genozip_open_z_file (char **z_filename)
{
    DataType z_data_type = txt_file->data_type;

    if (!(*z_filename) && !flag.to_stdout) {
        bool is_url = url_is_url (txt_file->name);
        const char *basename = is_url ? file_basename (txt_file->name, false, "", 0,0) : NULL;
        const char *local_txt_filename = basename ? basename : txt_file->name;

        unsigned fn_len = strlen (local_txt_filename);
        *z_filename = (char *)MALLOC (fn_len + 30); // add enough the genozip extension e.g. 23andme.genozip

        // if the file has an extension matching its type, replace it with the genozip extension, if not, just add the genozip extension
        const char *genozip_ext = file_exts[file_get_z_ft_by_txt_in_ft (txt_file->data_type, txt_file->type)];

        if (file_has_ext (local_txt_filename, file_exts[txt_file->type]))
            sprintf (*z_filename, "%.*s%s", (int)(fn_len - strlen (file_exts[txt_file->type])), local_txt_filename, genozip_ext); 
        else 
            sprintf (*z_filename, "%s%s", local_txt_filename, genozip_ext); 

        FREE (basename);
    }
    else if (flag.to_stdout) { // stdout
#ifdef _WIN32
        // this is because Windows redirection is in text (not binary) mode, meaning Windows edits the output stream...
        ASSINP (isatty(1), "%s: redirecting binary output is not supported on Windows, use --output instead", global_cmd);
#endif
        ASSINP (flag.force || !isatty(1), "%s: you must use --force to output a compressed file to the terminal", global_cmd);

        ASSERT0 (!flag.pair, "Error: cannot uncompress both paired end files while redirecting the output");
    } 

    z_file = file_open (*z_filename, flag.pair ? WRITEREAD : WRITE, Z_FILE, z_data_type);

    // note on BCF and CRAM: we used bcftools/samtools as an external compressor, so that genozip sees the text,
    // not binary, data of these files - the same as if the file were compressed with eg bz2
    if (z_file->data_type == DT_BAM) z_file->flags |= SEC_GENOZIP_HEADER_FL_TXT_IS_BIN; // compressed file is stored in binary form, and sizes are of the binary file
}

static void main_genozip (const char *txt_filename, 
                          char *z_filename,
                          bool is_first_file, bool is_last_file,
                          char *exec_name)
{
    SAVE_FLAGS;

    license_get(); // ask the user to register if she doesn't already have a license (note: only genozip requires registration - unzip,cat,ls do not)

    ASSINP (!z_filename || !url_is_url (z_filename), 
            "%s: output files must be regular files, they cannot be a URL: %s", global_cmd, z_filename);

    // get input file
    if (!txt_file) // open the file - possibly already open from main_load_reference
        txt_file = file_open (txt_filename, READ, TXT_FILE, DT_NONE); 

    // skip this file if its size is 0
    RETURNW (txt_file,, "Cannot compress file %s because its size is 0 - skipping it", txt_filename);

    flag.to_stdout = !txt_filename && !z_filename; // implicit setting of stdout by using stdin, unless -o was used

    if (!txt_file->file) goto done; // this is the case where multiple files are given in the command line, but this one is not compressible - we skip it

    stats_add_txt_name (txt_name);

    ASSERT0 (flag.bind || !z_file, "Error: expecting z_file to be NULL in non-bound mode");

    // get output FILE
    if (!z_file)  // skip if we're the second file onwards in bind mode - nothing to do
        main_genozip_open_z_file (&z_filename);
    
    flags_update_zip_one_file();

    zip_one_file (txt_file->basename, is_last_file);

    if (flag.show_stats && (!flag.bind || is_last_file)) stats_display();

    bool remove_txt_file = z_file && flag.replace && txt_filename;

    file_close (&txt_file, !is_last_file);

    if ((is_last_file || !flag.bind) && !flag.to_stdout && z_file) 
        file_close (&z_file, !is_last_file); 

    if (remove_txt_file) file_remove (txt_filename, true); 

    // test the compression, if the user requested --test
    if (flag.test && (!flag.bind || is_last_file)) main_test_after_genozip (exec_name, z_filename, is_last_file);

done:
    RESTORE_FLAGS;
}

static void main_list_dir(const char *dirname)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir (dirname);
    ASSINP (dir, "Error: failed to open directory: %s", strerror (errno));

    int ret = chdir (dirname);
    ASSINP (!ret, "Error: failed to chdir(%s)", dirname);

    while ((ent = readdir(dir))) 
        if (!file_is_dir (ent->d_name))  // don't go down subdirectories recursively
            main_genols (ent->d_name, false, dirname, true);
    
    closedir(dir);    

    ret = chdir ("..");
    ASSINP0 (!ret, "Error: failed to chdir(..)");
}

static inline DataType main_get_file_dt (const char *filename)
{   
    switch (command) {
        case ZIP: return txtfile_get_file_dt (filename);
        case PIZ: return zfile_get_file_dt (filename);
        default : return DT_NONE;
    }
}

static int main_sort_input_filenames (const void *fn1, const void *fn2)
{
    DataType dt1 = main_get_file_dt (*(char **)fn1);
    DataType dt2 = main_get_file_dt (*(char **)fn2);
    
    bool use_refhash1 = (dt1 == DT_FASTQ);
    bool use_refhash2 = (dt2 == DT_FASTQ);

    // refhash users - at the end of the list
    if (use_refhash1 != use_refhash2) return (int)use_refhash1 - (int)use_refhash2;

    // within refhash users, and within refhash non-users, sort by data type
    if (dt1 != dt2) return (int)dt1 - (int)dt2;

    // within files of the same data type, keep original order
    return fn1 - fn2;
}

static void main_load_reference (const char *filename, bool is_first_file, bool is_last_file)
{
    if (flag.reference != REF_EXTERNAL && flag.reference != REF_EXT_STORE) return;

    if (exe_type == EXE_GENOCAT && flag.show_stats) return; // we don't need the reference if we're just showing stats (ignore here if user provided it)
    
    int old_ref_use_aligner = flag.ref_use_aligner;
    DataType dt = main_get_file_dt (filename);
    flag.ref_use_aligner = (old_ref_use_aligner || dt == DT_FASTQ || dt == DT_FASTA) && primary_command == ZIP;

    // no need to load the reference if genocat just wants to see some sections (unless its genocat of the refernece file itself)
    if (flag.genocat_info_only && dt != DT_REF) return;

    // we also need the aligner if this is an unaligned SAM 
    if (!flag.ref_use_aligner && dt==DT_SAM && primary_command==ZIP) {

        // open here instead of in main_genozip
        txt_file = file_open (filename, READ, TXT_FILE, 0);

        // use the aligner if over 5 of the 100 first lines of the file are unaligned
        flag.ref_use_aligner = txt_file && txt_file->file && txtfile_test_data ('@', 100, 0.05, sam_zip_is_unaligned_line); 
    }

    RESET_VALUE (txt_file); // save and reset - for use by reference loader

    if (is_first_file)
        ref_load_external_reference (false, is_last_file); // also loads refhash if needed

    // Read the refhash and calculate the reverse compliment genome for the aligner algorithm - it was not used before and now it is
    else if (!old_ref_use_aligner && flag.ref_use_aligner) { 
        ref_generate_reverse_complement_genome();
        refhash_load_standalone();
    }

    RESTORE_VALUE (txt_file);
}

static void main_copy_command_line (int argc, char **argv)
{
    unsigned len=0, pw_len=0;
    const char *pw=0;

    for (int i=0; i < argc; i++)
        len += strlen (argv[i]) + 1; // +1 for seperator (' ' or '\0')

    command_line = CALLOC (len);

    if ((pw = crypt_get_password())) pw_len  = strlen (pw);

    for (int i=0; i < argc; i++) {

        unsigned arg_len = strlen (argv[i]);

        if (pw && !strcmp(argv[i], pw)) // "-p 123", "--pass 123" etc
            sprintf ((char*)&command_line[strlen(command_line)], "***%s", (i < argc-1 ? " ": "")); // hide password

        else if (pw && (arg_len >= pw_len + 2) &&  // check for -p123 or eg -fmp123
                    !strcmp (&argv[i][arg_len-pw_len], pw) && // not air-tight test, but good enough (eg "-ofilenamep123" will incorrectly trigger)
                    argv[i][0] == '-' &&
                    argv[i][arg_len-pw_len-1] == 'p')
            sprintf ((char*)&command_line[strlen(command_line)], "%.*s***%s", arg_len-pw_len, argv[i], (i < argc-1 ? " ": "")); // hide password

        else
            sprintf ((char*)&command_line[strlen(command_line)], "%s%s", argv[i], (i < argc-1 ? " ": ""));
    }
}

void TEST()
{
    FILE *fp = fopen ("seq", "rb");
    char *data = MALLOC (14236798);
    (void)!fread (data, 1, 14236798, fp); // (void)! to quieten compiler "warning: ignoring return value"
}


int main (int argc, char **argv)
{
    arch_initialize();
    buf_initialize(); 
    vb_initialize_evb();
    random_access_initialize();
    codec_initialize();

    //TEST();exit(0);

#ifdef _WIN32
    // lowercase argv[0] to allow case-insensitive comparison in Windows
    str_tolower (argv[0], argv[0]);
#else
    signal (SIGSEGV, main_sigsegv_handler);   // segmentation fault handler
#endif

    if      (strstr (argv[0], "genols"))    exe_type = EXE_GENOLS;
    else if (strstr (argv[0], "genocat"))   exe_type = EXE_GENOCAT;
    else if (strstr (argv[0], "genounzip")) exe_type = EXE_GENOUNZIP;
    else                                    exe_type = EXE_GENOZIP; // default
    
    global_cmd = file_basename (argv[0], true, "(executable)", NULL, 0); // global var

    bool is_short[256] = { 0 }; // indexed by character of short option.
    flags_init_from_command_line (argc, argv, is_short);
    main_copy_command_line (argc, argv); // can only be called after --password is processed

    // if command not chosen explicitly, use the default determined by the executable name
    if (command < 0) { 

        if (exe_type == EXE_GENOLS) command = LIST; // genols can be run without arguments
        
        // genozip with no input filename, no output filename, and no output or input redirection 
        // note: in docker stdin is a pipe even if going to a terminal. so we show the help even if
        // coming from a pipe. the user must use "-" to redirect from stdin
        else if (command == -1 && optind == argc && !flag.out_filename && 
                 (isatty(0) || arch_am_i_in_docker()) && isatty(1)) {
            // case: --register
            if (flag.do_register) {
                license_get();
                exit (0);
            }

            // case: requesting to display the reference: genocat --reference <ref-file> and optionally --regions
            if (exe_type == EXE_GENOCAT && flag.reference) 
                ref_load_external_reference (true, true);

            // otherwise: show help
            else
                main_print_help (false);

            return 0;
        }

        else if (exe_type == EXE_GENOUNZIP) command = PIZ;
        else if (exe_type == EXE_GENOCAT) { command = PIZ; flag.to_stdout = !flag.out_filename ; }
        else command = ZIP; // default 
    }

    unsigned num_files = argc - optind;

    flags_update (num_files, &argv[optind], is_short);

    // sort files by data type to improve VB re-using, and refhash-using files in the end to improve reference re-using
    qsort (&argv[optind], num_files, sizeof (argv[0]), main_sort_input_filenames);
    
    // determine how many threads we have - either as specified by the user, or by the number of cores
    if (flag.threads_str) {
        int ret = sscanf (flag.threads_str, "%u", &global_max_threads);
        ASSINP (ret == 1 && global_max_threads >= 1, "%s: %s requires an integer value of at least 1", global_cmd, OT("threads", "@"));
    }
    else global_max_threads = (double)arch_get_num_cores() * 1.4; // over-subscribe to keep all cores busy even when some threads are waiting on mutex or join
    
    // handle call commands except for ZIP, PIZ or LIST
    if (command == VERSION) { main_print_version();   return 0; }
    if (command == LICENSE) { license_display();      return 0; }
    if (command == HELP)    { main_print_help (true); return 0; }

    primary_command = command; 
    
    for (unsigned file_i=0; file_i < MAX (num_files, 1); file_i++) {

        char *next_input_file = optind < argc ? argv[optind++] : NULL;  // NULL means stdin
        
        if (next_input_file && !strcmp (next_input_file, "-")) next_input_file = NULL; // "-" is stdin too

        ASSINP (next_input_file || command != PIZ, 
                "%s: filename(s) required (redirecting from stdin is not possible)", global_cmd);

        ASSERTW (next_input_file || !flag.replace, "%s: ignoring %s option", global_cmd, OT("replace", "^")); 

        bool is_last_file = (file_i==num_files-1);

        main_load_reference (next_input_file, !file_i, is_last_file);
        
        switch (command) {
            case ZIP  : main_genozip (next_input_file, flag.out_filename, file_i==0, !next_input_file || is_last_file, argv[0]); break;
            case PIZ  : main_genounzip (next_input_file, flag.out_filename, is_last_file); break;           
            case LIST : main_genols (next_input_file, false, NULL, false); break;
            default   : ABORT ("%s: unrecognized command %c", global_cmd, command);
        }

        if (flag.pair) flag.pair = 3 - flag.pair; // alternate between PAIR_READ_1 and PAIR_READ_2
    }

    // if this is "list", finalize
    if (command == LIST) main_genols (NULL, true, NULL, false);

    return 0;
}
