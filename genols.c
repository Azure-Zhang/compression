// ------------------------------------------------------------------
//   genols.c
//   Copyright (C) 2019-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#include <fcntl.h>
#include <dirent.h>
#include <errno.h>
#include <sys/types.h>

#include "genozip.h"
#include "genols.h"
#include "file.h"
#include "buffer.h"
#include "flags.h"
#include "strings.h"
#include "zfile.h"
#include "vblock.h"
#include "endianness.h"

static void genols_list_dir(); // forward declaration

void genols (const char *z_filename, bool finalize, const char *subdir, bool recursive) 
{
    if (!finalize) {
        // no specific filename = show entire directory
        if (!z_filename) {
            genols_list_dir("."); 
            return;
        }

        char *last_c = (char *)&z_filename[strlen(z_filename)-1];
        if (*last_c == '/' 
#ifdef _WIN32
            || *last_c == '\\'
#endif
            ) {
            *last_c = 0; // remove trailing '/' or '\'
            genols_list_dir (z_filename);
            return;
        }
 
        // filename is a directory - show directory contents (but not recursively)
        if (!subdir && file_is_dir (z_filename)) {
            genols_list_dir (z_filename);
            return;
        }
    }

    static bool first_file = true;
    static unsigned files_listed=0, files_ignored=0;
    static int64_t total_uncompressed_len=0, total_compressed_len=0;

    const unsigned FILENAME_WIDTH = 40;

    const char *head_format = "\n%-7.7s %11s %10s %10s %6s %s  %*s %s\n";
    const char *foot_format = "\nTotal: %3u files    %10s %10s %5.*fX\n";
    const char *item_format = "%-7.7s %11s %10s %10s %5.*fX  %s  %s%s%*.*s %s\n";

    const char *head_format_bytes = "\n%-7.7s %11s %15s %15s %6s  %*s\n";
    const char *foot_format_bytes = "\nTotal: %3u files    %15s %15s %5.*fX\n";
    const char *item_format_bytes = "%-7.7s %11s %15s %15s %5.*fX  %s%s%*.*s\n";

    // we accumulate the string in str_buf and print in the end - so it doesn't get mixed up with 
    // warning messages regarding individual files
    static Buffer str_buf = EMPTY_BUFFER; 
    
    if (finalize) {
        if (files_listed > 1) {
            double ratio = total_compressed_len ? ((double)total_uncompressed_len / (double)total_compressed_len) : 0;

            if (flag.bytes) 
                bufprintf (evb, &str_buf, foot_format_bytes, files_listed, str_int_s (total_compressed_len).s, 
                           str_int_s (total_uncompressed_len).s, ratio < 100, ratio);
            else 
                bufprintf (evb, &str_buf, foot_format, files_listed, str_size(total_compressed_len).s, 
                           str_size(total_uncompressed_len).s, ratio < 100, ratio);
        }
        
        ASSERTW (!files_ignored, "Ignored %u file%s that %s not have a " GENOZIP_EXT " extension, or are invalid.", 
                 files_ignored, files_ignored==1 ? "" : "s", files_ignored==1 ? "does" : "do");
        
        goto finish;
    }

    if (first_file) {
        if (flag.bytes) 
            bufprintf (evb, &str_buf, head_format_bytes, "Type", "Lines", "Compressed", "Original", "Factor", -(int)FILENAME_WIDTH, "Name");
        else
            bufprintf (evb, &str_buf, head_format, "Type", "Lines", "Compressed", "Original", "Factor", " MD5 of original textual file    ", -(int)FILENAME_WIDTH, "Name", "Creation");
        
        first_file = false;
    }
    
    if (!file_has_ext (z_filename, GENOZIP_EXT) || !file_exists (z_filename)) {
        files_ignored++;
        goto finish;
    }

    bool is_subdir = subdir && (subdir[0] != '.' || subdir[1] != '\0');

    z_file = file_open (z_filename, READ, Z_FILE, 0); // open global z_file

    SectionHeaderGenozipHeader header;
    if (!zfile_read_genozip_header (&header))
        goto finish;

    double ratio = z_file->disk_size ? ((double)z_file->txt_data_so_far_bind / (double)z_file->disk_size) : 0;
    
    // TODO: have an option to print ref_file_name and ref_file_md5

    DataType dt = z_file->z_flags.txt_is_bin ? DTPZ (bin_type) : z_file->data_type;

    if (flag.bytes) 
        bufprintf (evb, &str_buf, item_format_bytes, dt_name (dt), str_uint_commas (z_file->num_lines).s, 
                   str_int_s (z_file->disk_size).s, str_int_s (z_file->txt_data_so_far_bind).s, ratio < 100, ratio, 
                   (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
                   is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH, TXT_FILENAME_LEN,
                   z_filename);
    else 
        bufprintf (evb, &str_buf, item_format, dt_name (dt), str_uint_commas (z_file->num_lines).s,
                   str_size (z_file->disk_size).s, str_size (z_file->txt_data_so_far_bind).s, ratio < 100, ratio, 
                   digest_display_ex (z_file->digest, DD_MD5_IF_MD5).s,
                   (is_subdir ? subdir : ""), (is_subdir ? "/" : ""),
                   is_subdir ? -MAX (1, FILENAME_WIDTH - 1 - strlen(subdir)) : -FILENAME_WIDTH, TXT_FILENAME_LEN,
                   z_filename, header.created);

    total_compressed_len   += z_file->disk_size;
    total_uncompressed_len += z_file->txt_data_so_far_bind;
    
    files_listed++;

    // if --list, OR if the user did genols on one file (not a directory), show bound components, if there are any
    if (flag.list || (!flag.multiple_files && !recursive && z_file->num_components >= 2)) {
        buf_add_string (evb, &str_buf, "Components:\n");
        Section sl_ent = NULL;
        uint64_t num_lines_count=0;
        while (sections_next_sec (&sl_ent, SEC_TXT_HEADER)) {
            SectionHeaderTxtHeader header = zfile_read_section_header (evb, sl_ent->offset, sl_ent->vblock_i, SEC_TXT_HEADER).txt_header;

            num_lines_count += BGEN64 (header.txt_num_lines);
            bufprintf (evb, &str_buf, item_format, "", str_uint_commas (BGEN64 (header.txt_num_lines)).s, "", 
                       str_size (BGEN64 (header.txt_data_size)).s, 
                       0, 0.0, digest_display_ex (header.digest_single, DD_MD5_IF_MD5).s, "", "",
                       -(int)FILENAME_WIDTH, TXT_FILENAME_LEN, header.txt_filename, "");
        }

        if (num_lines_count != z_file->num_lines)
            buf_add_string (evb, &str_buf, "\nNote: the difference between the file's num_lines and the total of its components' is the number of lines of the 1st component's header\n");
    }
    file_close (&z_file, false, false);

finish:
    if (!recursive) {
        buf_print (&str_buf, false);
        buf_free (&str_buf);
    }
}

static void genols_list_dir(const char *dirname)
{
    DIR *dir;
    struct dirent *ent;

    dir = opendir (dirname);
    ASSERT (dir, "failed to open directory: %s", strerror (errno));

    int ret = chdir (dirname);
    ASSERT (!ret, "failed to chdir(%s)", dirname);

    flag.multiple_files = true;

    while ((ent = readdir(dir))) 
        if (!file_is_dir (ent->d_name))  // don't go down subdirectories recursively
            genols (ent->d_name, false, dirname, true);
    
    closedir(dir);    

    ret = chdir ("..");
    ASSERT0 (!ret, "failed to chdir(..)");
}
