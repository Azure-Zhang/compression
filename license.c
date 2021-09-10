// ------------------------------------------------------------------
//   license.c
//   Copyright (C) 2020-2021 Black Paw Ventures Limited
//   Please see terms and conditions in the file LICENSE.txt

#ifdef _WIN32
#include <direct.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "genozip.h"
#include "strings.h"
#include "arch.h"
#include "profiler.h" // for TimeSpecType
#include "md5.h"
#include "url.h"
#include "license.h"
#include "version.h"
#include "buffer.h"
#include "flags.h"
#include "md5.h"
#include "website.h"
#include "version.h"
#include "file.h"

// these field names appear in the license file starting V12.0.7
#define LIC_FIELD_VERSION     "Genozip license version"
#define LIC_FIELD_INSTITUTION "License granted to"
#define LIC_FIELD_NAME        "Accepted by (name)"
#define LIC_FIELD_EMAIL       "Accepted by (email)"
#define LIC_FIELD_TIMESTAMP   "Timestamp of acceptance"
#define LIC_FIELD_IP          "IP address of acceptance"
#define LIC_FIELD_NUMBER      "License number"

#include "text_license.h"

static const char *license_filename = NULL;  // non-standard filename set with --licfile

static struct {
    bool initialized;
    bool has_details; // false if old-style license

    char name[256], email[256], institution[1024], ip[ARCH_IP_LEN], version[20];
    StrText timestamp;
    uint32_t license_num;
} rec = {};

static uint32_t license_calc_number (const Buffer *license_data)
{
    char data_no_ws[license_data->len];
    unsigned data_no_ws_len = str_remove_whitespace (license_data->data, license_data->len, data_no_ws);        

    return md5_do (data_no_ws, data_no_ws_len).words[0];
}

static void license_generate (Buffer *license_data)
{
    for (unsigned i=0; i < sizeof(license) / sizeof(char*); i++) {
        buf_add_string (evb, license_data, license[i]); // allocs one extra char
        NEXTENT (char, *license_data) = '\n';
    }

    bufprintf (evb, license_data, 
               LIC_FIELD_INSTITUTION": %s\n"
               LIC_FIELD_NAME": %s\n"
               LIC_FIELD_EMAIL": %s\n"
               LIC_FIELD_TIMESTAMP": %s\n"
               LIC_FIELD_IP": %s\n",
               rec.institution, rec.name, rec.email, rec.timestamp.s, rec.ip);
    
    rec.initialized = true;
    rec.has_details = true;
    rec.license_num = license_calc_number (license_data);
    strcpy (rec.version, GENOZIP_CODE_VERSION);

    bufprintf (evb, license_data, LIC_FIELD_NUMBER": %u\n", rec.license_num);
}

void license_set_filename (const char *filename)
{
    struct stat sb;
    ASSINP (!stat (filename, &sb), "Failed to access license file %s: %s", filename, strerror (errno));

    license_filename = filename;
}

static const char *get_license_filename (bool create_folder_if_needed)
{
    if (license_filename) return license_filename; // non-standard filename set with --licfile

#ifdef _WIN32
    ASSINP0 (getenv ("APPDATA"), "%s: cannot store license, because APPDATA env var is not defined");

    char folder[500];
    sprintf (folder, "%s/genozip", getenv ("APPDATA"));

    if (create_folder_if_needed) {
        int ret = _mkdir (folder); 
        ASSERT (ret >= 0 || errno == EEXIST, "failed to created the folder %s", folder);
    }

#else
    const char *folder = getenv ("HOME");
    ASSINP0 (folder, "%s: cannot calculate license file name, because $HOME env var is not defined");
#endif    

    char *filename = MALLOC (strlen(folder) + 50);
    sprintf (filename, "%s/.genozip_license", folder);

    return filename;
}

static const char *license_load_field (const char *field, unsigned n_lines, const char **lines, const unsigned *line_lens)
{
    unsigned field_len = strlen (field);

    for (int i=n_lines-1; i >= 0; i--)
        if (line_lens[i] > field_len+2 && !memcmp (lines[i], field, field_len) && lines[i][field_len] == ':' && lines[i][field_len+1] == ' ')
            return &lines[i][field_len+2];

    return ""; // not found
}

// IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
static void license_load (void)
{
    if (rec.initialized) return;

    const char *filename = get_license_filename (true);
    
    if (!file_exists (filename)) {
        flag.do_register = "";
        license_register ();
        return;
    }

    file_split_lines (filename, "license");
    
    // case: license generated by Genozip 12.0.6 or older
    if (n_lines == 1) {
        ASSINP (line_lens[0] >= 1 && IS_DIGIT (lines[0][0]), "failed to parse license file %s, please re-register with genozip --register", filename);
        rec.license_num = atoi (lines[0]);
    }

    // case: full license file
    else {
        char license_num_str[30] ="";
        #define COPY_FIELD(var,field) strncpy (var, license_load_field (field, n_lines, lines, line_lens), sizeof (var)-1)
        COPY_FIELD (rec.version, LIC_FIELD_VERSION);
        COPY_FIELD (rec.institution, LIC_FIELD_INSTITUTION);
        COPY_FIELD (rec.name, LIC_FIELD_NAME);
        COPY_FIELD (rec.email, LIC_FIELD_EMAIL);
        COPY_FIELD (rec.timestamp.s, LIC_FIELD_TIMESTAMP);
        COPY_FIELD (rec.ip, LIC_FIELD_IP);
        COPY_FIELD (license_num_str, LIC_FIELD_NUMBER);

        ASSINP (strlen (license_num_str) >= 1 && IS_DIGIT (license_num_str[0]), "failed to parse license file %s, please re-register with genozip --register", filename);
        rec.license_num = atoi (license_num_str);
        rec.has_details = true;

        data.len -= line_lens[n_lines-1] + 2;
        ASSINP0 (rec.license_num == license_calc_number (&data), "Invalid license. Please re-register with genozip --register");
    }

    rec.initialized = true;

    buf_destroy (&data);
}

static bool license_submit (char commerical, char update, const char *os, unsigned cores, const char *endianity, const char *user_host, const char *dist)
{
    // reference: https://stackoverflow.com/questions/18073971/http-post-to-a-google-form/47444396#47444396

    // FORM_ID is in the url when you preview your form
    #define PREFIX "https://docs.google.com/forms/d/e/1FAIpQLSc6pSBIOBsS5Pu-JNvfnLWV2Z1W7k-4f2pKRo5rTbiU_sCnIw/formResponse"
    
    /* To get entry IDs - in Chrome browser: 1. open form 2. click on eye icon to Preview 2. right-click Inspect 3. go to "console" tab 4. run this code:
    function loop(e){
    if(e.children)
    for(let i=0;i<e.children.length;i++){
        let c = e.children[i], n = c.getAttribute('name');
        if(n) console.log(`${c.getAttribute('aria-label')}: ${n}`);
        loop(e.children[i]);
     }
    }; loop(document.body);
    */

    // note: identical to register.sh
    char *url_format = PREFIX
                       "?entry.344252538=%.100s"
                       "&entry.926671216=%.50s"
                       "&entry.1734045469=%.50s"
                       "&entry.2009586582=%c"
                       "&entry.119966790=%c"
                       "&entry.81542373=%s"
                       "&entry.1668073218=%u"
                       "&entry.1943454647=%s"
                       "&entry.1763961212=%s"
                       "&entry.1655649315=%u"
                       "&entry.186159495=%s"
                       "&entry.1598028195=%s"
                       "&entry.1384715202=%s";

    char *institutionE = url_esc_non_valid_chars (rec.institution);
    char *nameE        = url_esc_non_valid_chars (rec.name);
    char *emailE       = url_esc_non_valid_chars (rec.email);
    char *osE          = url_esc_non_valid_chars (os);
    char *user_hostE   = url_esc_non_valid_chars (user_host);

    char url[sizeof (rec) + 200];
    sprintf (url, url_format, institutionE, nameE, emailE, commerical, update, osE, cores, rec.ip, user_hostE, rec.license_num, rec.version, dist, endianity);

    bool success = url_read_string (url, NULL, 0) >= 0;
    
    FREE (institutionE); FREE (nameE); FREE (emailE); FREE (osE); FREE (user_hostE);
    return success;
}

static bool license_verify_email (char *response, unsigned response_size, const char *unused)
{
    // sanity check that this is an email address
    return strlen (response) > 3 && strchr (response, '@') && strchr (response, '.');
}

static bool license_verify_name (char *response, unsigned response_size, const char *unused)
{
    if (!strchr (response, ' ')) {
        fprintf (stderr, "Please enter your full name\n");
        return false;
    }
    
    return true;
}

static void license_exit_if_not_confirmed (const char *response)
{
    if (response[0] == 'N') {
        fprintf (stderr, "\nYou have not registered. You may register at any time in the future.\n\nWishing you a wonderful day from the Genozip team! https://genozip.com\n");
        exit_ok();
    }
}

// UI flow to generate a new license for the user

// IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
void license_register (void)
{
    char confirm[100], commercial[100], update[100];
    const char *os, *dist, *endianity, *user_host;
    unsigned cores;

    str_split (flag.do_register, strlen (flag.do_register), 11, '|', field, true);
    str_nul_separate (field);

    // if stdin or stderr is redirected - we cannot ask the user an interactive question
    ASSINP0 (isatty(0) && isatty(2), "Use of genozip is free for non-commercial purposes, but requires registration. Please run: genozip --register.\n"
                                     "If you are unable to register (for example because this is a batch-job machine) please see: " WEBSITE_USING_ON_HPC);


    const char *filename = get_license_filename (true);

    if (!n_fields) {

        fprintf (stderr, "Welcome to genozip!\n\n"
                         "The use of genozip for non-commercial purposes (as defined in the license: "WEBSITE_LICENSE") is FREE, but requires registration.\n"
                         "If you are not sure whether your usage is considered non-commercial, please email "EMAIL_REGISTER"\n\n");

        if (file_exists (filename)) 
            str_query_user ("You are already registered. Are you sure you want to register again? (y or n) ", confirm, sizeof(confirm), str_verify_y_n, NULL);
        else 
            str_query_user ("Would you like to register now? ([y] or n) ", confirm, sizeof(confirm), str_verify_y_n, "Y");

        license_exit_if_not_confirmed (confirm);
    }

    file_remove (filename, true); // remove old license, if one exists

    if (n_fields) {
        strncpy (rec.institution, fields[0], sizeof(rec.institution)-1);
        strncpy (rec.name,        fields[1], sizeof(rec.name)-1);
        strncpy (rec.email,       fields[2], sizeof(rec.email)-1);
        strncpy (commercial,      fields[3], sizeof(commercial)-1);
        strncpy (update,          fields[4], sizeof(update)-1);
        strncpy (rec.ip,          fields[5], sizeof(rec.ip)-1);
        os        = fields[6];
        dist      = fields[7]; 
        endianity = fields[8];
        user_host = fields[9];
        cores     = atoi(fields[10]);
    }
    else {
        fprintf (stderr, "\nLicense details -\n");
    
        str_query_user ("\nInstitution / Company name: ", rec.institution, sizeof(rec.institution), str_verify_not_empty, NULL);

        str_query_user ("\nYour name: ", rec.name, sizeof(rec.name), license_verify_name, NULL);
        
        str_query_user ("\nYour email address: ", rec.email, sizeof(rec.email), license_verify_email, NULL);
        
        str_query_user ("\nDo you require a commercial license? If yes, we will contact you (this will not stop the registration now) (y or [n]) ", 
                        commercial, sizeof(commercial), str_verify_y_n, "N");

        str_query_user ("\nShall we update you by email when new features are added to genozip? ([y] or n) ", 
                        update, sizeof(update), str_verify_y_n, "Y");

        fprintf (stderr, "\n\nPlease read the terms and conditions of the license:\n\n"); 
        license_display(); 
        fprintf (stderr, "\n"); 

        str_query_user ("Do you accept the terms and conditions of the license? (y or n) ", confirm, sizeof(confirm), str_verify_y_n, NULL);
        license_exit_if_not_confirmed (confirm);

        os        = arch_get_os();
        dist      = DISTRIBUTION;
        cores     = arch_get_num_cores();
        endianity = arch_get_endianity();
        user_host = arch_get_user_host();
        memcpy (rec.ip, arch_get_ip_addr ("Failed to register the license"), ARCH_IP_LEN);
    }

    rec.timestamp = str_time();

    static Buffer license_data = EMPTY_BUFFER;
    license_generate (&license_data);

    if (!n_fields) {

        fprintf (stderr, "\nThank you. To complete your license registration, genozip will now submit the following information to the genozip licensing server:\n\n");

        // note: text needs to match scripts/register.sh
        fprintf (stderr, "=====================================================================\n");
        fprintf (stderr, LIC_FIELD_INSTITUTION": %s\n", rec.institution);
        fprintf (stderr, LIC_FIELD_NAME": %s\n", rec.name);
        fprintf (stderr, LIC_FIELD_EMAIL": %s\n", rec.email);
        fprintf (stderr, "Commercial: %s\n", commercial[0]=='Y' ? "Yes" : "No");
        fprintf (stderr, "Send new feature updates: %s\n", update[0]=='Y' ? "Yes" : "No");
        fprintf (stderr, "System info: OS=%s cores=%u endianity=%s IP=%s\n", os, cores, endianity, rec.ip);
        fprintf (stderr, "Username: %s\n", user_host);
        fprintf (stderr, "Genozip info: version=%s distribution=%s\n", GENOZIP_CODE_VERSION, dist);
        fprintf (stderr, "Genozip license number: %u\n", rec.license_num);
        fprintf (stderr, "I accept the terms and conditions of the Genozip license\n");
        fprintf (stderr, "=====================================================================\n\n");
        
        str_query_user ("Proceed with completing the registration? ([y] or n) ", confirm, sizeof(confirm), str_verify_y_n, "Y");
        license_exit_if_not_confirmed (confirm);
    }
        
    bool submitted = license_submit (commercial[0], update[0], os, cores, endianity, user_host, dist);

    ASSINP0 (submitted,
             "Failed to register the license, possibly because the Internet is not accessible or the registration server\n"
             "(which is hosted on a Google server) is not accessible. If this problem persists, you can register manually by\n"
             "sending an email to register@genozip.com - copy & paste the lines between the \"======\" into the email message.\n");

    ASSINP (file_put_data (filename, license_data.data, license_data.len), 
            "Failed to write license file %s: %s. If this is unexpected, email "EMAIL_SUPPORT" for help.", filename, strerror (errno));

    if (!n_fields) 
        fprintf (stderr, "\nSUCCESS. A Genozip license has been granted:\n"
                         "License type: Single User\nLicensee: %s\nFor use by: %s\n\n" 
                         "Documentation: " GENOZIP_URL "\n\n"
                         "Support: " EMAIL_SUPPORT "\n\n"
                         "Citing: " WEBSITE_PUBLICATIONS "\n\n", rec.institution, rec.name);

    buf_destroy (&license_data);
}

// IF YOU'RE CONSIDERING EDITING THIS CODE TO BYPASS THE REGISTRTION, DON'T! It would be a violation of the license,
// and might put you personally as well as your organization at legal and financial risk - see "Unauthorized use of Genozip"
// section of the license. Rather, please contact sales@genozip.com to discuss which license would be appropriate for your case.
uint32_t license_get_number (void)
{
    license_load();

    return rec.license_num;
}

const char *license_get_one_line (void)
{
    static char s[sizeof (rec) + sizeof (rec.name) + 200];

    sprintf (s, "License v%s granted to: %s for use by: %s accepted by: %s <%s> on %s from IP=%s", 
             rec.version, rec.institution, rec.name, rec.name, rec.email, rec.timestamp.s, rec.ip);

    return s;
}

bool license_has_details (void)
{
    return rec.has_details;
}

void license_display (void)
{
    const char *filename = get_license_filename (false);
    static Buffer license_data = {};
    
    if (file_exists (filename) && !flag.force) 
        file_get_file (evb, filename, &license_data, "license_data", true);

    // case: user has already accepted the license and it is new style license - display the license file
    if (license_data.len > 100) {
        str_split (license_data.data, license_data.len, 0, '\n', line, false);
        str_nul_separate (line);
        str_print_text (lines, n_lines-1, "", "\n\n", flag.lic_width);
    }
    
    // case: license not yet accepted or old style (up to 12.0.6) license - display the current version license
    else
        str_print_text (license, sizeof(license) / sizeof(char*), "", "\n\n", flag.lic_width);  // Makefile sets lic_width to a fixed width for Windows Installer and for Docs
}

