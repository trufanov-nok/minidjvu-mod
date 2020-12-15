#include "AppOptions.h"
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include "../../src/base/mdjvucfg.h" /* for i18n */
#ifdef _WIN32
#include <fileapi.h>
#include <direct.h> // _chdir()
#else
#include <unistd.h> // chdir()
#endif
#include <string.h>
#include <ctype.h>
#include <iostream>

// allocate a copy of string src to dest

void copy_str_alloc(char** dest, const char* src)
{
    if (*dest) {
        free(*dest);
    }
    int len = strlen(src);
    *dest = (char*) calloc(sizeof(char), len+1);
    strncpy(*dest, src, len);
}

inline int is_tiff(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".tiff")
            || mdjvu_ends_with_ignore_case(path, ".tif");
}


/*********************************************************
 *
 *                     Image Options
 *
 *********************************************************/


struct ImageOptions* image_options_create(struct ImageOptions* defaults)
{
    struct ImageOptions* opt = (struct ImageOptions*) malloc(sizeof(struct ImageOptions));

    if (defaults) {
        opt->clean = defaults->clean;
        opt->smooth = defaults->smooth;
        opt->erosion = defaults->erosion;
        opt->dpi = defaults->dpi;
        opt->dpi_specified = defaults->dpi_specified;
    } else {
        opt->clean = opt->smooth = opt->erosion = 0;
        opt->dpi = 300;
        opt->dpi_specified = 0;
    }

    return opt;
}

void image_options_destroy(struct ImageOptions* opt)
{
    free(opt);
}

/*********************************************************
 *
 *                     Input File
 *
 *********************************************************/

struct InputFile* input_file_create()
{
    struct InputFile* in = (struct InputFile*) malloc(sizeof(struct InputFile));
    in->page = 0;
    in->name = NULL;
    in->image_options = NULL;

    in->output_dpi = 300;
    in->chunk_id = NULL;
    in->output_size = 0;
    in->djbz = NULL;
    return in;
}

void input_file_destroy(InputFile* in)
{
    if (in->name) {
        free(in->name);
    }
    if (in->image_options) {
        image_options_destroy(in->image_options);
    }
    if (in->chunk_id) {
        MDJVU_FREEV(in->chunk_id);
    }

    free(in);
}

/*********************************************************
 *
 *                     File List
 *
 *********************************************************/

void file_list_init(struct FileList* fl)
{
    fl->size = fl->data_size = 0;
    fl->files = NULL;
}

void file_list_reserve(struct FileList *fl, int size)
{
    if (size > fl->data_size) {
        struct InputFile** new_data = (struct InputFile**) calloc(sizeof(struct InputFile*), size);
        memcpy(new_data, fl->files, sizeof(char*) * fl->data_size);
        if (fl->files) {
            free(fl->files);
        }
        fl->files = new_data;
        fl->data_size = size;
    }
}

void file_list_clear(struct FileList* fl, int own_files)
{
    if (own_files) {
        for (int i = 0; i < fl->size; i++) {
            input_file_destroy(fl->files[i]);
        }
    }

    if (fl->files) {
        free(fl->files);
        fl->files = NULL;
    }

}

void file_list_add_file_ref(struct FileList *fl, struct InputFile* file)
{
    if (fl->data_size <= fl->size) {
        file_list_reserve(fl, fl->data_size + 10);
    }
    fl->files[fl->size++] = file;
}

void file_list_add_filename_and_page(struct FileList *fl, const char* fname, int page, struct ImageOptions* options)
{
    struct InputFile* in = input_file_create();
    copy_str_alloc(&in->name, fname);
    in->page = page;
    if (options) {
        in->image_options = image_options_create(options);
    }
    file_list_add_file_ref(fl, in);
}

void file_list_add_filename(struct FileList *fl, const char* fname, struct ImageOptions* options)
{

#ifdef HAVE_LIBTIFF
    int tiff_cnt = 0;
    if (is_tiff(fname) && (tiff_cnt = mdjvu_get_tiff_page_count(fname)) > 1)
    {
        // multipage tiff is detected
        for (int i = 0; i < tiff_cnt; i++) {
            file_list_add_filename_and_page(fl, fname, i, options);
        }
    } else
#endif
    {
        file_list_add_filename_and_page(fl, fname, 0, options);
    }
}


inline void adjust_pages_interval(const char* fname, int tiff_cnt, int& pg_min, int& pg_max)
{
    if (tiff_cnt < pg_max+1) {
        pg_max = tiff_cnt;
        std::cerr << "WARNING: file \"" << fname << "\" is expected to be multipage and pages " <<
                     pg_min << " to " << pg_max << " are requested. But it contains only " <<
                     tiff_cnt << " pages.\nMissing pages are ignored.\n";
    }
    if (pg_min+1 > tiff_cnt) {
        std::cerr << "ERROR: file \"" << fname << "\" is expected to be multipage and pages " <<
                     pg_min << " to " << pg_max << " are requested. But it contains only " <<
                     tiff_cnt << " pages.\nProcess aborted\n";
        exit(-2);
    }
}

void file_list_add_filename_with_filter(struct FileList *fl, const char* fname, int pg_min, int pg_max, struct ImageOptions* options)
{
#ifndef HAVE_LIBTIFF
    if (pg_max || pg_min) {
        std::cerr << "ERROR: file \"" << fname << "\" is expected to be multipage and pages " <<
                     pg_min << " to " << pg_max << " are requested. But minidjvu-mod is built" \
                     " with no libtiff support.\n";
        exit(-2);
    }
#else
    int tiff_cnt = 0;
    if (is_tiff(fname) && (tiff_cnt = mdjvu_get_tiff_page_count(fname)) > 1)
    {
        // multipage tiff is detected

        adjust_pages_interval(fname, tiff_cnt, pg_min, pg_max);

        for (int i = pg_min; i <= pg_max; i++) {
            file_list_add_filename_and_page(fl, fname, i, options);
        }
    } else
#endif
    {
        file_list_add_filename_and_page(fl, fname, 0, options);
    }
}


int file_list_find(struct FileList* fl, const struct InputFile* in)
{
    for (int i = 0; i < fl->size; i++) {
        struct InputFile* in2 = fl->files[i];
        if (in->page == in2->page && strcmp(in2->name, in->name) == 0) {
            return i;
        }
    }
    return -1;
}

int file_list_find_by_filename_and_page(struct FileList* fl, const char* fname, int page)
{
    for (int i = 0; i < fl->size; i++) {
        struct InputFile* in = fl->files[i];
        if (in->page == page && strcmp(in->name, fname) == 0) {
            return i;
        }
    }
    return -1;
}


void file_list_add_ref(struct FileList* dest, struct FileList* src, const char* fname)
{
    int idx = -1;
#ifdef HAVE_LIBTIFF
    int tiff_cnt = 0;
    if (is_tiff(fname) && (tiff_cnt = mdjvu_get_tiff_page_count(fname)) > 1)
    {
        // multipage tiff is detected
        for (int i = 0; i < tiff_cnt; i++) {
            idx = file_list_find_by_filename_and_page(src, fname, i);
            if (idx == -1) {
                std::cerr << "ERROR: page " << idx << " of file \"" << fname << "\" is referenced by djbz but can't be found.\n";
                exit(-2);
            }
            file_list_add_file_ref(dest, src->files[idx]);
        }
    } else
#endif
    {
        idx = file_list_find_by_filename_and_page(src, fname, 0);
        if (idx == -1) {
            std::cerr << "ERROR: file \"" << fname << "\" is referenced by djbz but can't be found.\n";
            exit(-2);
        }
        file_list_add_file_ref(dest, src->files[idx]);
    }
}

void file_list_add_ref_with_filter(struct FileList* dest, struct FileList* src, const char* fname, int pg_min, int pg_max)
{
    int idx = -1;
#ifndef HAVE_LIBTIFF
    if (pg_max || pg_min) {
        std::cerr << "ERROR: file \"" << fname << "\" is expected to be multipage and pages " <<
                     pg_min << " to " << pg_max << " are referenced by djbz. But minidjvu-mod is built" \
                     " with no libtiff support.\n";
        exit(-2);
    }
#else
    int tiff_cnt = 0;
    if (is_tiff(fname) && (tiff_cnt = mdjvu_get_tiff_page_count(fname)) > 1)
    {
        // multipage tiff is detected

        adjust_pages_interval(fname, tiff_cnt, pg_min, pg_max);

        for (int i = pg_min; i <= pg_max; i++) {
            idx = file_list_find_by_filename_and_page(src, fname, i);
            if (idx == -1) {
                std::cerr << "ERROR: page " << idx << " of file \"" << fname << "\" is referenced by djbz but can't be found.\n";
                exit(-2);
            }

            file_list_add_file_ref(dest, src->files[idx]);
        }
    } else
#endif
    {
        idx = file_list_find_by_filename_and_page(src, fname, 0);
        if (idx == -1) {
            std::cerr << "ERROR: file \"" << fname << "\" is referenced by djbz but can't be found.\n";
            exit(-2);
        }

        file_list_add_file_ref(dest, src->files[idx]);
    }
}




/*********************************************************
 *
 *                     Djbz Options
 *
 *********************************************************/

struct DjbzOptions* djbz_setting_create(struct DjbzOptions* defaults)
{
    struct DjbzOptions* djbz = (struct DjbzOptions*) calloc(sizeof(struct DjbzOptions), 1);

    djbz->id = NULL;
    djbz->dict_suffix = NULL;
    djbz->averaging = defaults ? defaults->averaging : 0;
    djbz->aggression = defaults ? defaults->aggression : 100;
    djbz->classifier = defaults ? defaults->classifier : 3;
    djbz->no_prototypes = defaults ? defaults->no_prototypes : 0;
    djbz->erosion = defaults ? defaults->erosion : 0;
    if (defaults) {
        copy_str_alloc(&djbz->dict_suffix, defaults->dict_suffix);
    } else {
        copy_str_alloc(&djbz->dict_suffix, "djbz\0");
    }

    file_list_init(&djbz->file_list_ref);

    djbz->chunk_id = NULL;
    djbz->do_not_save = 0;
    djbz->output_size = 0;
    return djbz;
}

void djbz_setting_destroy(struct DjbzOptions* djbz)
{
    if (djbz->id) {
        free(djbz->id);
    }
    if (djbz->dict_suffix) {
        free(djbz->dict_suffix);
    }
    if (djbz->chunk_id) {
        MDJVU_FREEV(djbz->chunk_id);
    }
    // djbz don't own files in its filelist
    file_list_clear(&djbz->file_list_ref, 0);

    free(djbz);
}

/*********************************************************
 *
 *                     Djbz List
 *
 *********************************************************/

void djbz_list_init(struct DjbzList* dl)
{
    dl->size = dl->data_size = 0;
    dl->djbzs = NULL;
}

void djbz_list_reserve(struct DjbzList *dl, int size)
{
    if (size > dl->data_size) {
        struct DjbzOptions** new_data = (struct DjbzOptions**) calloc(sizeof(struct DjbzOptions*), size);
        memcpy(new_data, dl->djbzs, sizeof(char*) * dl->data_size);
        if (dl->djbzs) {
            free(dl->djbzs);
        }
        dl->djbzs = new_data;
        dl->data_size = size;
    }
}

void djbz_list_clear(struct DjbzList* dl)
{
    for (int i = 0; i < dl->size; i++) {
        djbz_setting_destroy(dl->djbzs[i]);
    }

    if (dl->djbzs) {
        free(dl->djbzs);
        dl->djbzs = NULL;
    }
}

void djbz_list_add_option(struct DjbzList* dl, struct DjbzOptions* djbz)
{
    const int size = dl->size;
    struct DjbzOptions** new_data = (struct DjbzOptions**) calloc(sizeof(struct DjbzOptions*), size + 1);
    memcpy(new_data, dl->djbzs, sizeof(struct DjbzOptions*) * size);
    if (dl->djbzs) {
        free(dl->djbzs);
    }
    dl->djbzs = new_data;
    dl->djbzs[dl->size++] = djbz;
}


/*********************************************************
 *
 *                  Chunk file operations
 *
 *********************************************************/

void chunk_file_create(struct ChunkFile* cf, char* fname)
{
    cf->file = NULL;
    cf->indirect_mode = fname ? 1: 0;
    cf->filename = fname;

    if (!cf->indirect_mode) {
#ifdef _WIN32
        cf->chunk_filename = (char*) MDJVU_MALLOCV(char, MAX_PATH+1);
        char pathname[MAX_PATH+1];
        if (GetTempPathA(MAX_PATH+1, pathname) == 0) {
            fprintf(stderr, _("Could not create a temporary file. (GetTempPathA)\n"));
            exit(-2);
        }
        char prefix[4] = {'d', 'j', 'v', 0 };
        if (GetTempFileNameA(pathname, prefix, 0, cf->filename) == 0) {
            fprintf(stderr, _("Could not create a temporary file. (GetTempFileNameA)\n"));
            exit(-2);
        }
#else
        cf->file = tmpfile(); // opens it too
        if (!cf->file) {
            fprintf(stderr, _("Unique filename cannot be generated or the unique file cannot be opened (errno: %s)\n"), strerror(errno));
            exit(-2);
        }
#endif
    }
}

void chunk_file_open(struct ChunkFile* cf)
{
    if (cf->indirect_mode) {
        cf->file = fopen(cf->filename, "w+b");
    } else {
#ifdef _WIN32
        cf->file = fopen(cf->chunk_filename, "w+b");
#else
        // already opened
#endif
    }

    if(!cf->file) {
        fprintf(stderr, _("Can't open file (errno: %s)\n"), strerror(errno));
        exit(-2);
    }
}

void chunk_file_close(struct ChunkFile* cf)
{
    int res = 0;
    if (cf->indirect_mode) {
        res = fclose(cf->file);
        cf->file = NULL;
    } else {
#ifdef _WIN32
        res = fclose(cf->file);
        cf->file = NULL;
#else
        // shouldn't be closed!
#endif
    }
    if (res != 0) {
        fprintf(stderr, _("Can't close file (errno: %s)\n"), strerror(errno));
        exit(-2);
    }
}

void chunk_file_destroy(struct ChunkFile* cf)
{
#ifdef _WIN32
    if (remove(cf->chunk_filename) != 0) {
        fprintf(stderr, _("Can't remove file (errno: %s)\n"), strerror(errno));
        exit(-2);
    }
    MDJVU_FREEV(cf->chunk_filename);
#else
    if (!cf->indirect_mode) {
        // not indirect
        if (fclose(cf->file) != 0) {
            fprintf(stderr, _("Can't close file (errno: %s)\n"), strerror(errno));
            exit(-2);
        }
#endif
    }

    cf->file = NULL;
}

/*********************************************************
 *
 *                  file name template routines
 *
 *********************************************************/

static int get_ext_delim_pos(const char *fname)
{
    size_t pos = strcspn(fname,".");
    size_t last = 0;

    while (last + pos != strlen(fname))
    {
        last += (pos + 1);
        pos = strcspn(fname + last,".");
    }
    return last;
}

static char *get_page_or_dict_name(const char *fname, int page)
{
    const int extpos = get_ext_delim_pos(fname);
    char * page_name = MDJVU_MALLOCV(char, extpos + 10);
    memset(page_name, 0, sizeof(char) *(extpos + 10));
    if (extpos > 0) {
        strncpy(page_name, fname, extpos-1);
    }
    if (page > 0) {
        sprintf(page_name + (extpos - 1),"#%03d.djvu", page+1);
    } else {
        strcat(page_name, ".djvu");
    }

    return(page_name);
}

static void replace_suffix(char *name, const char *suffix)
{
    int len = strlen(name);

    name[len-4] = '\0';
    strcat(name, suffix);
}

static const char *strip(const char *str, char sep)
{
    const char *t = strrchr(str, sep);
    if (t)
        return t + 1;
    else
        return str;
}


/* return path without a directory name */
static const char *strip_dir(const char *path)
{
    return path ? strip(strip(path, '\\'), '/') : NULL;
}



/*********************************************************
 *
 *                  Application options
 *
 *********************************************************/



void app_options_init(struct AppOptions* opts)
{
    /* options */
    opts->pages_per_dict = 10; /* 0 means infinity */
    opts->verbose = 0;
    opts->match = 0;
    opts->Match = 0;
    opts->report = 0;
    opts->warnings = 0;
    opts->indirect = 0;

#ifdef _OPENMP
    opts->max_threads = 0;
#endif

    opts->default_image_options = image_options_create(NULL);
    opts->default_djbz_options = djbz_setting_create(NULL);
    file_list_init(&opts->file_list);
    djbz_list_init(&opts->djbz_list);

    opts->output_file = NULL;
}

void app_options_free(struct AppOptions* opts)
{
    if (opts->output_file) {
        free(opts->output_file);
    }

    djbz_list_clear(&opts->djbz_list);
    file_list_clear(&opts->file_list, 1);


    free(opts->default_image_options);

    // don't free opts itself as we expect pointer to static object
    //    free(opts);
}

void app_options_set_djbz_suffix(struct AppOptions* opts, char* suffix)
{
    if (opts) {
        if (suffix) {
            copy_str_alloc(&opts->default_djbz_options->dict_suffix, suffix);
        } else {
            free(opts->default_djbz_options->dict_suffix);
            opts->default_djbz_options->dict_suffix = NULL;
        }
    }
}

void app_options_set_output_file(struct AppOptions* opts, const char* filename)
{
    int len = strlen(filename);
    if (filename && len > 0) {
        const char* file = filename;
        int idx = len-1;
        do {
            char c = filename[idx];
            if (c == '/' || c == '\\') {
                break;
            }
        } while (--idx >= 0);

        if (idx != -1) {
            file = filename+idx+1;
            char* path = NULL;
            copy_str_alloc(&path, filename);
            path[idx] = '\0';
#ifndef _WIN32
            if (chdir(path)==-1) {
#else
            if (_chdir(path)==-1) {
#endif
                std::cerr << "ERROR: can't change working directory to: " << path << "\n";
                exit(-2);
            }
            MDJVU_FREEV(path);
            copy_str_alloc(&opts->output_file, file);
        }

    }
}


void app_options_autocomplete_djbzs(struct AppOptions* opts)
{
    struct DjbzOptions * new_djbz = djbz_setting_create(opts->default_djbz_options);
    for(int i = 0; i < opts->file_list.size; i++) {
        struct InputFile* in = opts->file_list.files[i];
        if (in->djbz == NULL) {
            file_list_add_file_ref(&new_djbz->file_list_ref, in);
            in->djbz = new_djbz;

            if (new_djbz->file_list_ref.size >= opts->pages_per_dict) {
                djbz_list_add_option(&opts->djbz_list, new_djbz);
                new_djbz = djbz_setting_create(opts->default_djbz_options);
            }
        }
    }
    if (new_djbz->file_list_ref.size > 0) {
        djbz_list_add_option(&opts->djbz_list, new_djbz);
    } else {
        djbz_setting_destroy(new_djbz);
    }
}

void app_options_construct_chunk_ids(struct AppOptions* opts)
{
    // initialize elements with filenames as it can't be done in parallel blocks without critical sections
    for (int i = 0; i < opts->file_list.size; i++)
    {
        struct InputFile* in = opts->file_list.files[i];
        if (in->chunk_id) {
            MDJVU_FREEV(in->chunk_id);
        }
        in->chunk_id = get_page_or_dict_name(strip_dir(in->name), in->page);
        chunk_file_create(&in->chunk_file, opts->indirect ? in->chunk_id : NULL);
    }

    for (int i = 0; i < opts->djbz_list.size; i++) {
        struct DjbzOptions* djbz = opts->djbz_list.djbzs[i];
        if (djbz->chunk_id) {
            MDJVU_FREEV(djbz->chunk_id);
        }
        const char* suffix = djbz->dict_suffix ? djbz->dict_suffix : opts->default_djbz_options->dict_suffix;
        if (djbz->id) {
            const int len_id = strlen(djbz->id);
            const int len_suffix = strlen(suffix);
            djbz->chunk_id = MDJVU_MALLOCV(char, len_id + len_suffix + 2);
            sprintf(djbz->chunk_id, "%s.%s", djbz->id, suffix);
        } else {
            const char * path = djbz->file_list_ref.files[0]->chunk_id;
            const int size = strlen(path) + strlen(suffix) - 2;
            djbz->chunk_id = MDJVU_MALLOCV(char, size);
            strncpy(djbz->chunk_id, path, size);
            replace_suffix(djbz->chunk_id, suffix);
        }

        chunk_file_create(&djbz->chunk_file, opts->indirect ? djbz->chunk_id : NULL);
    }
}

