#ifndef APPOPTIONS_H
#define APPOPTIONS_H

#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h> // for NULL
#include <stdio.h>  // for FILE


struct ChunkFile
{
    FILE * file;
    char * filename;
    int indirect_mode;
};

struct ImageOptions
{
    int dpi_specified;
    int dpi;
    int smooth;
    int clean;
    int erosion;
};

struct DjbzOptions;

struct InputFile
{
    char* name;
    int page;
    // if null, use AppOptions::default_image_options
    struct ImageOptions* image_options;

    // <private>

    // output_dpi is the dpi that used for saving result image.
    // it's 300 by default. And may be overriden by default image options
    // or image options if they have dpi_specified flag.
    // if they aren't specified the encoder will try to determine dpi
    // from x resolution of tiff file and use it.
    int output_dpi;
    char * chunk_id;
    int output_size;
    struct ChunkFile chunk_file;
    struct DjbzOptions* djbz;
};

struct FileList
{
    int size;
    int data_size;
    struct InputFile** files;
};

struct DjbzOptions
{
    char * id;
    char* dict_suffix;
    int averaging;
    int aggression;
    int classifier;
    int no_prototypes;
    int erosion;
    struct FileList file_list_ref;
    // <private>
    char * chunk_id;
    int do_not_save;
    int output_size;
    struct ChunkFile chunk_file;
};

struct DjbzList
{
    int size;
    int data_size;
    struct DjbzOptions** djbzs;
};

struct AppOptions
{
    int pages_per_dict; /* 0 means infinity */
    int verbose;
    int match;
    int Match;
    int report;
    int warnings;
    int indirect;

    #ifdef _OPENMP
    int max_threads;
    #endif

    struct ImageOptions* default_image_options;
    struct FileList file_list;

    struct DjbzOptions* default_djbz_options;
    struct DjbzList djbz_list;

    char* output_file;
};

#ifdef __cplusplus
extern "C" {
#endif

void copy_str_alloc(char** dest, const char *src);

struct ImageOptions* image_options_create(struct ImageOptions* defaults);

void file_list_init(struct FileList* fl);
void file_list_reserve(struct FileList* fl, int size);
void file_list_clear(struct FileList *fl);
int  file_list_find(struct FileList* fl, const struct InputFile* in);
void file_list_add_file_ref(struct FileList* fl, struct InputFile *file);
void file_list_add_filename(struct FileList *fl, const char* fname, struct ImageOptions* options);
void file_list_add_filename_with_filter(struct FileList *fl, const char* fname, int pg_min, int pg_max, struct ImageOptions* options);
void file_list_add_ref(struct FileList* dest, struct FileList* src, const char* fname);
void file_list_add_ref_with_filter(struct FileList* dest, struct FileList* src, const char* fname, int pg_min, int pg_max);

struct DjbzOptions* djbz_setting_create(struct DjbzOptions *defaults);
void djbz_setting_destroy(struct DjbzOptions* sett);

void djbz_list_init(struct DjbzList* dl);
void djbz_list_reserve(struct DjbzList *dl, int size);
void djbz_list_clear(struct DjbzList* dl);
void djbz_list_add_option(struct DjbzList* dl, struct DjbzOptions* djbz);

void chunk_file_create(struct ChunkFile* cf, char* fname);
void chunk_file_open(struct ChunkFile* cf);
void chunk_file_close(struct ChunkFile* cf);
void chunk_file_destroy(struct ChunkFile* cf);


void app_options_init(struct AppOptions* opts);
void app_options_free(struct AppOptions* opts);
void app_options_set_djbz_suffix(struct AppOptions* opts, char* suffix);
void app_options_set_output_file(struct AppOptions* opts, const char* filename);

// search for input-files that aren't belong to any Djbz dict
// and assign them to a new Djbz according to default options.
void app_options_autocomplete_djbzs(struct AppOptions* opts);

void app_options_construct_chunk_ids(struct AppOptions* opts);

#ifdef __cplusplus
}
#endif

#endif // APPOPTIONS_H
