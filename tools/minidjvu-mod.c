/*
 * minidjvu-mod.c - an example of using the library
 */

#include <minidjvu-mod/minidjvu-mod.h>
#include "../src/base/mdjvucfg.h" /* for i18n, HAVE_LIBTIFF */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <locale.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifdef _WIN32
#include <fileapi.h>
#endif

#include "settings-reader/DjbzSetting.h"

/* TODO: remove duplicated code */


struct AppOptions
{
    int32 dpi;
    int32 pages_per_dict; /* 0 means infinity */
    int dpi_specified;
    int verbose;
    int smooth;
    int averaging;
    int match;
    int Match;
    int aggression;
    int classifier;
    int erosion;
    int clean;
    int report;
    int no_prototypes;
    int warnings;
    int indirect;
    const char* dict_suffix;
    #ifdef _OPENMP
    int max_threads;
    #endif
};

struct AppOptions options;

static void init_default_options() {
    /* options */
    options.dpi = 300;
    options.pages_per_dict = 10; /* 0 means infinity */
    options.dpi_specified = 0;
    options.verbose = 0;
    options.smooth = 0;
    options.averaging = 0;
    options.match = 0;
    options.Match = 0;
    options.aggression = 100;
    options.classifier = 3;
    options.erosion = 0;
    options.clean = 0;
    options.report = 0;
    options.no_prototypes = 0;
    options.warnings = 0;
    options.indirect = 0;
    options.dict_suffix = NULL;
    #ifdef _OPENMP
    options.max_threads = 0;
    #endif
}

///* options */
//int32 dpi = 300;
//int32 pages_per_dict = 10; /* 0 means infinity */
//int dpi_specified = 0;
//int verbose = 0;
//int smooth = 0;
//int averaging = 0;
//int match = 0;
//int Match = 0;
//int aggression = 100;
//int classifier = 3;
//int erosion = 0;
//int clean = 0;
//int report = 0;
//int no_prototypes = 0;
//int warnings = 0;
//int indirect = 0;
//const char* dict_suffix = NULL;
//#ifdef _OPENMP
//int max_threads = 0;
//#endif

/* ========================================================================= */

/* file name template routines (for multipage encoding) {{{ */

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

static char *get_page_or_dict_name(char **elements, int cnt, const char *fname)
{
    int i, extpos, same=-1;
    char *page_name, *pattern;
    
    extpos = get_ext_delim_pos(fname);
    page_name = MDJVU_MALLOCV(char, extpos + 10);
    memset(page_name,'\0',extpos + 6);
    if (extpos > 0)
        strncpy(page_name, fname, extpos-1);
    strcat(page_name, ".djvu");

    for (i=0; i<cnt; i++ )
    {
        if (strcmp(page_name, elements[i]) == 0)
        {
            same = i;
            break;
        }
    }
    
    if (same != -1)
    {
        int previdx=0, idx=0, res;
        /* discard the extension */
        page_name[extpos-1] = '\0';
        
        pattern = MDJVU_MALLOCV(char, extpos + 10);
        strcpy(pattern, page_name);
        strcat(pattern, "#%d.");

        for (i=same; i<cnt; i++)
        {
            if (mdjvu_ends_with_ignore_case(elements[i],".djvu"))
            {
                res = sscanf(elements[i],pattern,&idx);
                if (res && idx > previdx) previdx = idx;
            }
        }
        if (idx == 999)
        {
            fprintf(stderr, _("Cannot generate a unique name for %s\n"), fname);
            exit(1);
        }
        sprintf(page_name + (extpos - 1),"#%03d.djvu",idx+1);
        MDJVU_FREEV(pattern);
    }
    return(page_name);
}

static void replace_suffix(char *name, const char *suffix)
{
    int len = strlen(name);
    
    name[len-4] = '\0';
    strcat(name, suffix);
}

/* file name template routines (for multipage encoding) }}} */

/* ========================================================================= */

static void show_usage_and_exit(void)           /* {{{ */
{
    const char *what_it_does = _("encode/decode bitonal DjVu files");
    if (strcmp(MDJVU_VERSION, mdjvu_get_version()))
    {
        printf(_("minidjvu-mod - %s\n"), what_it_does);
        printf(_("Warning: program and library version mismatch:\n"));
        printf(_("    program version %s, library version %s.\n\n"), MDJVU_VERSION, mdjvu_get_version());

    }
    else
    {
        printf("minidjvu-mod %s - %s\n", MDJVU_VERSION, what_it_does);

    }
    printf(_("Usage:\n"));
    printf(_("single page encoding/decoding:\n"));
    printf(_("    minidjvu-mod [options] <input file> <output file>\n"));
    printf(_("multiple pages encoding:\n"));
    printf(_("    minidjvu-mod [options] <input file> ... <output file>\n"));
    printf(_("Formats supported:\n"));

    printf(_("    DjVu (single-page bitonal), PBM, Windows BMP"));
    if (mdjvu_have_tiff_support())
        printf(_(", TIFF.\n"));
    else
        printf(_("; TIFF support is OFF.\n"));

    printf(_("Options:\n"));
    printf(_("    -A, --Averaging:               compute \"average\" representatives\n"));
    printf(_("    -a <n>, --aggression <n>:      set aggression level (default 100)\n"));
    printf(_("    -C <n>, --Classifier <n>:      set symbols classifier mode (default 3)\n"));
    printf(_("                                   1 - behave similar to original one.\n"));
    printf(_("                                   2 - make additional efforts to achieve\n"));
    printf(_("                                       better compression. This require\n"));
    printf(_("                                       more CPU time and much more RAM\n"));
    printf(_("                                       as cache is allocated per thread.\n"));
    printf(_("                                   3 - similar to 2 but takes even more\n"));
    printf(_("                                       CPU time to achieve maximum level\n"));
    printf(_("                                       of document compression (same RAM).\n"));
    printf(_("                                   BE VERY CAREFUL with modes 2 and 3 as\n"));
    printf(_("                                       they can slow down your machine or\n"));
    printf(_("                                       overflow available RAM size. \n"));
    printf(_("                                       You may decrease number of threads\n"));
    printf(_("                                       to save some RAM at the cost of time.\n"));
    printf(_("                                       Or decrease pages per dict at a cost\n"));
    printf(_("                                       of filesize (not recommended).\n"));
    printf(_("    -c, --clean                    remove small black pieces\n"));
    printf(_("    -d <n> --dpi <n>:              set resolution in dots per inch\n"));
    printf(_("    -e, --erosion                  sacrifice quality to gain in size\n"));
    printf(_("    -i, --indirect:                generate an indirect multipage document\n"));
    printf(_("    -l, --lossy:                   use all lossy options (-s -c -m -e -A)\n"));
    printf(_("    -m, --match:                   match and substitute patterns\n"));
    printf(_("    -n, --no-prototypes:           do not search for prototypes\n"));
    printf(_("    -p <n>, --pages-per-dict <n>:  pages per dictionary (default 10)\n"));
    printf(_("    -r, --report:                  report multipage coding progress\n"));
    printf(_("    -s, --smooth:                  remove some badly looking pixels\n"));
#ifdef _OPENMP
    printf(_("    -t <n>, --threads-max <n>:     process pages assigned to a different\n"));
    printf(_("                                   dictionaries in up to N parallel threads.\n"));
    printf(_("                                   By default N is equal to the number of \n"));
    printf(_("                                   CPU cores in case there're 1 or 2 \n"));
    printf(_("                                   and number of CPU cores minus 1 otherwise\n"));
    printf(_("                                   Specify -t 1 to disable multithreading\n"));
#endif
    printf(_("    -u, --unbuffered:              unbuffered output to console\n"));
    printf(_("    -v, --verbose:                 print messages about everything\n"));
    printf(_("    -w, --warnings:                do not suppress TIFF warnings\n"));
    printf(_("    -X, --Xtension:                file extension for shared dictionary files\n"));
    printf(_("See the man page for detailed description of each option.\n"));
    
    printf(test());
   /// struct DjbzSetting** settings;
 //   read_settings("file::\\/tmp/sed.ex", settings);
    exit(2);
}                   /* }}} */

static int decide_if_bmp(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".bmp");
}

static int decide_if_djvu(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".djvu")
        || mdjvu_ends_with_ignore_case(path, ".djv");
}

static int decide_if_tiff(const char *path)
{
    return mdjvu_ends_with_ignore_case(path, ".tiff")
        || mdjvu_ends_with_ignore_case(path, ".tif");
}

/* ========================================================================= */

static mdjvu_image_t load_image(const char *path)
{
    mdjvu_error_t error;
    mdjvu_image_t image;

    if (options.verbose) printf(_("loading a DjVu page from `%s'\n"), path);
    image = mdjvu_load_djvu_page(path, &error);
    if (!image)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
    if (options.verbose)
    {
        printf(_("loaded; the page has %d bitmaps and %d blits\n"),
               mdjvu_image_get_bitmap_count(image),
               mdjvu_image_get_blit_count(image));
    }
    return image;
}

static mdjvu_matcher_options_t get_matcher_options(void)
{
    mdjvu_matcher_options_t m_options = NULL;
    if (options.match || options.Match)
    {
        m_options = mdjvu_matcher_options_create();
        mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_PITH_2);
        if (options.Match)
            mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_RAMPAGE);
        mdjvu_set_aggression(m_options, options.aggression);
    }
    return m_options;
}

static mdjvu_classify_options_t get_classify_options(void)
{
    mdjvu_classify_options_t m_options = NULL;
	m_options = mdjvu_classify_options_create();
    mdjvu_set_classifier(m_options, options.classifier);
	return m_options;
}

static void sort_and_save_image(mdjvu_image_t image, const char *path)
{
    mdjvu_error_t error;

    mdjvu_compression_options_t compr_opts = mdjvu_compression_options_create();
    mdjvu_matcher_options_t matcher_opts = get_matcher_options();
    mdjvu_set_classify_options(matcher_opts, get_classify_options());
    mdjvu_set_matcher_options(compr_opts, matcher_opts);

    mdjvu_set_clean(compr_opts, options.clean);
    mdjvu_set_verbose(compr_opts, options.verbose);
    mdjvu_set_no_prototypes(compr_opts, options.no_prototypes);
    mdjvu_set_averaging(compr_opts, options.averaging);
    mdjvu_compress_image(image, compr_opts);
    mdjvu_compression_options_destroy(compr_opts);

    if (options.verbose) printf(_("encoding to `%s'\n"), path);

    if (!mdjvu_save_djvu_page(image, path, NULL, &error, options.erosion))
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
}

static mdjvu_bitmap_t load_bitmap(const char *path, int tiff_idx)
{
    mdjvu_error_t error = NULL;
    mdjvu_bitmap_t bitmap;

    if (decide_if_bmp(path))
    {
        if (options.verbose) printf(_("loading from Windows BMP file `%s'\n"), path);
        bitmap = mdjvu_load_bmp(path, &error);
    }
    else if (decide_if_tiff(path))
    {
        if (options.verbose) printf(_("loading from TIFF file `%s'\n"), path);
        if (!options.warnings)
            mdjvu_disable_tiff_warnings();
        if (options.dpi_specified)
            bitmap = mdjvu_load_tiff(path, NULL, &error, tiff_idx);
        else
            bitmap = mdjvu_load_tiff(path, &options.dpi, &error, tiff_idx);
        if (options.verbose) printf(_("resolution is %d dpi\n"), options.dpi);
    }
    else if (decide_if_djvu(path))
    {
        mdjvu_image_t image = load_image(path);
        bitmap = mdjvu_render(image);
        mdjvu_image_destroy(image);
        if (options.verbose)
        {
            printf(_("bitmap %d x %d rendered\n"),
                   mdjvu_bitmap_get_width(bitmap),
                   mdjvu_bitmap_get_height(bitmap));
        }
    }
    else
    {
        if (options.verbose) printf(_("loading from PBM file `%s'\n"), path);
        bitmap = mdjvu_load_pbm(path, &error);
    }

    if (!bitmap)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }

    if (options.smooth)
    {
        if (options.verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

    return bitmap;
}

static void save_bitmap(mdjvu_bitmap_t bitmap, const char *path)
{
    mdjvu_error_t error;
    int result;

    if (decide_if_bmp(path))
    {
        if (options.verbose) printf(_("saving to Windows BMP file `%s'\n"), path);
        result = mdjvu_save_bmp(bitmap, path, options.dpi, &error);
    }
    else if (decide_if_tiff(path))
    {
        if (options.verbose) printf(_("saving to TIFF file `%s'\n"), path);
        if (!options.warnings)
            mdjvu_disable_tiff_warnings();
        result = mdjvu_save_tiff(bitmap, path, &error);
    }
    else
    {
        if (options.verbose) printf(_("saving to PBM file `%s'\n"), path);
        result = mdjvu_save_pbm(bitmap, path, &error);
    }

    if (!result)
    {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
}

/* ========================================================================= */

static void decode(int argc, char **argv)
{
    mdjvu_image_t image;    /* a sequence of blits (what is stored in DjVu) */
    mdjvu_bitmap_t bitmap;  /* the result                                   */

    if (options.verbose) {
        printf(_("\nDECODING\n"));
        printf(_("________\n\n"));
    }

    image = load_image(argv[1]);
    bitmap = mdjvu_render(image);
    mdjvu_image_destroy(image);

    if (options.verbose)
    {
        printf(_("bitmap %d x %d rendered\n"),
               mdjvu_bitmap_get_width(bitmap),
               mdjvu_bitmap_get_height(bitmap));
    }

    if (options.smooth)
    {
        if (options.verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

    save_bitmap(bitmap, argv[2]);
    mdjvu_bitmap_destroy(bitmap);
}


static mdjvu_image_t split_and_destroy(mdjvu_bitmap_t bitmap)
{
    mdjvu_image_t image;
    if (options.verbose) printf(_("splitting the bitmap into pieces\n"));
    image = mdjvu_split(bitmap, options.dpi, /* options:*/ NULL);
    mdjvu_bitmap_destroy(bitmap);
    if (options.verbose)
    {
        printf(_("the split image has %d pieces\n"),
                mdjvu_image_get_blit_count(image));
    }
    if (options.clean)
    {
        if (options.verbose) printf(_("cleaning\n"));
        mdjvu_clean(image);
        if (options.verbose)
        {
            printf(_("the cleaned image has %d pieces\n"),
                    mdjvu_image_get_blit_count(image));
        }
    }
    return image;
}


static void encode(int argc, char **argv)
{
    mdjvu_bitmap_t bitmap;
    mdjvu_image_t image;

    if (options.verbose) {
        printf(_("\nENCODING\n"));
        printf(_("________\n\n"));
    }

    bitmap = load_bitmap(argv[1], 0);

    image = split_and_destroy(bitmap);
    sort_and_save_image(image, argv[2]);
    mdjvu_image_destroy(image);
}


/* Filtering is nondjvu->nondjvu job. */
static void filter(int argc, char **argv)
{
    mdjvu_bitmap_t bitmap;

    if (options.verbose) {
        printf(_("\nFILTERING\n"));
        printf(_("_________\n\n"));
    }

    bitmap = load_bitmap(argv[1], 0);
    save_bitmap(bitmap, argv[2]);
    mdjvu_bitmap_destroy(bitmap);
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
    return strip(strip(path, '\\'), '/');
}


static void multipage_encode(int n, char **pages, char *outname, uint32 multipage_tiff)
{
    int ndicts = (options.pages_per_dict <= 0)? 1 :
                                        (n % options.pages_per_dict > 0) ?  (int) fabs((float)n/options.pages_per_dict) + 1:
                                                                    (int) fabs((float)n/options.pages_per_dict);
    const int el_size = n + ndicts;
    char **elements = MDJVU_MALLOCV(char *, el_size);
    int  *sizes     = MDJVU_MALLOCV(int, el_size);
    mdjvu_compression_options_t compr_opts;
    mdjvu_error_t error;

    FILE *f;
    FILE ** tfs = MDJVU_MALLOCV(FILE *, ndicts);
#ifdef _WIN32
    char* tempfilenames = (char*) MDJVU_MALLOCV(char, ndicts*(MAX_PATH+1));
    memset(tempfilenames, 0, ndicts*(MAX_PATH+1));
#endif

    options.match = 1;

    if (!decide_if_djvu(outname)) {
        fprintf(stderr, _("when encoding many pages, output file must be DjVu\n"));
        exit(1);
    }

    if (!options.indirect) {
#ifdef _WIN32
        char pathname[MAX_PATH+1];
        if (GetTempPathA(MAX_PATH+1, pathname) == 0) {
            fprintf(stderr, _("Could not create a temporary file. (GetTempPathA)\n"));
            exit(1);
        }
        char prefix[4] = {'d', 'j', 'v', 0 };
#endif

        for (int i = 0; i < ndicts; i++) {
#ifndef _WIN32
            tfs[i] = tmpfile();
#else
            char* tempfilename = tempfilenames + i*(MAX_PATH+1);
            if (GetTempFileNameA(pathname, prefix, 0, tempfilename) == 0) {
                fprintf(stderr, _("Could not create a temporary file. (GetTempFileNameA)\n"));
                exit(1);
            }
            tfs[i] = fopen(tempfilename, "w+b");
#endif
            if (!tfs[i])
            {
                fprintf(stderr, _("Could not create a temporary file\n"));
                exit(1);
            }
        }
    }

    if (options.verbose) {
        printf(_("\nMULTIPAGE ENCODING\n"));
        printf(_("__________________\n\n"));
        printf(_("%d pages total\n"), n);
    }

    compr_opts = mdjvu_compression_options_create();
    mdjvu_matcher_options_t m_opt = get_matcher_options();
    mdjvu_set_classify_options(m_opt, get_classify_options());
    mdjvu_set_matcher_options(compr_opts, m_opt);

    mdjvu_set_clean(compr_opts, options.clean);
    mdjvu_set_verbose(compr_opts, options.verbose);
    mdjvu_set_no_prototypes(compr_opts, options.no_prototypes);
    mdjvu_set_report(compr_opts, options.report);
    mdjvu_set_averaging(compr_opts, options.averaging);
    mdjvu_set_report_total_pages(compr_opts, n);

    /* compressing */
    if (options.pages_per_dict <= 0) options.pages_per_dict = n;
    if (options.pages_per_dict > n) options.pages_per_dict = n;

#ifdef _OPENMP
    if (!options.max_threads) {
        if (omp_get_num_procs() > 2)
            omp_set_num_threads( omp_get_num_procs() - 1 );
    } else {
        omp_set_num_threads( options.max_threads );
    }
#endif

    // initialize elements with filenames as it can't be done in parallel blocks without critical sections
    int el = 0;
    for (int j = 0; j < ndicts; j++) {
        const int32 pages_compressed = j*options.pages_per_dict;
        const int32 pages_to_compress = pages_compressed + options.pages_per_dict > n ? n - pages_compressed
                                                                                      : options.pages_per_dict;

        char * path = get_page_or_dict_name(elements, el, strip_dir(pages[multipage_tiff ? 0 : pages_compressed]));
        char * dict_name = MDJVU_MALLOCV(char, strlen(path) + strlen(options.dict_suffix) - 2);
        strcpy(dict_name, path);
        replace_suffix(dict_name, options.dict_suffix);
        elements[el++] = dict_name;
        for (int i = 0; i < pages_to_compress; i++)
        {
            if (i > 0) {
                path = get_page_or_dict_name(elements, el, strip_dir(pages[multipage_tiff ? 0 : pages_compressed + i]));
            }
            elements[el++] = path;
        }
    }

    int block;
    int processed_pages = 0;
// no need to check _OPENMP as unsupported pragmas are ignored
#pragma omp parallel for schedule(static, 1) shared(processed_pages)
    for (block = 0; block < ndicts; block++) {
        mdjvu_image_t *images = MDJVU_MALLOCV(mdjvu_image_t, options.pages_per_dict);
        int32 pages_compressed = block*options.pages_per_dict;
        int32 pages_to_compress = pages_compressed + options.pages_per_dict > n ? n - pages_compressed : options.pages_per_dict;
        int el = pages_compressed + block;

        mdjvu_set_report_start_page(compr_opts, pages_compressed + 1);


        mdjvu_bitmap_t bitmap;
        for (int i = 0; i < pages_to_compress; i++)
        {
            if (multipage_tiff)
                bitmap = load_bitmap(pages[0], pages_compressed + i);
            else
                bitmap = load_bitmap(pages[pages_compressed + i], 0);
            images[i] = split_and_destroy(bitmap);
            if (options.report)
                printf(_("Loading: %d of %d completed\n"), pages_compressed + i + 1, n);
        }

        mdjvu_image_t dict = mdjvu_compress_multipage(pages_to_compress, images, compr_opts);

        const char * dict_name = elements[el];
        if (!options.indirect)
            sizes[el] = mdjvu_file_save_djvu_dictionary(dict, (mdjvu_file_t) tfs[block], 0, &error, options.erosion);
        else
            sizes[el] = mdjvu_save_djvu_dictionary(dict, dict_name, &error, options.erosion);

        if (!sizes[el])
        {
            fprintf(stderr, "%s: %s\n", dict_name, mdjvu_get_error_message(error));
            exit(1);
        }

        el++;

        for (int i = 0; i < pages_to_compress; i++, el++)
        {
            const char * path = elements[el];

            if (options.verbose)
                printf(_("saving page #%d into %s using dictionary %s\n"), pages_compressed + i + 1, path, dict_name);

            if (!options.indirect)
                sizes[el] = mdjvu_file_save_djvu_page(images[i], (mdjvu_file_t) tfs[block], strip_dir(dict_name), 0, &error, options.erosion);
            else
                sizes[el] = mdjvu_save_djvu_page(images[i], path, strip_dir(dict_name), &error, options.erosion);
            if (!sizes[el])
            {
                fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
                exit(1);
            }

            mdjvu_image_destroy(images[i]);
            if (options.report) {
                printf(_("Saving: %d of %d completed\n"), pages_compressed + i + 1, n);
                processed_pages++;
                float res = 100.0*processed_pages/n;
                printf(_("[%02d."), (int)res); //ensure dot as delimiter as in C locale
                printf(_("%02d%%]\n"), (int)(100*(res - (int)res)));
            }
        }
        mdjvu_image_destroy(dict);
        //        pages_compressed += pages_to_compress;
        MDJVU_FREEV(images);
    } //  #pragma omp parallel

    if (!options.indirect)
    {
        f = fopen(outname, "wb");
        if (!f)
        {
            fprintf(stderr, "%s: %s\n", outname, (const char *) mdjvu_get_error(mdjvu_error_fopen_write));
            exit(1);
        }
        mdjvu_files_save_djvu_dir(elements, sizes, el_size, (mdjvu_file_t) f, (mdjvu_file_t*) tfs, ndicts, &error);
        for (int i = 0; i < ndicts; i++) {
            fclose(tfs[i]);
#ifdef _WIN32
            char* tempfilename = tempfilenames + i*(MAX_PATH+1);
            remove(tempfilename);
#endif
        }
        MDJVU_FREEV(tfs);
#ifdef _WIN32
        MDJVU_FREEV(tempfilenames);
#endif
        fclose(f);
    }
    else
        mdjvu_save_djvu_dir(elements, sizes, el_size, outname, &error);

    for (int i=0; i < el_size; i++) MDJVU_FREEV(elements[i]);
    MDJVU_FREEV(elements);
    MDJVU_FREEV(sizes);

    /* destroying */
    mdjvu_compression_options_destroy(compr_opts);
}

/* same_option(foo, "opt") returns 1 in three cases:
 *
 *      foo is "o" (first letter of opt)
 *      foo is "opt"
 *      foo is "-opt"
 */
static int same_option(const char *option, const char *s)
{
    if (option[0] == s[0] && !option[1]) return 1;
    if (!strcmp(option, s)) return 1;
    if (option[0] == '-' && !strcmp(option + 1, s)) return 1;
    return 0;
}

static int process_options(int argc, char **argv)
{
    int i;
    for (i = 1; i < argc && argv[i][0] == '-'; i++)
    {
        char *option = argv[i] + 1;
        if (same_option(option, "verbose"))
            options.verbose = 1;
        else if (same_option(option, "smooth"))
            options.smooth = 1;
        else if (same_option(option, "match"))
            options.match = 1;
        else if (same_option(option, "Match"))
            options.Match = 1;
        else if (same_option(option, "no-prototypes"))
            options.no_prototypes = 1;
        else if (same_option(option, "erosion"))
            options.erosion = 1;
        else if (same_option(option, "clean"))
            options.clean = 1;
        else if (same_option(option, "warnings"))
            options.warnings = 1;
        else if (same_option(option, "report"))
            options.report = 1;
        else if (same_option(option, "Averaging"))
            options.averaging = 1;
        else if (same_option(option, "lossy"))
        {
            options.smooth = 1;
            options.match = 1;
            options.erosion = 1;
            options.clean = 1;
            options.averaging = 1;
        }
        else if (same_option(option, "Lossy"))
        {
            options.smooth = 1;
            options.Match = options.match = 1;
            options.erosion = 1;
            options.clean = 1;
            options.averaging = 1;
        }
        else if (same_option(option, "pages-per-dict"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.pages_per_dict = atoi(argv[i]);
            if (options.pages_per_dict < 0)
            {
                fprintf(stderr, _("bad --pages-per-dict value\n"));
                exit(2);
            }
        }
        else if (same_option(option, "dpi"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.dpi = atoi(argv[i]);
            options.dpi_specified = 1;
            if (options.dpi < 20 || options.dpi > 2000)
            {
                fprintf(stderr, _("bad resolution\n"));
                exit(2);
            }
        }
        else if (same_option(option, "aggression"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.aggression = atoi(argv[i]);
            options.match = 1;
        }
        else if (same_option(option, "Classifier"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.classifier = atoi(argv[i]);
        }
        else if (same_option(option, "Xtension"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.dict_suffix = argv[i];
        }
        else if (same_option(option, "indirect"))
            options.indirect = 1;
#ifdef _OPENMP
        else if (same_option(option, "threads-max"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.max_threads = atoi(argv[i]);
        }
#endif
        else if (same_option(option, "unbuffered"))
        {
            setbuf(stdout, NULL);
            setbuf(stderr, NULL);
        }
        else
        {
            fprintf(stderr, _("unknown option: %s\n"), argv[i]);
            exit(2);
        }
    }
    return i;
}

int main(int argc, char **argv)
{
    init_default_options();

    int arg_start;
#ifdef HAVE_LIBTIFF
    int tiff_cnt;
#endif

    setlocale(LC_ALL, "");
#ifdef HAVE_GETTEXT
    bindtextdomain("minidjvu-mod", LOCALEDIR);
    textdomain("minidjvu-mod");
#endif


    arg_start = process_options(argc, argv);
    if ( options.dict_suffix == NULL ) options.dict_suffix = "djbz";

    argc -= arg_start - 1;
    argv += arg_start - 1;

    if (argc < 3)
        show_usage_and_exit();

    if (argc > 3)
    {
        multipage_encode(argc - 2, argv + 1, argv[argc - 1], 0);
    }
#ifdef HAVE_LIBTIFF
    else if (decide_if_tiff(argv[1]) && (tiff_cnt = mdjvu_get_tiff_page_count(argv[1])) > 1 )
    {
        multipage_encode(tiff_cnt, argv + 1, argv[argc - 1], 1);
    }
#endif
    else if (decide_if_djvu(argv[2]))
    {
        encode(argc, argv);
    }
    else
    {
        if (decide_if_djvu(argv[1]))
            decode(argc, argv);
        else
            filter(argc, argv);
    }

    if (options.verbose) printf("\n");
	#ifndef NDEBUG
		if (alive_bitmap_counter)
           printf(_("alive_bitmap_counter = %d\n"), alive_bitmap_counter);
    #endif
    return 0;
}
