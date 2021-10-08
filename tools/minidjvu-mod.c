/*
 * minidjvu-mod.c - an example of using the library
 */

#include <minidjvu-mod/minidjvu-mod.h>
#include "../src/base/mdjvucfg.h" /* for i18n */
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

#include "settings-reader/AppOptions.h"
#include "settings-reader/SettingsReaderAdapter.h"

/* TODO: remove duplicated code */


struct AppOptions options;


/* ========================================================================= */

static void show_usage_and_exit(void)           /* {{{ */
{
    const char *what_it_does = _("encodes/decodes bitonal DjVu files");
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
    printf(_("single-page encoding/decoding:\n"));
    printf(_("    minidjvu-mod [options] <input file> <output file>\n"));
    printf(_("multiple-page encoding:\n"));
    printf(_("    minidjvu-mod [options] <input files> ... <output file>\n"));
    printf(_("Formats supported:\n"));

    printf(_("    DjVu (single-page bitonal), PBM, Windows BMP"));
    if (mdjvu_have_tiff_support())
        printf(_(", TIFF.\n"));
    else
        printf(_("; TIFF support is OFF.\n"));

    printf(_("Options:\n"));
    printf(_("    -A, --Averaging:               compute \"average\" representatives\n"));
    printf(_("    -a <n>, --aggression <n>:      set aggression level (default 100)\n"));
    printf(_("    -c, --clean:                   remove small black pieces\n"));
    printf(_("    -d <n>, --dpi <n>:             set resolution in dots per inch\n"));
    printf(_("    -e, --erosion:                 sacrifice quality to gain in size\n"));
    printf(_("    -i, --indirect:                generate an indirect multipage document\n"));
    printf(_("    -j, --jb2:                     save pages as jb2 chunks instead of djvu.\n"));
    printf(_("                                   implies indirect mode.\n"));
    printf(_("    -l, --lossy:                   use all lossy options (-s -c -m -e -A)\n"));
    printf(_("    -m, --match:                   match and substitute patterns\n"));
    printf(_("    -n, --no-prototypes:           do not search for prototypes\n"));
    printf(_("    -p <n>, --pages-per-dict <n>:  pages per dictionary (default 10)\n"));
    printf(_("    -r, --report:                  report multipage coding progress\n"));
    printf(_("    -s, --smooth:                  remove some bad-looking pixels\n"));
    printf(_("    -S <settings-file>:            read document settings from <settings-file>.\n"));
    printf(_("                                   <settings-file> must be a .txt file.\n"));
    printf(_("                                   Some command-line options may be overridden.\n"));
    printf(_("                                   Details can be found in documentation (man pages).\n"));
#ifdef _OPENMP
    printf(_("    -t <n>, --threads-max <n>:     process pages assigned to different\n"));
    printf(_("                                   dictionaries in up to N parallel threads.\n"));
    printf(_("                                   By default N is equal to the number of \n"));
    printf(_("                                   CPU cores if there are 1 or 2 \n"));
    printf(_("                                   and number of CPU cores minus 1 otherwise.\n"));
    printf(_("                                   Specify -t 1 to disable multithreading\n"));
#endif
    printf(_("    -u, --unbuffered:              unbuffered output to console\n"));
    printf(_("    -v, --verbose:                 print messages about everything\n"));
    printf(_("    -w, --warnings:                do not suppress TIFF warnings\n"));
    printf(_("    -X, --Xtension:                file extension for shared dictionary files.\n"));
    printf(_("                                   Default value is \"djbz\"\n"));
    printf(_("See the man page for detailed description of each option.\n"));

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

static mdjvu_matcher_options_t get_matcher_options(struct DjbzOptions* djbz)
{
    mdjvu_matcher_options_t m_options = NULL;
    if (options.match || options.Match)
    {
        m_options = mdjvu_matcher_options_create();
        mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_PITH_2);
        if (options.Match)
            mdjvu_use_matcher_method(m_options, MDJVU_MATCHER_RAMPAGE);
        mdjvu_set_aggression(m_options, djbz? djbz->aggression : options.default_djbz_options->aggression);
    }
    return m_options;
}

static void sort_and_save_image(mdjvu_image_t image, const char *path, const struct InputFile* in)
{
    mdjvu_error_t error;

    mdjvu_compression_options_t compr_opts = mdjvu_compression_options_create();
    mdjvu_matcher_options_t matcher_opts = get_matcher_options(NULL);
    mdjvu_set_matcher_options(compr_opts, matcher_opts);

    mdjvu_set_verbose(compr_opts, options.verbose);
    mdjvu_set_no_prototypes(compr_opts, options.default_djbz_options->no_prototypes);
    mdjvu_set_averaging(compr_opts, options.default_djbz_options->averaging);
    mdjvu_compress_image(image, compr_opts);
    mdjvu_compression_options_destroy(compr_opts);

    if (options.verbose) printf(_("encoding to `%s'\n"), path);

    const struct ImageOptions* img_opts = in->image_options ? in->image_options : options.default_image_options;

    int res;
    if (!options.save_as_chunk) {
        res = mdjvu_save_djvu_page(image, path, NULL, &error, img_opts->erosion);
    } else {
        res = mdjvu_save_jb2(image, path, &error, img_opts->erosion);
    }

    if (!res) {
        fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
        exit(1);
    }
}

static mdjvu_bitmap_t load_bitmap(struct InputFile* in)
{
    mdjvu_error_t error = NULL;
    mdjvu_bitmap_t bitmap;

    int detect_dpi = 1;
    in->output_dpi = options.default_image_options->dpi;
    if (options.default_image_options->dpi_specified) { // default dpi is overwritten
        detect_dpi = 0;
    }
    if (in->image_options && in->image_options->dpi_specified) { // image dpi is overwritten
        in->output_dpi = in->image_options->dpi;
        detect_dpi = 0;
    }

    const struct ImageOptions* img_opts = in->image_options ? in->image_options : options.default_image_options;

    if (img_opts->is_virtual) {
        return NULL;
    } else if (decide_if_bmp(in->name))
    {
        if (options.verbose) printf(_("loading from Windows BMP file `%s'\n"), in->name);
        bitmap = mdjvu_load_bmp(in->name, &error);
    }
    else if (decide_if_tiff(in->name))
    {
        if (options.verbose) printf(_("loading from TIFF file `%s'\n"), in->name);
        if (!options.warnings)
            mdjvu_disable_tiff_warnings();
        int dpi_from_file = -1;
        bitmap = mdjvu_load_tiff(in->name, detect_dpi ? &dpi_from_file : NULL, &error, in->page);
        if (detect_dpi && dpi_from_file != -1) { // we tried to read dpi from file
            if (dpi_from_file < 20 || dpi_from_file > 2000) {
                if (options.verbose) printf(_("Warning: image %s reports incorrect DPI value (%d) and default resolution %d will be used\n"), in->name, dpi_from_file, in->output_dpi);
            } else {
                in->output_dpi = dpi_from_file;
            }
        }

        if (options.verbose) printf(_("resolution is %d dpi\n"), in->output_dpi);
    }
    else if (decide_if_djvu(in->name))
    {
        mdjvu_image_t image = load_image(in->name);
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
        if (options.verbose) printf(_("loading from PBM file `%s'\n"), in->name);
        bitmap = mdjvu_load_pbm(in->name, &error);
    }

    if (!bitmap)
    {
        fprintf(stderr, "%s: %s\n", in->name, mdjvu_get_error_message(error));
        exit(1);
    }


    if (img_opts->smooth)
    {
        if (options.verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

    return bitmap;
}

static void save_bitmap(mdjvu_bitmap_t bitmap, const char *path, const struct InputFile* in)
{
    mdjvu_error_t error;
    int result;

    if (decide_if_bmp(path))
    {
        if (options.verbose) printf(_("saving to Windows BMP file `%s'\n"), path);
        result = mdjvu_save_bmp(bitmap, path, in->output_dpi, &error);
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
static void decode()
{
    mdjvu_image_t image;    /* a sequence of blits (what is stored in DjVu) */
    mdjvu_bitmap_t bitmap;  /* the result                                   */

    if (options.verbose) {
        printf(_("\nDECODING\n"));
        printf(_("________\n\n"));
    }

    image = load_image(options.file_list.files[0]->name);
    bitmap = mdjvu_render(image);
    mdjvu_image_destroy(image);

    if (options.verbose)
    {
        printf(_("bitmap %d x %d rendered\n"),
               mdjvu_bitmap_get_width(bitmap),
               mdjvu_bitmap_get_height(bitmap));
    }

    const struct InputFile* in = options.file_list.files[0];
    const struct ImageOptions* img_opts = in->image_options ? in->image_options : options.default_image_options;

    if (img_opts->smooth)
    {
        if (options.verbose) printf(_("smoothing the bitmap\n"));
        mdjvu_smooth(bitmap);
    }

    save_bitmap(bitmap, options.output_file, in);
    mdjvu_bitmap_destroy(bitmap);
}


static mdjvu_image_t split_and_destroy(mdjvu_bitmap_t bitmap, const struct InputFile* in)
{
    mdjvu_image_t image;
    if (options.verbose) printf(_("splitting the bitmap into pieces\n"));

    const struct ImageOptions* img_opts = in->image_options ? in->image_options : options.default_image_options;
    if (!bitmap && img_opts->is_virtual) {
        image = mdjvu_image_create(img_opts->virtual_w, img_opts->virtual_h);
        mdjvu_image_set_resolution(image, in->output_dpi);

        if (options.verbose)
        {
            printf(_("creating virtual empty image(%d, %d, dpi: %d) in place of file `%s`\n"),
                   img_opts->virtual_w, img_opts->virtual_h, in->output_dpi, in->name);
        }

        return image;
    }
    image = mdjvu_split(bitmap, in->output_dpi, /* options:*/ NULL);
    mdjvu_bitmap_destroy(bitmap);
    if (options.verbose)
    {
        printf(_("the split image has %d pieces\n"),
               mdjvu_image_get_blit_count(image));
    }

    if (img_opts->clean)
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

static void encode()
{
    mdjvu_bitmap_t bitmap;
    mdjvu_image_t image;

    if (options.verbose) {
        printf(_("\nENCODING\n"));
        printf(_("________\n\n"));
    }

    struct InputFile* in = options.file_list.files[0];

    bitmap = load_bitmap(in);
    image = split_and_destroy(bitmap, in);
    if (options.save_as_chunk) {
        replace_suffix(options.output_file, "jb2");
    }
    sort_and_save_image(image, options.output_file, in);
    mdjvu_image_destroy(image);
}


/* Filtering is nondjvu->nondjvu job. */
static void filter()
{
    mdjvu_bitmap_t bitmap;

    if (options.verbose) {
        printf(_("\nFILTERING\n"));
        printf(_("_________\n\n"));
    }

    struct InputFile* in = options.file_list.files[0];

    bitmap = load_bitmap(in);
    save_bitmap(bitmap, options.output_file, in);
    mdjvu_bitmap_destroy(bitmap);
}

static void print_progress(double progress)
{
    printf(_("[%02d."), (int)progress); //ensure dot as delimiter as in C locale
    printf(_("%02d%%]\n"), (int)(100*(progress - (int)progress)));
}


static void multipage_encode()
{

    if (!decide_if_djvu(options.output_file)) {
        fprintf(stderr, _("when encoding many pages, output file must be DjVu\n"));
        exit(1);
    }

    options.match = 1;

    mdjvu_error_t error;

    if (options.verbose) {
        printf(_("\nMULTIPAGE ENCODING\n"));
        printf(_("__________________\n\n"));
        printf(_("%d pages total\n"), options.file_list.size);
    }

    /* compressing */

#ifdef _OPENMP
    if (!options.max_threads) {
        if (omp_get_num_procs() > 2)
            omp_set_num_threads( omp_get_num_procs() - 1 );
    } else {
        omp_set_num_threads( options.max_threads );
    }
#endif

    int djbz_idx;
    double processed_pages = 0;
    // no need to check _OPENMP as unsupported pragmas are ignored
#pragma omp parallel for schedule(static, 1) shared(processed_pages)
    for (djbz_idx = 0; djbz_idx < options.djbz_list.size; djbz_idx++)
    {
        struct DjbzOptions* const djbz = options.djbz_list.djbzs[djbz_idx];

        mdjvu_compression_options_t compr_opts = mdjvu_compression_options_create();
        mdjvu_matcher_options_t m_opt = get_matcher_options(djbz);
        mdjvu_set_matcher_options(compr_opts, m_opt);

        mdjvu_set_verbose(compr_opts, options.verbose);
        mdjvu_set_no_prototypes(compr_opts, djbz->no_prototypes);
        mdjvu_set_report(compr_opts, options.report);
        mdjvu_set_averaging(compr_opts, djbz->averaging);
        mdjvu_set_report_total_pages(compr_opts, options.file_list.size);

        mdjvu_image_t *images = MDJVU_MALLOCV(mdjvu_image_t, djbz->file_list_ref.size);
        int32 pages_compressed = 0;
        for (int i = 0; i < djbz_idx; i++) {
            pages_compressed += options.djbz_list.djbzs[i]->file_list_ref.size;
        }

        int el = pages_compressed + djbz_idx;

        mdjvu_set_report_start_page(compr_opts, pages_compressed + 1);


        mdjvu_bitmap_t bitmap;
        for (int i = 0; i < djbz->file_list_ref.size; i++)
        {
            struct InputFile* in = djbz->file_list_ref.files[i];
            bitmap = load_bitmap(in);
            images[i] = split_and_destroy(bitmap, in);
            if (options.report) {
                printf(_("Loading: %d of %d completed\n"), pages_compressed + i + 1, options.file_list.size);
                processed_pages += 0.3;
                print_progress(100.0*processed_pages/options.file_list.size);
            }
        }

        mdjvu_image_t dict = mdjvu_compress_multipage(djbz->file_list_ref.size, images, compr_opts);

        if (mdjvu_image_get_bitmap_count(dict) == 0) {
            // do not save empty Djbz (Djbz might be empty for ex., if page list contains only 1 page)
            djbz->do_not_save = 1; // mark element to be skipped at mdjvu_save_djvu_dir
        } else {
            chunk_file_open(&djbz->chunk_file);
            djbz->output_size = mdjvu_file_save_djvu_dictionary(dict, (mdjvu_file_t) djbz->chunk_file.file, 0, &error, djbz->erosion);
            chunk_file_close(&djbz->chunk_file);

            if (!djbz->output_size)
            {
                fprintf(stderr, "%s: %s\n", djbz->chunk_id, mdjvu_get_error_message(error));
                exit(1);
            }
        }

        if (options.report) {
            processed_pages += djbz->file_list_ref.size*0.3;
            print_progress(100.0*processed_pages/options.file_list.size);
        }


        for (int i = 0; i < djbz->file_list_ref.size; i++, el++)
        {
            struct InputFile * in = djbz->file_list_ref.files[i];
            const char * path = in->chunk_id;

            if (options.verbose) {
                if (!djbz->do_not_save) {
                    printf(_("saving page #%d into %s using dictionary %s\n"), pages_compressed + i + 1, path, djbz->chunk_id);
                } else {
                    printf(_("saving page #%d into %s omitting dictionary %s as it's empty\n"), pages_compressed + i + 1, path, djbz->chunk_id);
                }
            }

            const struct ImageOptions* const img_opts = in->image_options ? in->image_options : options.default_image_options;

            const char * dict_chunk = djbz->chunk_id;
            // if chunk is NULL the pages are saved without reference to djbz
            if (djbz->do_not_save) {
                dict_chunk = NULL;
            } else {
                int b = mdjvu_image_get_bitmap_count(images[i]);
                int n = mdjvu_image_get_blit_count(images[i]);
                if (n+b == 0) {
                    dict_chunk = NULL;
                }
            }

            chunk_file_open(&in->chunk_file);
            if (!options.save_as_chunk) {
                in->output_size = mdjvu_file_save_djvu_page(images[i], (mdjvu_file_t) in->chunk_file.file, dict_chunk, 0, &error, img_opts->erosion);
            } else {
                int pos = ftell((FILE *) in->chunk_file.file);
                mdjvu_file_save_jb2(images[i], (mdjvu_file_t) in->chunk_file.file, &error, img_opts->erosion);
                in->output_size = ftell((FILE *) in->chunk_file.file) - pos;
            }
            chunk_file_close(&in->chunk_file);

            if (!in->output_size)
            {
                fprintf(stderr, "%s: %s\n", path, mdjvu_get_error_message(error));
                exit(1);
            }

            mdjvu_image_destroy(images[i]);
            if (options.report) {
                printf(_("Saving: %d of %d completed\n"), pages_compressed + i + 1, options.file_list.size);
                processed_pages += 0.4;
                int val = 100.0*processed_pages / options.file_list.size;
				if (val > 100) {
					val = 100;
				}
				print_progress(val);
            }
        }
        mdjvu_image_destroy(dict);
        MDJVU_FREEV(images);
        mdjvu_compression_options_destroy(compr_opts);
    } //  #pragma omp parallel

    // Saving document directory
    // Let's construct page chunks in a right order.
    int el_size = options.file_list.size + options.djbz_list.size; // max
    char **elements = MDJVU_MALLOCV(char *, el_size);
    int  *sizes     = MDJVU_MALLOCV(int, el_size);
    FILE **files = NULL;
    if (!options.indirect) {
        files    = MDJVU_MALLOCV(FILE *, el_size);
    }

    // assuming options.djbz_list is ordered by min_idx in SettingsReader::constructChunkIDs()
    int el = 0;
    for (int i = 0; i < options.file_list.size; i++) {
        struct InputFile * in = options.file_list.files[i];

        if (in->djbz && !in->djbz->do_not_save) {
            in->djbz->do_not_save = 1;
            elements[el] = in->djbz->chunk_id;
            sizes[el] = in->djbz->output_size;
            if (!options.indirect) {
                chunk_file_open(&in->djbz->chunk_file);
                files[el] = in->djbz->chunk_file.file;
            }
            el++;
        }

        elements[el] = in->chunk_id;
        sizes[el] = in->output_size;
        if (!options.indirect) {
            chunk_file_open(&in->chunk_file);
            files[el] = in->chunk_file.file;
        }
        el++;
    }
    el_size = el;


    if (!options.indirect) {
        FILE* f = fopen(options.output_file, "wb");
        if (!f) {
            fprintf(stderr, "%s: %s\n", options.output_file, (const char *) mdjvu_get_error(mdjvu_error_fopen_write));
            exit(1);
        }

        mdjvu_files_save_djvu_dir(elements, sizes, el_size, (mdjvu_file_t) f, (mdjvu_file_t*) files, el_size, &error);
        fclose(f);

        for (int i = 0; i < options.djbz_list.size; i++) {
            struct DjbzOptions* djbz = options.djbz_list.djbzs[i];
            chunk_file_close(&djbz->chunk_file);
            chunk_file_destroy(&djbz->chunk_file);
        }

        for (int i = 0; i < options.file_list.size; i++) {
            struct InputFile* in = options.file_list.files[i];
            chunk_file_close(&in->chunk_file);
            chunk_file_destroy(&in->chunk_file);
        }

    } else {
        mdjvu_save_djvu_dir(elements, sizes, el_size, options.output_file, &error);
    }

    if (options.report) {
        print_progress(100.0);
    }
    MDJVU_FREEV(files);
    MDJVU_FREEV(elements);
    MDJVU_FREEV(sizes);
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

static void process_options(int argc, char **argv)
{
    int i;
    int settings_file_idx = -1;

    if (argc > 0) {
        // initialize output_file first as we may need to
        // change working dir during that
        // to create chunk files for indirect mode
        app_options_set_output_file(&options, argv[argc - 1]);
        argc--;
    } else {
        show_usage_and_exit();
    }

    for (i = 1; i < argc && argv[i][0] == '-'; i++)
    {
        char *option = argv[i] + 1;
        if (same_option(option, "verbose"))
            options.verbose = 1;
        else if (same_option(option, "smooth"))
            options.default_image_options->smooth = 1;
        else if (same_option(option, "match"))
            options.match = 1;
        else if (same_option(option, "Match"))
            options.Match = 1;
        else if (same_option(option, "no-prototypes"))
            options.default_djbz_options->no_prototypes = 1;
        else if (same_option(option, "erosion"))
            options.default_image_options->erosion = 1;
        else if (same_option(option, "clean"))
            options.default_image_options->clean = 1;
        else if (same_option(option, "warnings"))
            options.warnings = 1;
        else if (same_option(option, "report"))
            options.report = 1;
        else if (same_option(option, "Averaging"))
            options.default_djbz_options->averaging = 1;
        else if (same_option(option, "lossy"))
        {
            options.default_image_options->smooth = 1;
            options.match = 1;
            options.default_image_options->erosion = 1;
            options.default_image_options->clean = 1;
            options.default_djbz_options->averaging = 1;
        }
        else if (same_option(option, "Lossy"))
        {
            options.default_image_options->smooth = 1;
            options.Match = options.match = 1;
            options.default_image_options->erosion = 1;
            options.default_image_options->clean = 1;
            options.default_djbz_options->averaging = 1;
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
            options.default_image_options->dpi = atoi(argv[i]);
            options.default_image_options->dpi_specified = 1;
            if ( (options.default_image_options->dpi < 20) ||
                 (options.default_image_options->dpi > 2000) )
            {
                fprintf(stderr, _("bad resolution\n"));
                exit(2);
            }
        }
        else if (same_option(option, "aggression"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            options.default_djbz_options->aggression = atoi(argv[i]);
            options.match = 1;
        }
        else if (same_option(option, "Xtension"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            app_options_set_djbz_suffix(&options, argv[i]);
        }
        else if (same_option(option, "indirect"))
            options.indirect = 1;
        else if (same_option(option, "jb2"))
            options.save_as_chunk = 1;
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
        else if (same_option(option, "Settings"))
        {
            i++;
            if (i == argc) show_usage_and_exit();
            settings_file_idx = i;
        }
        else
        {
            fprintf(stderr, _("unknown option: %s\n"), argv[i]);
            exit(2);
        }
    }

    if (options.save_as_chunk && !options.indirect) {
        // imply indirect mode.
        options.indirect = 1;
    }


    for (; i < argc; i++) {
        file_list_add_filename(&options.file_list, argv[i], argv[i], 0);
    }

    if (settings_file_idx != -1) {
        // read settings file as a last step
        if (!read_app_options_from_file(argv[settings_file_idx], &options)) {
            exit(3);
        }
    }

    app_options_autocomplete_djbzs(&options);
    app_options_construct_chunk_ids(&options);

    if (!options.output_file || options.file_list.size == 0)
        show_usage_and_exit();
}

int main(int argc, char **argv)
{
    setlocale(LC_ALL, "");
#ifdef HAVE_GETTEXT
    bindtextdomain("minidjvu-mod", LOCALEDIR);
    textdomain("minidjvu-mod");
#endif

    app_options_init(&options);

    process_options(argc, argv);


    if (options.file_list.size > 1)
    {
        multipage_encode();
    } else if (decide_if_djvu(options.output_file)) {
        encode();
    } else if (decide_if_djvu(options.file_list.files[0]->name)) {
        decode();
    } else {
        filter();
    }

    if (options.verbose) printf("\n");
#ifndef NDEBUG
    if (alive_bitmap_counter)
        printf(_("alive_bitmap_counter = %d\n"), alive_bitmap_counter);
#endif

    app_options_free(&options);
    return 0;
}
