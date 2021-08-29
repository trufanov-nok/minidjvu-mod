/*
 * compress.c - recommended sequences of algorithms invocations to compress
 */


#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

struct MinidjvuCompressionOptions
{
    int clean;
    int verbose;
    int no_prototypes;
    int averaging;
    int report;
    int report_start_page;
    int report_total_pages;
    mdjvu_matcher_options_t matcher_options;
};

MDJVU_IMPLEMENT mdjvu_compression_options_t mdjvu_compression_options_create()
{
    mdjvu_compression_options_t opt = (mdjvu_compression_options_t)
        malloc(sizeof(struct MinidjvuCompressionOptions));
    mdjvu_init();
    opt->clean = 0;
    opt->verbose = 0;
    opt->report = 0;
    opt->averaging = 0;
    opt->no_prototypes = 0;
    opt->matcher_options = NULL;
    return opt;
}

MDJVU_IMPLEMENT void mdjvu_compression_options_destroy(mdjvu_compression_options_t opt)
{
    if (opt->matcher_options)
        mdjvu_matcher_options_destroy(opt->matcher_options);
    free(opt);
}

MDJVU_IMPLEMENT void mdjvu_set_matcher_options(mdjvu_compression_options_t opt, mdjvu_matcher_options_t v)
{
    if (opt->matcher_options)
        mdjvu_matcher_options_destroy(opt->matcher_options);
    opt->matcher_options = v;
}

MDJVU_IMPLEMENT void mdjvu_set_clean(mdjvu_compression_options_t opt, int v)
    {opt->clean = v;}
MDJVU_IMPLEMENT void mdjvu_set_verbose(mdjvu_compression_options_t opt, int v)
    {opt->verbose = v;}
MDJVU_IMPLEMENT void mdjvu_set_averaging(mdjvu_compression_options_t opt, int v)
    {opt->averaging = v;}
MDJVU_IMPLEMENT void mdjvu_set_no_prototypes(mdjvu_compression_options_t opt, int v)
    {opt->no_prototypes = v;}
MDJVU_IMPLEMENT void mdjvu_set_report(mdjvu_compression_options_t opt, int v)
    {opt->report = v;}
MDJVU_IMPLEMENT void mdjvu_set_report_start_page(mdjvu_compression_options_t opt, int v)
    {opt->report_start_page = v;}
MDJVU_IMPLEMENT void mdjvu_set_report_total_pages(mdjvu_compression_options_t opt, int v)
    {opt->report_total_pages = v;}

static void find_substitutions(mdjvu_image_t image,
                                              struct MinidjvuCompressionOptions *opt)
{
    mdjvu_matcher_options_t m_opt = opt->matcher_options;
    int32 i, n = mdjvu_image_get_bitmap_count(image);
    int32 *tags = (int32 *) malloc(n * sizeof(int32));
    int32 max_tag = mdjvu_classify_bitmaps(image, tags, m_opt, /* centers_needed: */ opt->averaging);
    mdjvu_bitmap_t *representatives = (mdjvu_bitmap_t *)
        calloc(max_tag + 1 /* cause starts with 1 */, sizeof(mdjvu_bitmap_t));
    int32 *cx = (int32 *) malloc(n * sizeof(int32));
    int32 *cy = (int32 *) malloc(n * sizeof(int32));

    if (!mdjvu_image_has_substitutions(image))
       mdjvu_image_enable_substitutions(image);
    if (!opt->averaging)
    {
        for (i = 0; i < n; i++)
        {
            if (!representatives[tags[i]])
                representatives[tags[i]] = mdjvu_image_get_bitmap(image, i);
        }
    }
    else
    {
        /* Average */
        mdjvu_bitmap_t *sources = (mdjvu_bitmap_t *) calloc(n, sizeof(mdjvu_bitmap_t));
        for (i = 1; i <= max_tag; i++)
        {
            int j, sources_found = 0;
            mdjvu_bitmap_t rep;
    
            for (j = 0; j < n; j++)
            {
                if (tags[j] == i)
                    sources[sources_found++] = mdjvu_image_get_bitmap(image, j);
            }
    
            for (j = 0; j < sources_found; j++)
                mdjvu_image_get_center(image, sources[j], &cx[j], &cy[j]);
    
            rep = mdjvu_average(sources, sources_found, cx, cy);
            mdjvu_image_add_bitmap(image, rep);
            mdjvu_image_set_substitution(image, rep, rep);
            representatives[i] = rep;
        }
        
        mdjvu_image_disable_centers(image);
        free(sources);
    }
    assert(mdjvu_image_check_indices(image));
    free(cx);
    free(cy);


    for (i = 0; i < n; i++)
    {
        if (!tags[i]) continue; /* that's for images with no-subst flag */
        mdjvu_image_set_substitution(image,
                                     mdjvu_image_get_bitmap(image, i),
                                     representatives[tags[i]]);
    }

    free(representatives);
    free(tags);
}


static int32 count_prototypes(mdjvu_image_t image)
{
    int32 i, s = 0, n = mdjvu_image_get_bitmap_count(image);
    for (i = 0; i < n; i++)
    {
        mdjvu_bitmap_t bmp = mdjvu_image_get_bitmap(image, i);
        if (mdjvu_image_get_prototype(image, bmp))
            s++;
    }
    return s;
}


MDJVU_IMPLEMENT void mdjvu_compress_image(mdjvu_image_t image, mdjvu_compression_options_t opt)
{
    mdjvu_compression_options_t options;
    if (opt)
        options = opt;
    else
        options = mdjvu_compression_options_create();

    if (options->verbose) puts(_("deciding what pieces are letters"));
    mdjvu_calculate_not_a_letter_flags(image);

    if (options->verbose) puts(_("sorting blits"));
    mdjvu_sort_blits(image);
    if (options->verbose) puts(_("sorting bitmaps"));
    mdjvu_image_sort_bitmaps(image);

    if (options->matcher_options)
    {
        if (options->verbose) puts(_("matching patterns"));
        find_substitutions(image, options);
        
        if (options->verbose) puts(_("adjusting substitution coordinates"));
        mdjvu_adjust(image);
        
        if (options->verbose) puts(_("removing unused bitmaps"));
        mdjvu_image_remove_unused_bitmaps(image);
        
        if (options->verbose)
        {
            printf(_("the image now has %d bitmaps\n"),
                   mdjvu_image_get_bitmap_count(image));
        }
        
        if (options->averaging)
        {
           if (options->verbose) puts(_("sorting bitmaps (again)"));
           mdjvu_image_sort_bitmaps(image);
        }
    }

    if (options->no_prototypes)
    {
        mdjvu_image_enable_prototypes(image);
    }
    else
    {
        if (options->verbose) puts(_("finding prototypes"));
        mdjvu_find_prototypes(image);
        if (options->verbose)
        {
            printf(_("%d bitmaps have prototypes\n"),
                   count_prototypes(image));
        }
    }

    if (!opt)
        mdjvu_compression_options_destroy(options);
}


/* -----------------------   MULTIPAGE STUFF   ------------------------------ */


/* Fill the dictionary with the bitmaps chosen with dictionary_flags.
 * The `representatives' array is filled with new dictionary letters as needed.
 *
 * Arguments:
 *      max_tag          - obviously, the maximal tag;         tag 0 is unused
 *      representatives  - [0..max_tag] array of letters;      tag 0 is unused
 *      dictionary_flags - [0..max_tag] array of 1/0 flags;    tag 0 is unused
 *
 * Returns:
 *      the newly created dictionary (a 0 x 0 image to be destroyed later)
 *
 * FIXME: cloning bitmaps is a waste of memory - better to detach them instead.
 */

static mdjvu_image_t get_dictionary(int32 max_tag,
                                    mdjvu_bitmap_t *representatives)
{
    mdjvu_image_t dictionary = mdjvu_image_create(0,0); /* 0 x 0 image */
    int32 tag;
    for (tag = 1; tag <= max_tag; tag++)
    {
        mdjvu_bitmap_t rep;
        if ((rep = representatives[tag]) == NULL) continue;
        mdjvu_bitmap_t clone = mdjvu_bitmap_clone(rep);
        representatives[tag] = clone;
        mdjvu_image_add_bitmap(dictionary, clone);
    }
    return dictionary;
}


/* -------------------------------------------------------------------------- */

/* Set `substitution' pointers in each page to representatives.
 *
 * Arguments:
 *      n                   - number of pages
 *      pages
 *      total_bitmaps_count - sum of all bitmap counts on each page
 *      tags                - an array of tags (to index `representatives')
 *      representatives     - [0..max_tag] array of bitmaps (tag 0 unused)
 */

static void set_substitutions(int n,
                              mdjvu_image_t *pages,
                              int32 total_bitmaps_count,
                              int32 *tags,
                              mdjvu_bitmap_t *representatives)
{
    int page_number;
    int32 total_bitmaps_passed = 0;
    for (page_number = 0; page_number < n; page_number++)
    {
        mdjvu_image_t page = pages[page_number];
        int32 bitmap_count = mdjvu_image_get_bitmap_count(page);
        int32 i; /* index of bitmap in a page */

        if (!mdjvu_image_has_substitutions(page))
            mdjvu_image_enable_substitutions(page);

        for (i = 0; i < bitmap_count; i++)
        {
            int32 tag = tags[total_bitmaps_passed++];
            if (!tag) continue; /* skip non-substitutable bitmaps */
            mdjvu_image_set_substitution(page,
                                         mdjvu_image_get_bitmap(page, i),
                                         representatives[tag]);
        }
    }
}

/* -------------------------------------------------------------------------- */

static void report_classify(void *param, int page_completed)
{
    mdjvu_compression_options_t r = (mdjvu_compression_options_t) param;
    if (r->report)
    {
        printf(_("Classification: %d of %d completed\n"), r->report_start_page + page_completed,
                                                       r->report_total_pages);
    }
}

static void report_prototypes(void *param, int page_completed)
{
    mdjvu_compression_options_t r = (mdjvu_compression_options_t) param;
    if (r->report)
    {
        printf(_("Prototype search: %d of %d completed\n"), r->report_start_page + page_completed,
                                                         r->report_total_pages);
    }
}


MDJVU_FUNCTION mdjvu_image_t mdjvu_compress_multipage(int n, mdjvu_image_t *pages, mdjvu_compression_options_t options)
{
    mdjvu_image_t dictionary = NULL;
    int i;
    int32 total_bitmaps_count, max_tag;
    mdjvu_bitmap_t *representatives;
    int32 *tags;
    int32 *npatterns;
    unsigned char *dictionary_flags;

    total_bitmaps_count = 0;
    for (i = 0; i < n; i++)
    {
        total_bitmaps_count += mdjvu_image_get_bitmap_count(pages[i]);

        if (options->verbose) printf(_("deciding what pieces are letters in page #%d\n"), i);
        mdjvu_calculate_not_a_letter_flags(pages[i]);

        if (options->verbose) printf(_("sorting letters in page #%d\n"), i);
        mdjvu_sort_blits(pages[i]);
        mdjvu_image_sort_bitmaps(pages[i]);

        if (!mdjvu_image_has_substitutions(pages[i]))
            mdjvu_image_enable_substitutions(pages[i]);
    }

    tags = MDJVU_MALLOCV(int32, total_bitmaps_count);
    if (options->report) printf(_("started classification\n"));
    max_tag = mdjvu_multipage_classify_bitmaps
        (n, total_bitmaps_count, pages, tags,
         ((struct MinidjvuCompressionOptions *) options)->matcher_options,
         report_classify, options, options->averaging);
    if (options->report) printf(_("finished classification\n"));

    dictionary_flags = (unsigned char *) malloc((max_tag + 1));
    representatives = (mdjvu_bitmap_t *)
        malloc((max_tag + 1) * sizeof(mdjvu_bitmap_t));

    npatterns = MDJVU_MALLOCV(int32, n);
    for (i = 0; i < n; i++)
    {
        npatterns[i] = mdjvu_image_get_bitmap_count(pages[i]);
    }
    mdjvu_multipage_get_dictionary_flags(n, npatterns,
                               max_tag, tags, dictionary_flags);
    MDJVU_FREEV(npatterns);

    if (!options->averaging)
    {
        mdjvu_multipage_choose_representatives(n, pages, max_tag, tags, representatives, dictionary_flags);
        dictionary = get_dictionary(max_tag, representatives);

    }
    else
        dictionary = mdjvu_multipage_choose_average_representatives(
            n, pages, total_bitmaps_count, max_tag, tags, representatives, dictionary_flags);

    for (i = 0; i < n; i++)
        mdjvu_image_set_dictionary(pages[i], dictionary);

    set_substitutions(n, pages, total_bitmaps_count, tags, representatives);

    mdjvu_multipage_adjust(dictionary, n, pages);
    for (i = 0; i < n; i++)
        mdjvu_image_remove_unused_bitmaps(pages[i]);

    if (options->no_prototypes) {
        for (i = 0; i < n; i++) {
            // it seems this required to suppress prototype search for local sjbz later
            mdjvu_image_enable_prototypes(pages[i]);
        }
    } else {
        if (options->report) printf(_("started prototype search\n"));
        mdjvu_multipage_find_prototypes(dictionary, n, pages,
                                        report_prototypes, options);
        if (options->report) printf(_("finished prototype search\n"));
    }

    free(dictionary_flags);
    free(representatives);
    MDJVU_FREEV(tags);

    return dictionary;
}
