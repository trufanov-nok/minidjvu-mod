/*
 * delegate.c - choosing a representative over a class of letters
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

MDJVU_IMPLEMENT void mdjvu_multipage_choose_representatives
        (int32 npages,
         mdjvu_image_t *pages,
         int32 max_tag,
         int32 *tags,
         mdjvu_bitmap_t *representatives,
         unsigned char *dictionary_flags)
{
    int page_number;
    int32 total_bitmaps_passed = 0;
    memset(representatives, 0, (max_tag + 1) * sizeof(mdjvu_bitmap_t));
    for (page_number = 0; page_number < npages; page_number++)
    {
        mdjvu_image_t page = pages[page_number];
        int32 bitmap_count = mdjvu_image_get_bitmap_count(page);
        int32 i; /* index of bitmap in a page */

        for (i = 0; i < bitmap_count; i++)
        {
            int32 tag = tags[total_bitmaps_passed++];
            if (!tag || !dictionary_flags[tag]) continue; /* skip non-substitutable bitmaps */
            if (!representatives[tag])
                representatives[tag] = mdjvu_image_get_bitmap(page, i);
        }
    }
}

MDJVU_IMPLEMENT mdjvu_image_t mdjvu_multipage_choose_average_representatives
        (int32 npages,
         mdjvu_image_t *pages,
         int32 total_count,
         int32 max_tag,
         int32 *tags,
         mdjvu_bitmap_t *representatives,
         unsigned char *dictionary_flags)
{
    int32 page_number, tag;
    mdjvu_bitmap_t *sources;
    int32 *cx, *cy;
    mdjvu_image_t dictionary = mdjvu_image_create(0,0); /* 0 x 0 image */

    memset(representatives, 0, (max_tag + 1) * sizeof(mdjvu_bitmap_t));

    sources = (mdjvu_bitmap_t *) malloc(total_count * sizeof(mdjvu_bitmap_t));
    cx = (int32 *) malloc(total_count * sizeof(int32));
    cy = (int32 *) malloc(total_count * sizeof(int32));

    for (tag = 1; tag <= max_tag; tag++)
    {
        int32 sources_found = 0, total_bitmaps_passed = 0;
        mdjvu_bitmap_t rep;
        if (!dictionary_flags[tag]) continue;

        memset(sources, 0, total_count * sizeof(mdjvu_bitmap_t));
        memset(cx, 0, total_count * sizeof(int32));
        memset(cy, 0, total_count * sizeof(int32));
        for (page_number = 0; page_number < npages; page_number++)
        {
            mdjvu_image_t page = pages[page_number];
            int32 bitmap_count = mdjvu_image_get_bitmap_count(page);
            int32 i; /* index of bitmap in a page */
    
            for (i = 0; i < bitmap_count; i++)
            {
                if (tags[total_bitmaps_passed++] == tag)
                {
                    sources[sources_found] = mdjvu_image_get_bitmap(page, i);
                    mdjvu_image_get_center(page, sources[sources_found], &cx[sources_found], &cy[sources_found]);
                    sources_found++;
                }
            }
        }

        if (sources_found)
        {
            rep = mdjvu_average(sources, sources_found, cx, cy);
            representatives[tag] = rep;
        }
    }
    free(cx);
    free(cy);
    free(sources);

    for (page_number = 0; page_number < npages; page_number++)
        mdjvu_image_disable_centers(pages[page_number]);

    for (tag = 1; tag <= max_tag; tag++)
    {
        mdjvu_bitmap_t rep;
        
        if (!dictionary_flags[tag] || (rep = representatives[tag]) == NULL) continue;
        mdjvu_image_add_bitmap(dictionary, rep);
    }
    return dictionary;
}
