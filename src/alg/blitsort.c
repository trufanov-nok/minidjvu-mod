/*
 * blitsort.c - sorting blits in an image
 */

/* The algorithm is taken from DjVuLibre with minor (stupid) changes. */


#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h>



typedef struct
{
    int32 left, top, right, bottom;
    int32 original_index;
} BlitPassport;


static int compare_top_edges_downward(const void *p1, const void *p2)
{
    const BlitPassport *bp1 = (const BlitPassport *) p1;
    const BlitPassport *bp2 = (const BlitPassport *) p2;
    int32 d;
    d = bp1->top - bp2->top;
    if (d) return d;
    d = bp1->left - bp2->left;
    if (d) return d;
    return 0;
}

static int compare_left_edges_rightward(const void *p1, const void *p2)
{
    const BlitPassport *bp1 = (const BlitPassport *) p1;
    const BlitPassport *bp2 = (const BlitPassport *) p2;
    int32 d;
    d = bp1->left - bp2->left;
    if (d) return d;
    d = bp1->top - bp2->top;
    if (d) return d;
    return 0;
}

static int compare_integers_reversed(const void *p1, const void *p2)
{
    if (sizeof(int) == sizeof(int32)) /* hopefully compilers will optimize it */
    {
        return * (const int *) p2 - * (const int *) p1;
    }
    else
    {
        int32 d = * (const int32 *) p2 - * (const int32 *) p1;
        return d > 0 ?  1 :
               d < 0 ? -1 : 0;
    }
}


MDJVU_IMPLEMENT void mdjvu_sort_blits(mdjvu_image_t img)
{
    /* We're going to sort only blits with `is_a_letter' flag set. */

    int32 char_blit_count = 0;
    int32 blit_count, i, j, maxtopchange, ccno;
    BlitPassport *bps;
    int32 *bottoms, *passport_of_blit;

    /* Count letter blits */
    blit_count = mdjvu_image_get_blit_count(img);
    for (i = 0; i < blit_count; i++)
    {
        mdjvu_bitmap_t bmp = mdjvu_image_get_blit_bitmap(img, i);
        if (!mdjvu_image_get_not_a_letter_flag(img, bmp))
            char_blit_count++;
    }

    if (char_blit_count < 2) return;

    /* Allocate `bps' and `bottoms' arrays */
    bps = (BlitPassport *) malloc(char_blit_count * sizeof(BlitPassport));
    bottoms = (int32 *) malloc(char_blit_count * sizeof(int32));

    /* Fill in `bps' with character blit passports */
    j = 0;
    for (i = 0; i < blit_count; i++)
    {
        mdjvu_bitmap_t bmp = mdjvu_image_get_blit_bitmap(img, i);
        if (!mdjvu_image_get_not_a_letter_flag(img, bmp))
        {
            int32 x = bps[j].left = mdjvu_image_get_blit_x(img, i);
            int32 y = bps[j].top  = mdjvu_image_get_blit_y(img, i);;
            bps[j].right  = x + mdjvu_bitmap_get_width(bmp)  - 1;
            bps[j].bottom = y + mdjvu_bitmap_get_height(bmp) - 1;
            bps[j].original_index = i;
            j++;
        }
    }

    /* Sort the BlitPassports list in top-to-bottom order. */
    qsort(bps, char_blit_count, sizeof(BlitPassport),
          &compare_top_edges_downward);

    /* Subdivide the ccarray list roughly into text lines [LYB] */
    /* Determine maximal top deviation */
    maxtopchange = mdjvu_image_get_width(img) / 40;
    if (maxtopchange < 32) maxtopchange = 32;

    /* Loop until processing all ccs */
    ccno = 0;
    while (ccno < char_blit_count)  /* ccno will be increasing constantly */
    {
        /* Gather first line approximation */
        int32 sublist_top    = bps[ccno].top;
        int32 sublist_bottom = bps[ccno].bottom;

        int32 nccno;

        /* nccno will be at least ccno + 1,
         * or otherwise we're hung.
         */
        for (nccno = ccno; nccno < char_blit_count; nccno++)
        {
            int32 bottom;
            if (bps[nccno].top > sublist_bottom) break;
            if (bps[nccno].top > sublist_top + maxtopchange) break;
            bottom = bps[nccno].bottom;
            bottoms[nccno - ccno] = bottom;
            if (bottom > sublist_bottom)
                sublist_bottom = bottom;
        }

        /* If more than one candidate cc for the line */
        if (nccno > ccno + 1)
        {
            /* Compute median bottom */
            int32 bottom;
            qsort(bottoms, nccno - ccno, sizeof(int32),
                  &compare_integers_reversed);
            bottom = bottoms[ (nccno - ccno - 1) / 2 ];

            /* Compose final line */
            for (nccno = ccno; nccno < char_blit_count; nccno++)
                if (bps[nccno].top > bottom)
                    break;

            /* Sort final line */
            qsort(bps + ccno, nccno - ccno, sizeof(BlitPassport),
                  &compare_left_edges_rightward);
        }

        /* Next line */
        ccno = nccno;
    }

    /* Permute the blits according to `bps' */
    passport_of_blit = (int32 *) malloc(blit_count * sizeof(int32));
    for (i = 0; i < blit_count; i++)
        passport_of_blit[i] = -1;
    for (i = 0; i < char_blit_count; i++)
        passport_of_blit[bps[i].original_index] = i;

    /* We'll maintain that bps[i].original_index points to the same blit */
    for (i = 0; i < char_blit_count; i++)
    {
        int32 blit_to_put_here = bps[i].original_index;
        mdjvu_image_exchange_blits(img, blit_to_put_here, i);
        if (passport_of_blit[i] != -1)
            bps[passport_of_blit[i]].original_index = blit_to_put_here;
        passport_of_blit[blit_to_put_here] = passport_of_blit[i];
    }

    free(passport_of_blit);
    free(bps);
    free(bottoms);
}
