/*
 * patterns.c - pattern matching algorithm
 */

/* This is `patterns.c', the unit that handles pattern matching.
 * Its task is only to compare pairs of images, not to classify a set of them.
 * And this has absolutely nothing to do with choosing a cross-coding prototype.
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include "bitmaps.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <endian.h>

#define TIMES_TO_THIN 1
#define TIMES_TO_THICKEN 1

#define SIGNATURE_SIZE 32


typedef struct
{
    double pithdiff1_threshold;
    double pithdiff2_threshold;
    double shiftdiff1_threshold;
    double shiftdiff2_threshold;
    double shiftdiff3_threshold;
    int aggression;
    int method;
} Options;

/* These are hand-tweaked parameters of this classifier. */

static const double pithdiff1_veto_threshold       = 23;
static const double pithdiff2_veto_threshold       = 4;
static const double shiftdiff1_veto_threshold      = 1000;
static const double shiftdiff2_veto_threshold      = 1500;
static const double shiftdiff3_veto_threshold      = 2000;

static const double size_difference_threshold = 10;
static const double mass_difference_threshold = 10;

static const double shiftdiff1_falloff        = .9;
static const double shiftdiff2_falloff        = 1;
static const double shiftdiff3_falloff        = 1.15;

static void interpolate(Options *opt, const double *v1, const double *v2,
                        int l, int r, int x)
{
    double w1 = ((double)(r - x)) / (r - l); /* weights */
    double w2 = 1 - w1;
    opt->pithdiff1_threshold  = v1[0] * w1 + v2[0] * w2;
    opt->pithdiff2_threshold  = v1[1] * w1 + v2[1] * w2;
    opt->shiftdiff1_threshold = v1[2] * w1 + v2[2] * w2;
    opt->shiftdiff2_threshold = v1[3] * w1 + v2[3] * w2;
    opt->shiftdiff3_threshold = v1[4] * w1 + v2[4] * w2;
}


/* Sets `aggression' for pattern matching.
 * Lower values are safer, bigger values produce smaller files.
 */

MDJVU_IMPLEMENT void mdjvu_set_aggression(mdjvu_matcher_options_t opt, int level)
{
    const double set200[5] = {30,     3, 200,  200, 15};
    const double set100[5] = {10,   0.9, 100,  100,  5};
    const double   set0[5] = {0,     0,   0,    0,   0};

    if (level < 0) level = 0;

    ((Options *) opt)->aggression = level;

    if (level > 100)
        interpolate((Options *) opt, set100, set200, 100, 200, level);
    else
        interpolate((Options *) opt, set0,   set100, 0,   100, level);
}

/* ========================================================================== */

MDJVU_IMPLEMENT mdjvu_matcher_options_t mdjvu_matcher_options_create(void)
{
    mdjvu_matcher_options_t options = (mdjvu_matcher_options_t) MALLOC1(Options);

    mdjvu_init();
    mdjvu_set_aggression(options, 100);
    ((Options *) options)->method = 0;
    return options;
}

MDJVU_IMPLEMENT void mdjvu_use_matcher_method(mdjvu_matcher_options_t opt, int method)
{
    ((Options *) opt)->method |= method;
}

MDJVU_IMPLEMENT void mdjvu_matcher_options_destroy(mdjvu_matcher_options_t opt)
{
    Options * options = (Options *) opt;
    FREE1(options);
}

/* ========================================================================== */

typedef unsigned char byte;

typedef struct ComparableImageData
{
    int32 lossless; // if set on the only meaningful field is bitmap
    mdjvu_bitmap_t bitmap; // NULL if not lossless
    byte **pixels; /* 0 - purely white, 255 - purely black (inverse to PGM!) */
    byte **pith2_inner;
    byte **pith2_outer;
//    byte **pith2_inner_old;
//    byte **pith2_outer_old;
    int32 width, height, mass;
    int32 mass_center_x, mass_center_y;
    byte signature[SIGNATURE_SIZE];  /* for shiftdiff 1 and 3 tests */
    byte signature2[SIGNATURE_SIZE]; /* for shiftdiff 2 test */
} Image;



/* Each image pair undergoes simple tests (dimensions and mass)
 * and at most five more advanced tests.
 * Each test may end up with three outcomes: veto (-1), doubt (0) and match(1).
 * Images are equivalent if and only if
 *     there was no `veto'
 *     and there was at least one `match'.
 */


/* We check whether images' dimensions are different
 *     no more than by size_difference_threshold percent.
 * Return value is usual: veto (-1) or doubt (0).
 * Mass checking was introduced by Leon Bottou.
 */

static int simple_tests(Image *i1, Image *i2)
{
    int32 w1 = i1->width, h1 = i1->height, m1 = i1->mass;
    int32 w2 = i2->width, h2 = i2->height, m2 = i2->mass;

    if (100.* w1 > (100.+ size_difference_threshold) * w2) return -1;
    if (100.* w2 > (100.+ size_difference_threshold) * w1) return -1;
    if (100.* h1 > (100.+ size_difference_threshold) * h2) return -1;
    if (100.* h2 > (100.+ size_difference_threshold) * h1) return -1;
    if (100.* m1 > (100.+ mass_difference_threshold) * m2) return -1;
    if (100.* m2 > (100.+ mass_difference_threshold) * m1) return -1;

    return 0;
}


#define USE_PITHDIFF 1
#define USE_SHIFTDIFF_1 1
#define USE_SHIFTDIFF_2 1
#define USE_SHIFTDIFF_3 1


/* Computing distance by comparing pixels {{{ */

/* This function compares two images pixel by pixel.
 * The exact way to compare pixels is defined by two functions,
 *     compare_row and compare_with_white.
 * Both functions take pointers to byte rows and their length.
 *
 * Now images are aligned by mass centers.
 * Code needs some clarification, yes...
 */
static int32 distance_by_pixeldiff_functions_by_shift(Image *i1, Image *i2,
    int32 (*compare_row)(byte *, byte *, int32),
    int32 (*compare_1_with_white)(byte *, int32),
    int32 (*compare_2_with_white)(byte *, int32),
    int32 ceiling,
    int32 shift_x, int32 shift_y) /* of i1's coordinate system with respect to i2 */
{
    int32 w1 = i1->width, w2 = i2->width, h1 = i1->height, h2 = i2->height;
    int32 min_y = shift_y < 0 ? shift_y : 0;
    int32 right1 = shift_x + w1;
    int32 max_y_plus_1 = h2 > shift_y + h1 ? h2 : shift_y + h1;
    int32 i;
    int32 min_overlap_x = shift_x > 0 ? shift_x : 0;
    int32 max_overlap_x_plus_1 = w2 < right1 ? w2 : right1;
    int32 min_overlap_x_for_i1 = min_overlap_x - shift_x;
    int32 max_overlap_x_plus_1_for_i1 = max_overlap_x_plus_1 - shift_x;
    int32 overlap_length = max_overlap_x_plus_1 - min_overlap_x;
    int32 score = 0;

    if (overlap_length <= 0) return INT32_MAX;

    for (i = min_y; i < max_y_plus_1; i++)
    {
        int32 y1 = i - shift_y;

        /* calculate difference in the i-th line */


        if (i < 0 || i >= h2)
        {
            /* calculate difference of i1 with white */
            score += compare_1_with_white(i1->pixels[y1], w1);
        }
        else if (i < shift_y || i >= shift_y + h1)
        {
            /* calculate difference of i2 with white */
            score += compare_2_with_white(i2->pixels[i], w2);
        }
        else
        {
            /* calculate difference in a line where the bitmaps overlap */
            score += compare_row(i1->pixels[y1] + min_overlap_x_for_i1,
                                 i2->pixels[i] + min_overlap_x,
                                 overlap_length);


            /* calculate penalty for the left margin */
            if (min_overlap_x > 0)
                score += compare_2_with_white(i2->pixels[i], min_overlap_x);
            else
                score += compare_1_with_white(i1->pixels[y1], min_overlap_x_for_i1);

            /* calculate penalty for the right margin */
            if (max_overlap_x_plus_1 < w2)
            {
                score += compare_2_with_white(
                    i2->pixels[i] + max_overlap_x_plus_1,
                    w2 - max_overlap_x_plus_1);
            }
            else
            {
                score += compare_1_with_white(
                     i1->pixels[y1] + max_overlap_x_plus_1_for_i1,
                     w1 - max_overlap_x_plus_1_for_i1);

            }
        }

        if (score >= ceiling) {
            return INT32_MAX;
        }
    }
    return score;
}

static int32 distance_by_pixeldiff_functions(Image *i1, Image *i2,
    int32 (*compare_row)(byte *, byte *, int32),
    int32 (*compare_1_with_white)(byte *, int32),
    int32 (*compare_2_with_white)(byte *, int32),
    int32 ceiling)
{
    byte **p1, **p2;
    int32 w1, w2, h1, h2;
    int32 shift_x, shift_y; /* of i1's coordinate system with respect to i2 */
    /*int32 s = 0, i, i_start, i_cap;
    int32 right_margin_start, right_margin_width;*/

    /* make i1 to be narrower than i2 */
    if (i1->width > i2->width)
    {
        Image *img = i1;
        i1 = i2;
        i2 = img;
    }

    w1 = i1->width; h1 = i1->height; p1 = i1->pixels;
    w2 = i2->width; h2 = i2->height; p2 = i2->pixels;

    /* (shift_x, shift_y) */
    /*     is what should be added to i1's coordinates to get i2's coordinates. */
    shift_x = (w2 - w2/2) - (w1 - w1/2); /* center favors right */
    shift_y = h2/2 - h1/2;               /* center favors top */

    shift_x = i2->mass_center_x - i1->mass_center_x;
    if (shift_x < 0)
        shift_x = (shift_x - MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;
    else
        shift_x = (shift_x + MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;

    shift_y = i2->mass_center_y - i1->mass_center_y;
    if (shift_y < 0)
        shift_y = (shift_y - MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;
    else
        shift_y = (shift_y + MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;

    return distance_by_pixeldiff_functions_by_shift(
        i1, i2, compare_row, compare_1_with_white, compare_2_with_white,
        ceiling, shift_x, shift_y);
}

/* Computing distance by comparing pixels }}} */
/* inscribed framework penalty counting {{{ */

/* (Look at `frames.c' to see what it's all about) */

#if USE_PITHDIFF

/* If the framework of one letter is inscribed into another and vice versa,
 *     then those letters are probably equivalent.
 * That's the idea...
 * Counting penalty points here for any pixel
 *     that's framework in one image and white in the other.
 */

static int32 pithdiff_compare_row(byte *row1, byte *row2, int32 n)
{
    int32 i, s = 0;
    for (i = 0; i < n; i++) {
        if (row1[i] == 0xFF || row2[i] == 0xFF)
            s += fabs(row1[i] - row2[i]);
    }
    return s;
}

static int32 pithdiff_compare_with_white(byte *row, int32 n)
{
    int32 i, s = 0;
    for (i = 0; i < n; i++) if (row[i] == 255) s += 255;
    return s;
}

static int32 pithdiff_distance(Image *i1, Image *i2, int32 ceiling)
{
    return distance_by_pixeldiff_functions(i1, i2,
            &pithdiff_compare_row,
            &pithdiff_compare_with_white,
            &pithdiff_compare_with_white,
            ceiling);
}

static int pithdiff_equivalence(Image *i1, Image *i2, double threshold, int32 dpi)
{
    int32 perimeter = i1->width + i1->height + i2->width + i2->height;
    int32 ceiling = (int32) (pithdiff1_veto_threshold * dpi * perimeter / 100);
    int32 d = pithdiff_distance(i1, i2, ceiling);
    if (d == INT32_MAX) return -1;
    if (d < threshold * dpi * perimeter / 100) return 1;
    return 0;
}

#endif /* if USE_PITHDIFF */

/* inscribed framework penalty counting }}} */

/* shift signature comparison {{{ */

/* Just finding the square of a normal Euclidean distance between vectors
 * (but with falloff)
 */

#if USE_SHIFTDIFF_1 || USE_SHIFTDIFF_2 || USE_SHIFTDIFF_3
static int shiftdiff_equivalence(byte *s1, byte *s2, double falloff, double veto, double threshold)
{
    int i, delay_before_falloff = 1, delay_counter = 1;
    double penalty = 0;
    double weight = 1;

    for (i = 1; i < SIGNATURE_SIZE; i++) /* kluge: ignores the first byte */
    {
        int difference = s1[i] - s2[i];
        penalty += difference * difference * weight;
        if (!--delay_counter)
        {
            weight *= falloff;
            delay_counter = delay_before_falloff <<= 1;
        }
    }

    if (penalty >= veto * SIGNATURE_SIZE) return -1;
    if (penalty <= threshold * SIGNATURE_SIZE) return 1;
    return 0;
}
#endif
/* shift signature comparison }}} */

/* Finding mass center {{{ */

static void get_mass_center(unsigned char **pixels, int32 w, int32 h,
                     int32 *pmass_center_x, int32 *pmass_center_y)
{
    double x_sum = 0, y_sum = 0, mass = 0;
    int32 i, j;

    for (i = 0; i < h; i++)
    {
        unsigned char *row = pixels[i];
        for (j = 0; j < w; j++)
        {
            unsigned char pixel = row[j];
            x_sum += pixel * j;
            y_sum += pixel * i;
            mass  += pixel;
        }
    }

    *pmass_center_x = (int32) (x_sum * MDJVU_CENTER_QUANT / mass);
    *pmass_center_y = (int32) (y_sum * MDJVU_CENTER_QUANT / mass);
}

/* Finding mass center }}} */


#if __BYTE_ORDER == __BIG_ENDIAN
inline size_t swap_t(size_t* val, int size) {
    if (size >= sizeof (size_t))
        return *val;

    size_t res = 0;
    memcpy(&res, val, size);
    return res;
}
#elif __BYTE_ORDER == __LITTLE_ENDIAN
inline size_t swap_t(size_t* val, unsigned int size) {
    size_t res_buf = 0;
    memcpy(&res_buf, val, size);

    unsigned char * a = (unsigned char *) &res_buf;
    for (unsigned int i = 0; i < (sizeof (size_t)/2); i++) {
        unsigned char t = a[i];
        a[i] = a[sizeof (size_t) - i - 1];
        a[sizeof (size_t) - i - 1] = t;
    }
    return res_buf;
}
#endif

static void sweep(unsigned char **pixels, unsigned char **source, int w, int h)
{
    /* pretty same implementation as in mdjvu_smooth() */

    const size_t int_len_in_bits = sizeof (size_t)*8;
    const size_t len = (w + (int_len_in_bits -1) ) / int_len_in_bits;
    const size_t tail_len = (w % int_len_in_bits) ? ((w % int_len_in_bits) + 7) >> 3 : sizeof (size_t);

    const size_t mask1 = (~(size_t)0x0) << 1; //0b11111..110
    const size_t mask2 = (size_t)0x01 << (int_len_in_bits-1); //0b100000.00
    const size_t mask3 = mask2 >> 1; //0b01000..00
    const size_t mask4 = mask3 >> 1; //0b00100..00
    const size_t mask5 = mask1 & ~mask2; //0b01111..10

    for (int y = 0; y < h; y++) {
        size_t *r_p = (size_t *) pixels[y]; /* result    */
        size_t *u_p = (y > 0) ? (size_t *) source[y-1] : NULL;
        size_t *t_p = (size_t *) source[y];
        size_t *l_p = (y+1 < h) ? (size_t *) source[y+1] : NULL;

        size_t u_buf = 0, t_buf = 0, l_buf = 0;
        size_t u_val = 0, t_val = 0, l_val = 0;
        size_t u_cur = 0, t_cur = 0, l_cur = 0;

        for (unsigned int i = 0; i < len; i++) {
            if (u_p) {
                u_cur = swap_t(u_p++, i==len-1?tail_len:sizeof (size_t));
                u_val = u_buf | (u_cur >> 2);
                u_buf = u_cur << (int_len_in_bits - 2);
            }
            if (l_p) {
                l_cur = swap_t(l_p++, i==len-1?tail_len:sizeof (size_t));
                l_val = l_buf | (l_cur >> 2);
                l_buf = l_cur << (int_len_in_bits - 2);
            }

            t_cur = swap_t(t_p++, i==len-1?tail_len:sizeof (size_t));
            t_val = t_buf | (t_cur >> 2);
            t_buf = t_cur << (int_len_in_bits - 2);

            size_t res = u_val | (t_val << 1) | t_val | (t_val >> 1) | l_val;

            size_t tail = res & mask3;
            size_t head = res & mask4;

            if (tail && i) {
                // for i == 0 tail is always false and this is not called
                // we access last byte instead of size_t* to not mess with int endiannes
                *(((unsigned char *)r_p)-1) |= 1; // last bit is always 0 bcs of mask5
            }


            res = u_cur | (t_cur << 1) | t_cur | (t_cur >> 1) | l_cur;

            if (i != len-1) {
                res &= mask5;
            }

            if (head) {
                res |= mask2;
            }

            res = swap_t(&res, /*i==len-1?-1*tail_len:*/sizeof (size_t));
            memcpy(r_p++, &res, (i==len-1)?tail_len:sizeof (size_t));
        }
    }

}

static unsigned char **quick_thin(unsigned char **pixels, int w, int h, int N)
{
    const int row_size = (w+7) >> 3;

    unsigned char **aux = mdjvu_create_2d_array(row_size, h);
    memcpy(aux[0], pixels[0], row_size*h);
    unsigned char **buf = mdjvu_create_2d_array(row_size, h);
    memset(buf[0], 0, row_size*h);

    invert_bitmap(aux, w, h);

    while (N--)
    {
        sweep(buf, aux, w, h);
        if (N) {
            assign_unpacked_bitmap(aux, buf, w, h);
        }
    }

    invert_bitmap(buf, w, h);

    mdjvu_destroy_2d_array(aux);
    return buf;
}

static unsigned char **quick_thicken(unsigned char **pixels, int w, int h, int N)
{
    int r_w = w + N * 2;
    int r_h = h + N * 2;

    const int row_size = (r_w+7) >> 3;

    unsigned char **aux = mdjvu_create_2d_array(row_size, r_h);
    memset(aux[0], 0, row_size*r_h);
    assign_unpacked_bitmap_with_shift(aux, pixels, w, h, N);
    unsigned char **buf = mdjvu_create_2d_array(row_size, r_h);
    memset(buf[0], 0, row_size*r_h);


    while (N--)
    {
        sweep(buf, aux, r_w, r_h);
        if (N) {
            assign_unpacked_bitmap(aux, buf, r_w, r_h);
        }
    }

    mdjvu_destroy_2d_array(aux);
    return buf;
}



static void sweep_old(unsigned char **pixels, unsigned char **source, int w, int h)
{
    int x, y;

    for (y = 0; y < h; y++)
    {
        unsigned char *row    = pixels[y];
        unsigned char *srow   = source[y];
        unsigned char *supper = source[y-1];
        unsigned char *slower = source[y+1];

        for (x = 0; x < w; x++)
        {
            row[x] = (  supper[x] |
                        srow[x-1] | srow[x] | srow[x+1] |
                        slower[x] );
        }
    }
}

static unsigned char **quick_thin_old(unsigned char **pixels, int w, int h, int N)
{
    unsigned char **aux = provide_margins(pixels, w, h, 1);
    unsigned char **buf = allocate_bitmap_with_white_margins(w, h);

    clear_bitmap(buf, w, h);
    invert_bitmap_old(aux, w, h, 0);

    while (N--)
    {
        sweep_old(buf, aux, w, h);
        assign_bitmap(aux, buf, w, h);
    }

    invert_bitmap_old(buf, w, h, 0);

    free_bitmap_with_margins(aux);
    return buf;
}

static unsigned char **quick_thicken_old(unsigned char **pixels, int w, int h, int N)
{
    int r_w = w + N * 2;
    int r_h = h + N * 2;
    int y, x;
    unsigned char **aux = allocate_bitmap_with_white_margins(r_w, r_h);
    unsigned char **buf = allocate_bitmap_with_white_margins(r_w, r_h);

    clear_bitmap(buf, r_w, r_h);
    clear_bitmap(aux, r_w, r_h);

    for (y = 0; y < h; y++) {
        memcpy(aux[y + N] + N, pixels[y], w);
    }

    while (N--)
    {
        sweep_old(buf, aux, r_w, r_h);
        if (N) {
            assign_bitmap(aux, buf, r_w, r_h);
        }
    }

    free_bitmap_with_margins(aux);
    return buf;
}


#ifndef NO_MINIDJVU
mdjvu_pattern_t mdjvu_pattern_create(mdjvu_matcher_options_t opt, mdjvu_bitmap_t bitmap, int32 enforce_lossless)
{
    mdjvu_init();

    Options *m_opt = (Options *) opt;
    Image *img = MALLOC1(Image);
    
    enforce_lossless |= !m_opt->aggression;
    img->lossless = enforce_lossless;
    img->bitmap = enforce_lossless ? bitmap : NULL;
    if (enforce_lossless) {
        img->pixels = img->pith2_inner = img->pith2_outer = NULL;
        img->mass = img->mass_center_x = img->mass_center_y = 0;
        return (mdjvu_pattern_t) img;
    }
    
    int32 w = mdjvu_bitmap_get_width(bitmap);
    int32 h = mdjvu_bitmap_get_height(bitmap);

    img->width = w;
    img->height = h;
    img->pixels = allocate_bitmap(w, h);
    mdjvu_bitmap_unpack_all(bitmap, img->pixels);
    img->mass = mdjvu_bitmap_get_mass(bitmap);

    mdjvu_soften_pattern(img->pixels, img->pixels, w, h);

    get_mass_center(img->pixels, w, h,
                    &img->mass_center_x, &img->mass_center_y);
    mdjvu_get_gray_signature(img->pixels, w, h,
                             img->signature, SIGNATURE_SIZE);

    mdjvu_get_black_and_white_signature(img->pixels, w, h,
                                        img->signature2, SIGNATURE_SIZE);

    //  the !m_opt->aggression is interpreted as lossless now
    // if (!m_opt->aggression)
    // {
    //     free_bitmap(img->pixels);
    //     img->pixels = NULL;
    // }

    if (m_opt->method & MDJVU_MATCHER_PITH_2)
    {
        img->pith2_inner = quick_thin( mdjvu_bitmap_access_packed_data(bitmap), w, h, TIMES_TO_THIN);
        img->pith2_outer = quick_thicken( mdjvu_bitmap_access_packed_data(bitmap), w, h, TIMES_TO_THICKEN);

//        byte **pixels = mdjvu_create_2d_array(w, h);
//        mdjvu_bitmap_unpack_all(bitmap, pixels);
//        img->pith2_inner_old = quick_thin_old(pixels, w, h, TIMES_TO_THIN); //quick_thin_old(pixels, w, h, TIMES_TO_THIN);
//        img->pith2_outer_old = quick_thicken_old(pixels, w, h, TIMES_TO_THICKEN);
//        mdjvu_destroy_2d_array(pixels);
//        byte **pixels = mdjvu_create_2d_array(w, h);
//        mdjvu_bitmap_unpack_all(bitmap, pixels);
//        byte ** temp =  quick_thin_old(pixels, w, h, TIMES_TO_THIN); //quick_thin_old(pixels, w, h, TIMES_TO_THIN);
//        mdjvu_destroy_2d_array(pixels);
//        mdjvu_bitmap_t tt = mdjvu_bitmap_create(w/*+TIMES_TO_THICKEN*2*/, h/*+TIMES_TO_THICKEN*2*/);
//        mdjvu_bitmap_pack_all(tt, temp);
//        free_bitmap_with_margins(temp);

//        // comp
//        int rr = 0;


//        for (int y = 0; y < h/*+TIMES_TO_THICKEN*2*/; y++) {
//            unsigned char* c1 = mdjvu_bitmap_access_packed_row(tt, y);
//            unsigned char* c2 = img->pith2_inner[y];
//            for (int x = 0; x < (w/*+TIMES_TO_THICKEN*2*/+7)>>3; x++) {
//                if (c1[x] != c2[x]) {
//                    fprintf(stderr, "err %u %u: %u %u %u\n",c1[x],c2[x], y, x, w+TIMES_TO_THICKEN*2);
//                    int sadfa = pow(1,43);
//                    rr = sadfa;
//                }
//            }

//        }


//        mdjvu_bitmap_destroy(tt);

        assert(img->pith2_inner);
        assert(img->pith2_outer);
    }
    else
    {
        img->pith2_inner = NULL;
        img->pith2_outer = NULL;
    }

    return (mdjvu_pattern_t) img;
}
#endif


/* get a center (in 1/MDJVU_CENTER_QUANT pixels; defined in the header for image) */
MDJVU_IMPLEMENT void mdjvu_pattern_get_center(mdjvu_pattern_t p, int32 *cx, int32 *cy)
{
    *cx = ((Image *) p)->mass_center_x;
    *cy = ((Image *) p)->mass_center_y;
}


// Generate a lookup table for 8 bit integers
#define B2(n) n, n + 1, n + 1, n + 2
#define B4(n) B2(n), B2(n + 1), B2(n + 1), B2(n + 2)
#define B6(n) B4(n), B4(n + 1), B4(n + 1), B4(n + 2)

// Lookup table that store the sum of bits for all uchar values

const unsigned char lookup_table[256] = { B6(0), B6(1), B6(1), B6(2) };


inline int32 pith2_calc_val(size_t val) {
    unsigned char* v = (unsigned char*) &val;
    int32 s = 0;
    for (int i = 0; i < sizeof(size_t); i++)
        s += lookup_table[v[i]];
    return s;
}

inline size_t pith2_row_subset_op(size_t val_a, size_t val_b, char inverted) {
    return (inverted)    ?    val_b & ~val_a    :    val_a & ~val_b;
}

static int32 pith2_row_subset(byte *A, int32 pos_a, byte *B, int32 pos_b, int32 w)
{
    A += pos_a / 8; pos_a %= 8;
    B += pos_b / 8; pos_b %= 8;

    char inv = 0;
    if (pos_a < pos_b) {
        byte* t = A; A = B; B = t;
        int32 p = pos_a; pos_a = pos_b; pos_b = p;
        inv = 1;
    }

    const int size_t_len_bits = sizeof(size_t)*8;
    const unsigned char shift_right = pos_a - pos_b; // shift_right is < 8

    const size_t mask = ~(size_t)0;
    const size_t start_mask = mask >> pos_a;
    const size_t end_mask   = mask << (size_t_len_bits - ( (pos_a + w) % size_t_len_bits)) % size_t_len_bits;

    int len_a = (pos_a + w) / size_t_len_bits;
    int len_b = (pos_b + w) / size_t_len_bits;
    int tail_a = (((pos_a + w) % size_t_len_bits) + 7) >> 3;
    int tail_b = (((pos_b + w) % size_t_len_bits) + 7) >> 3;
    if (!tail_a) {
        len_a--; tail_a = sizeof(size_t);
    }
    if (!tail_b) {
        len_b--; tail_b = sizeof(size_t);
    }

    size_t * a = (size_t *) A;
    size_t * b = (size_t *) B;

    size_t val_a, val_b, buf = 0;

    int first = 1;
    int32 s = 0;
    while (len_a >= 0) { // len_a >= len_b
        if (len_a > 0) {
            val_a = swap_t(a++, sizeof(size_t));
        } else if (len_a == 0 && tail_a) {
            val_a = swap_t(a, tail_a);
        } else val_a = 0;

        if (len_b > 0) {
            val_b = swap_t(b++, sizeof(size_t));
        } else if (len_b == 0 && tail_b) {
            val_b = swap_t(b, tail_b);
        } else val_b = 0;

        if (shift_right) {
            size_t t = val_b << (size_t_len_bits - shift_right);
            val_b = buf | (val_b >> shift_right);
            buf = t;
        }

        size_t val =  pith2_row_subset_op(val_a, val_b, inv);

        if (first) {
            val &= start_mask;
            first = 0;
        }
        if (len_a == 0) {
            val &= end_mask;
        }

        len_a--; len_b--;

        s += pith2_calc_val( val );
    }

//    if (!(len_a+len_b)) {
//        val_a = swap_t(a, tail_a);
//        val_b = swap_t(b, tail_b);
//        val_a &= start_mask_a & end_mask_a;
//        val_b &= start_mask_b & end_mask_b;
//        val_b >>= shift_right;
//        val_a = pith2_row_subset_op(val_a, val_b, inv);
//        return pith2_calc_val(val_a) * 255;
//    }

//    val_a = *a++ & start_mask_a;
//    val_b = *b++ & start_mask_b;

//    val_a = swap_t(&val_a, sizeof(size_t));
//    val_b = swap_t(&val_b, sizeof(size_t));

//    buf = val_b << (sizeof(size_t) - shift_right);
//    val_b >>= shift_right;

//    int32 s = pith2_calc_val( pith2_row_subset_op(val_a, val_b, inv) );

//    for (int i = 1; i < len -1; i++) {
//         val_a = swap_t(a++, sizeof(size_t));
//         val_b = swap_t(b++, sizeof(size_t));
//         size_t t = val_b << (sizeof(size_t) - shift_right);
//         val_b = buf | (val_b >> shift_right);
//         buf = t;
//         s += pith2_calc_val( pith2_row_subset_op(val_a, val_b, inv) );
//    }

//    val_a = val_b = 0;
//    if (tail_a)
//    memcpy(&val_a, a, tail_a);
//    if (tail_b)
//    memcpy(&val_b, b, tail_b);
//    val_a &= end_mask_a;
//    val_b &= end_mask_b;
//    val_a = swap_t(&val_a, sizeof(size_t));
//    val_b = swap_t(&val_b, sizeof(size_t));
//    val_b = buf | (val_b >> shift_right);
//    val_a = pith2_row_subset_op(val_a, val_b, inv);
//    s += pith2_calc_val(val_a);

    return s * 255;
}

static int32 pith2_row_has_black(byte *row, int32 start_idx, int32 length)
{
    if (!length) return 0;

    row += start_idx / 8;
    start_idx %= 8;
    const unsigned char start_mask = 0xFF >> start_idx;
    const unsigned char end_mask   =  0xFF << ((8 - ( (start_idx + length) % 8)) % 8);
    const int32 len = ((start_idx + length + 7) >> 3) - 1;

    if (len) {
        int32 i, s = lookup_table[row[0] & start_mask];
        for (i = 1; i < len; i++) {
            s += lookup_table[row[i]];
        }
        s += lookup_table[row[len] & end_mask];
        return s * 255;
    }
    return lookup_table[*row & start_mask & end_mask] * 255;
}

static int32 pith2_row_subset_old(byte *A, byte *B, int32 length)
{
    int32 i, s = 0;
    for (i = 0; i < length; i++)
    {
        if (A[i] & !B[i])
            s += 255;
    }
    return s;
}

static int32 pith2_row_has_black_old(byte *row, int32 length)
{
    int32 i, s = 0;
    for (i = 0; i < length; i++)
    {
        if (row[i])
            s += 255;
    }
    return s;
}

static int32 pith2_return_0(byte *A, int32 length)
{
    return 0;
}


static int pith2_is_subset(mdjvu_pattern_t ptr1, mdjvu_pattern_t ptr2, double threshold, int32 dpi)
{
    Image *img1 = (Image *) ptr1;
    Image *img2 = (Image *) ptr2;
    Image ptr1_inner;
    Image ptr2_outer;
    int32 perimeter = img1->width + img1->height + img2->width + img2->height;
    int32 ceiling = (int32) (pithdiff2_veto_threshold * dpi * perimeter / 100);
    int32 d = 0;

    ptr1_inner.pixels = img1->pith2_inner;
    assert(img1->pith2_inner);
    ptr1_inner.width  = img1->width;
    ptr1_inner.height = img1->height;
    ptr1_inner.mass_center_x = img1->mass_center_x;
    ptr1_inner.mass_center_y = img1->mass_center_y;

    ptr2_outer.pixels = img2->pith2_outer;
    assert(img2->pith2_outer);
    ptr2_outer.width  = img2->width  + TIMES_TO_THICKEN*2;
    ptr2_outer.height = img2->height + TIMES_TO_THICKEN*2;
    ptr2_outer.mass_center_x = img2->mass_center_x + MDJVU_CENTER_QUANT;
    ptr2_outer.mass_center_y = img2->mass_center_y + MDJVU_CENTER_QUANT;


    Image *i1 = &ptr1_inner;
    Image *i2 = &ptr2_outer;

//    int32 score2 = distance_by_pixeldiff_functions(&ptr1_inner, &ptr2_outer,
//        &pith2_row_subset_old,
//        &pith2_row_has_black_old,
//        &pith2_return_0,
//        ceiling);


//    ptr1_inner.pixels = img1->pith2_inner;
//    ptr2_outer.pixels = img2->pith2_outer;

    int32 shift_x, shift_y; /* of i1's coordinate system with respect to i2 */
    int32 w1, w2, h1, h2;
    {
        /* make i1 to be narrower than i2 */
        if (i1->width > i2->width)
        {
            Image* img = i1;
            i1 = i2;
            i2 = img;
        }

        w1 = i1->width; h1 = i1->height;
        w2 = i2->width; h2 = i2->height;

        /* (shift_x, shift_y) */
        /*     is what should be added to i1's coordinates to get i2's coordinates. */
        shift_x = (w2 - w2/2) - (w1 - w1/2); /* center favors right */
        shift_y = h2/2 - h1/2;               /* center favors top */

        shift_x = i2->mass_center_x - i1->mass_center_x;
        if (shift_x < 0)
            shift_x = (shift_x - MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;
        else
            shift_x = (shift_x + MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;

        shift_y = i2->mass_center_y - i1->mass_center_y;
        if (shift_y < 0)
            shift_y = (shift_y - MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;
        else
            shift_y = (shift_y + MDJVU_CENTER_QUANT / 2) / MDJVU_CENTER_QUANT;
    }

//    distance_by_pixeldiff_functions_by_shift(
//            i1, i2, compare_row, compare_1_with_white, compare_2_with_white,
//            ceiling, shift_x, shift_y);

    int32 score = 0;

    {
        int32 min_y = shift_y < 0 ? shift_y : 0;
        int32 right1 = shift_x + w1;
        int32 max_y_plus_1 = h2 > shift_y + h1 ? h2 : shift_y + h1;
        int32 i;
        int32 min_overlap_x = shift_x > 0 ? shift_x : 0;
        int32 max_overlap_x_plus_1 = w2 < right1 ? w2 : right1;
        int32 min_overlap_x_for_i1 = min_overlap_x - shift_x;
        int32 max_overlap_x_plus_1_for_i1 = max_overlap_x_plus_1 - shift_x;
        int32 overlap_length = max_overlap_x_plus_1 - min_overlap_x;

        if (overlap_length <= 0) return -1;

        for (i = min_y; i < max_y_plus_1; i++)
        {
            int32 y1 = i - shift_y;

            /* calculate difference in the i-th line */

            if (i < 0 || i >= h2)
            {
                /* calculate difference of i1 with white */
                score += pith2_row_has_black(i1->pixels[y1], 0, w1);
            }
            else if (i >= shift_y && i < shift_y + h1)
            {
                /* calculate difference in a line where the bitmaps overlap */
                score += pith2_row_subset(i1->pixels[y1], min_overlap_x_for_i1,
                                          i2->pixels[i],  min_overlap_x,
                                          overlap_length);


                /* calculate penalty for the left margin */
                if (min_overlap_x <= 0) {
                    score += pith2_row_has_black(i1->pixels[y1], 0, min_overlap_x_for_i1);
                }

                /* calculate penalty for the right margin */
                if (max_overlap_x_plus_1 >= w2) {
                    score += pith2_row_has_black(
                                i1->pixels[y1], max_overlap_x_plus_1_for_i1,
                                w1 - max_overlap_x_plus_1_for_i1);

                }
            }

            if (score >= ceiling) return -1;
        }
    }


    if (score < threshold * dpi * perimeter / 100) return 1;
    return 0;
}

/* Requires `opt' to be non-NULL */
static int compare_patterns(mdjvu_pattern_t ptr1, mdjvu_pattern_t ptr2,/*{{{*/
                            int32 dpi, Options *opt)

{
    Image *i1 = (Image *) ptr1, *i2 = (Image *) ptr2;

    // check if lossless compression is enforced    
    if (i1->lossless != i2->lossless) {
        return -1;
    } else if (i1->lossless) {
        // just perform size check and memcmp()
        return mdjvu_bitmap_match(i1->bitmap, i2->bitmap) ? 1: -1;
    }

    // lossless is false
    
    
    int i, state = 0; /* 0 - unsure, 1 - equal unless veto */

    if (simple_tests(i1, i2)) return -1;

    #if USE_SHIFTDIFF_1
        i = shiftdiff_equivalence(i1->signature, i2->signature,
            shiftdiff1_falloff, shiftdiff1_veto_threshold, opt->shiftdiff1_threshold);
        if (i == -1) return -1;
        state |= i;
    #endif

    #if USE_SHIFTDIFF_2
        i = shiftdiff_equivalence(i1->signature2, i2->signature2,
            shiftdiff2_falloff, shiftdiff2_veto_threshold, opt->shiftdiff2_threshold);
        if (i == -1) return -1;
        state |= i;
    #endif

    #if USE_SHIFTDIFF_3
        i = shiftdiff_equivalence(i1->signature, i2->signature,
            shiftdiff3_falloff, shiftdiff3_veto_threshold, opt->shiftdiff3_threshold);
        if (i == -1) return -1;
        state |= i;
    #endif

    i = pith2_is_subset(ptr1, ptr2, opt->pithdiff2_threshold, dpi);
    if (i < 1) return i;
    i = pith2_is_subset(ptr2, ptr1, opt->pithdiff2_threshold, dpi);
    if (i < 1) return i;

    if (opt->method & MDJVU_MATCHER_RAMPAGE)
        return 1;

    #if USE_PITHDIFF
        if (opt->aggression > 0)
        {
            i = pithdiff_equivalence(i1, i2, opt->pithdiff1_threshold, dpi);
            if (i == -1) return 0; /* pithdiff has no right to veto at upper level */
            state |= i;
        }
    #endif

    #if 0
        if (opt->aggression > 0)
        {
            i = softdiff_equivalence(i1, i2, opt->softdiff_threshold, dpi);
            if (i == -1) return 0;  /* softdiff has no right to veto at upper level */
            state |= i;
        }
    #endif

    return state;
}/*}}}*/

MDJVU_IMPLEMENT int mdjvu_match_patterns(mdjvu_pattern_t ptr1, mdjvu_pattern_t ptr2,
                    int32 dpi, mdjvu_matcher_options_t options)
{
    Options *opt;
    int result;
    if (options)
        opt = (Options *) options;
    else
        opt = (Options *) mdjvu_matcher_options_create();

    result = compare_patterns(ptr1, ptr2, dpi, opt);

    if (!options)
        mdjvu_matcher_options_destroy((mdjvu_matcher_options_t) opt);

    return result;
}

MDJVU_IMPLEMENT int mdjvu_pattern_mem_size(mdjvu_pattern_t p)
{
   Image *img = (Image *) p;
   const int size_of_pointers_map = img->height * sizeof(unsigned char *);

   int res = sizeof(Image);
   if (img->pixels) res += img->width * img->height + size_of_pointers_map;

   const int row_size = (img->width + 7) >> 3;
   if (img->pith2_inner) res += row_size * img->height + size_of_pointers_map;
   if (img->pith2_outer) res += row_size * img->height + size_of_pointers_map;
   return res;
}

MDJVU_IMPLEMENT void mdjvu_pattern_destroy(mdjvu_pattern_t p)/*{{{*/
{
    Image *img = (Image *) p;
    if (img->pixels)
        free_bitmap(img->pixels);

    if (img->pith2_inner)
        mdjvu_destroy_2d_array(img->pith2_inner);

    if (img->pith2_outer)
        mdjvu_destroy_2d_array(img->pith2_outer);

//    if (img->pith2_inner_old)
//        free_bitmap_with_margins(img->pith2_inner_old);

//    if (img->pith2_outer_old)
//        free_bitmap_with_margins(img->pith2_outer_old);

    FREE1(img);
}/*}}}*/
