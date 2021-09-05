/*
 * matcher.h - matching patterns
 */

#ifndef MDJVU_MATCHER_H
#define MDJVU_MATCHER_H


typedef struct MinidjvuMatcherOptions *mdjvu_matcher_options_t;

MDJVU_FUNCTION mdjvu_matcher_options_t mdjvu_matcher_options_create(void);
MDJVU_FUNCTION void mdjvu_set_aggression(mdjvu_matcher_options_t, int level);

/* "matcher methods" (bitmask, not enum) */
#define MDJVU_MATCHER_DEFAULT 0
#define MDJVU_MATCHER_PITH_2  1
#define MDJVU_MATCHER_RAMPAGE 2

/* turn method on (|=) */
MDJVU_FUNCTION void mdjvu_use_matcher_method(mdjvu_matcher_options_t, int method);

MDJVU_FUNCTION void mdjvu_matcher_options_destroy(mdjvu_matcher_options_t);


/* To get an image ready for comparisons, one have to `prepare' it.
 * A prepared image is called a `pattern' here.
 * A pattern is mostly opaque except that its center may be retrieved.
 */


/* the struct itself is not defined in this header */
typedef struct MinidjvuPattern *mdjvu_pattern_t;


/* Allocate a pattern and calculate all necessary information.
 * Memory consumption is byte per pixel + constant (with default matcher).
 * The pattern would be independent on the bitmap given.
 *     (that is, you can destroy the bitmap immediately)
 */
#ifndef NO_MINIDJVU
MDJVU_FUNCTION mdjvu_pattern_t mdjvu_pattern_create(mdjvu_matcher_options_t, mdjvu_bitmap_t);
#endif

/* Return size of a pattern in memory in bytes */

MDJVU_FUNCTION int mdjvu_pattern_mem_size(mdjvu_pattern_t p);

/* Destroy the pattern. */

MDJVU_FUNCTION void mdjvu_pattern_destroy(mdjvu_pattern_t);


/* get a center (in 1/MDJVU_CENTER_QUANT pixels; defined in the header for image) */
MDJVU_FUNCTION void mdjvu_pattern_get_center(mdjvu_pattern_t, int32 *cx, int32 *cy);

/* Compare patterns.
 * Returns
 * +1 if images are considered equivalent,
 * -1 if they are considered totally different (just to speed up things),
 *  0 if unknown, but probably different.
 * Exchanging the order of arguments should not change the outcome.
 * If you have found that A ~ B and B ~ C,
 *     then you may assume A ~ C regardless of this function's result.
 *
 * Options may be NULL.
 */

MDJVU_FUNCTION int mdjvu_match_patterns(mdjvu_pattern_t, mdjvu_pattern_t,
                                        int32 dpi,
                                        mdjvu_matcher_options_t);


/* Auxiliary functions used in pattern matcher (TODO: comment them) */

/* `result' and `pixels' may be the same array */
MDJVU_FUNCTION void mdjvu_soften_pattern(unsigned char **result,
    unsigned char **pixels, int32 w, int32 h);

MDJVU_FUNCTION void mdjvu_get_gray_signature(
    unsigned char **data, int32 w, int32 h,
    unsigned char *result, int32 size);

MDJVU_FUNCTION void mdjvu_get_black_and_white_signature(
    unsigned char **data, int32 w, int32 h,
    unsigned char *result, int32 size);

#endif /* MDJVU_PATTERNS_H */
