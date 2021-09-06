/*
 * classify.h - classifying patterns
 */


#ifndef MDJVU_CLASSIFY_H
#define MDJVU_CLASSIFY_H

/* Classifies a set of patterns.
 * result - array of tags ranging from 1 to return value,
 *    and 0 for those cells which were NULL (yes, NULLs are permitted).
 * Every tag has at least one pattern to which it's attached.
 * Equally tagged images are classified equivalent.
 */

MDJVU_FUNCTION int32 mdjvu_classify_patterns
    (mdjvu_pattern_t *, int32 *result, int32 n, int32 dpi,
     mdjvu_matcher_options_t, int verbose);

#ifndef NO_MINIDJVU /* that's for DjVuLibre */

/* Special tag 0 is reserved for bitmaps marked "no-substitution".
 * If centers_needed, also extract bitmap centers from the patterns.
 */
MDJVU_FUNCTION int32 mdjvu_classify_bitmaps
    (mdjvu_image_t, int32 *result, mdjvu_matcher_options_t, int centers_needed, int verbose);



/* MULTIPAGE CLASSIFICATION */

/* npages - number of pages
 * total_npatterns - total number of patterns (must be sum over npatterns)
 * npatterns[i] - number of patterns on the i-th page
 * dpi[i] - resolution of the i-th page
 * result[i] - i-th tag; tags from all pages are put consecutively
 *
 * return value - maximal tag
 *
 * XXX: dpi is not correctly handled
 */

MDJVU_FUNCTION int32 mdjvu_multipage_classify_patterns
	(int32 npages, int32 total_npatterns, const int32 *npatterns, mdjvu_pattern_t **,
	 int32 *result, const int32 *dpi, mdjvu_matcher_options_t,
     void (*report)(void *, int), int verbose);

MDJVU_FUNCTION int32 mdjvu_multipage_classify_bitmaps
    (int32 npages, int32 total_npatterns, mdjvu_image_t *,
     int32 *result, mdjvu_matcher_options_t,
     void (*report)(void *, int), int centers_needed, int verbose);


/* Decide what bitmaps will be put into the dictionary (by tag).
 * This implementation simply chooses tags which occur more than in one page.
 *
 * Arguments:
 *      dictionary_flags - [0..max_tag] array of 1/0 flags; it's the result
 */

MDJVU_FUNCTION void mdjvu_multipage_get_dictionary_flags
	(int32 npages, const int32 *npatterns, int32 max_tag,
	 const int32 *tags, unsigned char *dictionary_flags);



#endif /* NO_MINIDJVU */

#endif /* MDJVU_CLASSIFY_H */
