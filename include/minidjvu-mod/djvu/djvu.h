/*
 * djvu.h - functions to load from single-page DjVuBitonal files
 */

#ifndef MDJVU_DJVU_H
#define MDJVU_DJVU_H

#include "iff.h"

MDJVU_FUNCTION int mdjvu_locate_jb2_chunk(mdjvu_file_t file, int32 *plength, mdjvu_error_t *);

/*
 * BUG: resolution is not loaded
 */

MDJVU_FUNCTION mdjvu_image_t mdjvu_file_load_djvu_page(mdjvu_file_t file, mdjvu_error_t *);
MDJVU_FUNCTION mdjvu_image_t mdjvu_load_djvu_page(const char *path, mdjvu_error_t *);

/*
 * 1 - success, 0 - failure
 * After mdjvu_file_save_djvu_page() the file cursor is before the JB2 chunk.
 */
MDJVU_FUNCTION int mdjvu_files_save_djvu_dir(char **elements, int *sizes, int n,
                                             mdjvu_file_t file, mdjvu_file_t* tempfiles, int num_tempfiles, mdjvu_error_t *perr);
MDJVU_FUNCTION int mdjvu_file_save_djvu_dir( char **elements, int *sizes, int n,
                                             mdjvu_file_t file, mdjvu_file_t tmpfile, mdjvu_error_t *perr);
MDJVU_FUNCTION int mdjvu_save_djvu_dir(char **elements, int *sizes, int n, const char *path, mdjvu_error_t *perr);

MDJVU_FUNCTION int mdjvu_file_save_djvu_page(mdjvu_image_t, mdjvu_file_t, const char *dict_name,
                                             int insert_magic, mdjvu_error_t *perr, int erosion);
MDJVU_FUNCTION int mdjvu_save_djvu_page(mdjvu_image_t image, const char *path, const char *dict_name, mdjvu_error_t *perr, int erosion);

MDJVU_FUNCTION int mdjvu_file_save_djvu_dictionary(mdjvu_image_t, mdjvu_file_t,
                                             int insert_magic, mdjvu_error_t *, int erosion);
MDJVU_FUNCTION int mdjvu_save_djvu_dictionary(mdjvu_image_t image, const char *path, mdjvu_error_t *, int erosion);

MDJVU_FUNCTION void mdjvu_write_dirm_bundled(char **elements, int *sizes, int n, mdjvu_file_t f, mdjvu_error_t *perr);
MDJVU_FUNCTION void mdjvu_write_dirm_indirect(char **elements, int *sizes, int n, mdjvu_file_t f, mdjvu_error_t *perr);

/* Writes DjVu INFO chunk, as described in DjVu 2 Spec., 6.4.2, page 7.
 * Gamma and version numbers are always default.
 * This function does not write the chunk header.
 */
MDJVU_FUNCTION void mdjvu_write_info_chunk(mdjvu_file_t, mdjvu_image_t image);

/* Reads DjVu INFO chunk, as described in DjVu 2 Spec., 6.4.2, page 7.
 * Gamma and version numbers are dropped; width, height and dpi are returned.
 * Pointers may be NULL.
 * This function assumes the cursor is positioned to the start of chunk.
 */
MDJVU_FUNCTION void mdjvu_read_info_chunk(mdjvu_file_t,
                                          int32 *w, int32 *h, int32 *dpi);

#endif
