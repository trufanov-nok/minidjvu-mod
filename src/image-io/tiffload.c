/*
 * tiffload.c - loading TIFF bitmaps
 */

#include "../base/mdjvucfg.h"

#ifdef HAVE_LIBTIFF
    #include <tiffio.h>
    #define MDJVU_USE_TIFFIO
#endif
#include <minidjvu-mod/minidjvu-mod.h>

#include <stdlib.h>
#include <string.h>

#ifdef HAVE_LIBTIFF

static mdjvu_bitmap_t load_tiff(const char *path, int32 *presolution, mdjvu_error_t *perr, uint32 idx)
{
    uint16 photometric;
    uint32 w, h;
    uint16 bits_per_sample = 0, samples_per_pixel = 0;
    float dpi;
    mdjvu_bitmap_t result;
    tsize_t scanline_size;
    unsigned char *scanline;
    uint32 i;

    TIFF *tiff = TIFFOpen(path, "r");
    for ( i=0; tiff != NULL && i<idx; i++ )
    {
        if (!TIFFReadDirectory(tiff))
            break;
    }

    *perr = NULL;
    if (!tiff || i<idx)
    {
        *perr = mdjvu_get_error(mdjvu_error_fopen_read);
        return NULL;
    }

    /* test if bitonal */
    TIFFGetFieldDefaulted(tiff, TIFFTAG_BITSPERSAMPLE, &bits_per_sample);
    TIFFGetFieldDefaulted(tiff, TIFFTAG_SAMPLESPERPIXEL, &samples_per_pixel);
    if (bits_per_sample != 1 || samples_per_pixel != 1)
    {
        *perr = mdjvu_get_error(mdjvu_error_corrupted_tiff);
        TIFFClose(tiff);
        return NULL;
    }

    /* photometric */
    photometric = PHOTOMETRIC_MINISWHITE;
    TIFFGetFieldDefaulted(tiff, TIFFTAG_PHOTOMETRIC, &photometric);

    /* image size */
    if (!TIFFGetFieldDefaulted(tiff, TIFFTAG_IMAGEWIDTH, &w)
     || !TIFFGetFieldDefaulted(tiff, TIFFTAG_IMAGELENGTH, &h))
    {
        *perr = mdjvu_get_error(mdjvu_error_corrupted_tiff);
        TIFFClose(tiff);
        return NULL;
    }

    /* get the resolution */
    if (presolution) {
        *presolution = -1;
        if (TIFFGetFieldDefaulted(tiff, TIFFTAG_XRESOLUTION, &dpi))
        {
            *presolution = (int32) dpi;
        }
    }

    result = mdjvu_bitmap_create(w, h);

    scanline_size = TIFFScanlineSize(tiff);

    if (scanline_size < mdjvu_bitmap_get_packed_row_size(result))
    {
        *perr = mdjvu_get_error(mdjvu_error_corrupted_tiff);
        TIFFClose(tiff);
        mdjvu_bitmap_destroy(result);
        return NULL;
    }

    scanline = (unsigned char *) malloc(scanline_size);

    for (i = 0; i < h; i++)
    {
        if (TIFFReadScanline(tiff, (tdata_t)scanline, i, 0) < 0)
        {
            *perr = mdjvu_get_error(mdjvu_error_corrupted_tiff);
            TIFFClose(tiff);
            free(scanline);
            mdjvu_bitmap_destroy(result);
            return NULL;
        }

        if (photometric != PHOTOMETRIC_MINISWHITE)
        {
            /* invert the row */
            int32 k;
            int32 s = (int32) scanline_size;
            for (k = 0; k < s; k++)
                scanline[k] = ~scanline[k];
        }

        /* clear the padding bits */
        if (w & 7)
            scanline[scanline_size - 1] &= ~(0xFF >> (w & 7));

        memcpy(mdjvu_bitmap_access_packed_row(result, i),
               scanline,
               mdjvu_bitmap_get_packed_row_size(result));
    }

    free(scanline);

    TIFFClose(tiff);
    return result;
}

MDJVU_IMPLEMENT uint32 mdjvu_get_tiff_page_count(const char *path)
{
    int dircount = 0;
    TIFF* tif = TIFFOpen(path, "r");

    /* a "directory" is a page in a multipage tiff */

    if ( tif ) {
        do {
            dircount++;
        } while ( TIFFReadDirectory(tif) );
        TIFFClose(tif);
    }
    return dircount;
}

#endif /* HAVE_LIBTIFF */

MDJVU_IMPLEMENT mdjvu_bitmap_t mdjvu_load_tiff(const char *path, int32 *presolution, mdjvu_error_t *perr, uint32 idx)
{
    #ifdef HAVE_LIBTIFF
        return load_tiff(path, presolution, perr, idx);
    #else
        *perr = mdjvu_get_error(mdjvu_error_tiff_support_disabled);
        return NULL;
    #endif
}
