/*
 * jb2save.cpp - saving JB2 raw streams
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include "jb2coder.h"
#include <stdlib.h>

static void encode_bitmap_with_erosion_support
    (JB2Encoder &jb2, mdjvu_bitmap_t bitmap, mdjvu_bitmap_t proto, int erosion)
{
    if (erosion)
    {
        mdjvu_bitmap_t flip_candidates = mdjvu_get_erosion_mask(bitmap);
        jb2.encode(bitmap, proto, flip_candidates);
        mdjvu_bitmap_destroy(flip_candidates);
    }
    else
        jb2.encode(bitmap, proto);
}

static int open_bitmap_record(mdjvu_image_t img, int32 index,
    bool with_blit, int32 *table, int32 &library_size, JB2Encoder &jb2,
    mdjvu_error_t *perr, int erosion)
{
    mdjvu_image_t dictionary = mdjvu_image_get_dictionary(img);
    int32 d = 0; /* shared dictionary size */
    if (dictionary) d = mdjvu_image_get_bitmap_count(dictionary);

    table[index] = -2;
    mdjvu_bitmap_t bitmap = mdjvu_image_get_bitmap(img, index);
    mdjvu_bitmap_t proto = mdjvu_image_get_prototype(img, bitmap);
    if (proto && mdjvu_image_has_bitmap(img, proto))
        proto = mdjvu_image_get_substitution(img, proto);

    if (!proto)
    {
        // encode directly

        if (with_blit)
            jb2.open_record(jb2_new_symbol_add_to_image_and_library);
        else
            jb2.open_record(jb2_new_symbol_add_to_library_only);

        encode_bitmap_with_erosion_support(jb2, bitmap, NULL, erosion);
    }
    else
    {
        if (mdjvu_image_has_bitmap(img, proto))
        {
            int32 proto_index = mdjvu_bitmap_get_index(proto);

            if (table[proto_index] == -2)
            {
                if (perr) *perr = mdjvu_get_error(mdjvu_error_recursive_prototypes);
                return 0; // prototype recursion error
            }

            if (table[proto_index] == -1)
            {
                // prototype is not yet encoded - encode it
                int32 bl, bt, bw, bh;
                mdjvu_bitmap_get_bounding_box(bitmap, &bl, &bt, &bw, &bh);
                assert(!bl);
                assert(!bt);
                assert(bw == mdjvu_bitmap_get_width(bitmap));
                assert(bh == mdjvu_bitmap_get_height(bitmap));

                int result = open_bitmap_record
                    (img, proto_index, false, table, library_size, jb2, perr, erosion);
                if (!result) return 0;
                jb2.close_record();
            }
        }

        if (with_blit)
            jb2.open_record(jb2_matched_symbol_with_refinement_add_to_image_and_library);
        else
            jb2.open_record(jb2_matched_symbol_with_refinement_add_to_library_only);
        jb2.matching_symbol_index.set_interval(0, d + library_size - 1);
        if (mdjvu_image_has_bitmap(img, proto))
        {
            int32 encoded_index = table[mdjvu_bitmap_get_index(proto)];
            jb2.zp.encode(d + encoded_index, jb2.matching_symbol_index);
        }
        else
        {
            assert(mdjvu_image_has_bitmap(dictionary, proto));
            jb2.zp.encode(mdjvu_image_get_dictionary_index(dictionary, proto),
                          jb2.matching_symbol_index);
        }

        encode_bitmap_with_erosion_support(jb2, bitmap, proto, erosion);
    }
    table[index] = library_size++;
    return 1;
}

MDJVU_IMPLEMENT int mdjvu_file_save_jb2_dictionary(mdjvu_image_t image, mdjvu_file_t f, mdjvu_error_t *perr, int erosion)
{
    if (!mdjvu_image_has_prototypes(image))
        mdjvu_find_prototypes(image);

    if (perr) *perr = NULL;
    int32 n = mdjvu_image_get_bitmap_count(image);
    JB2Encoder jb2((FILE *) f);
    ZPEncoder &zp = jb2.zp;

    /* opening record */
    jb2.open_record(jb2_start_of_image);
        zp.encode(0, jb2.image_size);
        zp.encode(0, jb2.image_size);
        zp.encode(0, jb2.eventual_image_refinement); /* seems to be unused */
    jb2.close_record();

    /* Now let's start. */
    int32 library_size = 0;

    // The library table keeps indices of shapes in the encoded library.
    // If the shape wasn't yet encoded, the value is -1.
    // If the shape is to be encoded, the value is -2.
    int32 *library_table = (int32 *) malloc(n * sizeof(int32));
    int32 i;
    for (i = 0; i < n; i++) library_table[i] = -1;

    for (i = 0; i < n; i++)
    {
        if (library_table[i] == -1)
        {
            // we have not encoded the bitmap yet
            // encode it!
            if (!open_bitmap_record(image, i, false, library_table,
                                                     library_size,
                                                     jb2, perr, erosion))
            {
                free(library_table);
                return 0;
            }
            jb2.close_record();
        }
    }

    /* end of data record */
    jb2.open_record(jb2_end_of_data);
    jb2.close_record();

    mdjvu_image_enable_dictionary_indices(image);
    for (i = 0; i < n; i++)
    {
        int32 dict_index = library_table[i];
        assert(dict_index >= 0);
        mdjvu_image_set_dictionary_index(image,
                                         mdjvu_image_get_bitmap(image, i),
                                         dict_index);
    }
    free(library_table);
    return 1;
}

MDJVU_IMPLEMENT int mdjvu_file_save_jb2(mdjvu_image_t image, mdjvu_file_t f, mdjvu_error_t *perr, int erosion)
{
    if (!mdjvu_image_has_prototypes(image))
        mdjvu_find_prototypes(image);

    if (perr) *perr = NULL;

    mdjvu_image_t dictionary = mdjvu_image_get_dictionary(image);

    int32 n = mdjvu_image_get_bitmap_count(image);
    int32 b = mdjvu_image_get_blit_count(image);
    JB2Encoder jb2((FILE *) f);
    ZPEncoder &zp = jb2.zp;

    int32 d = 0;
    if (dictionary && n > 0)
    {
        d = mdjvu_image_get_bitmap_count(dictionary);
        jb2.open_record(jb2_require_dictionary_or_reset);
            zp.encode(d, jb2.required_dictionary_size);
        jb2.close_record();
    }

    /* opening record */
    jb2.open_record(jb2_start_of_image);
        zp.encode(mdjvu_image_get_width(image),  jb2.image_size);
        zp.encode(mdjvu_image_get_height(image), jb2.image_size);
        zp.encode(0, jb2.eventual_image_refinement); /* seems to be unused */
    jb2.close_record();

    /* Now let's start. */
    int32 library_size = 0;

    // The library table keeps indices of shapes in the encoded library.
    // If the shape wasn't yet encoded, the value is -1.
    // If the shape is to be encoded, the value is -2.
    int32 *library_table = (int32 *) malloc(n * sizeof(int32));
    int32 i;
    for (i = 0; i < n; i++) library_table[i] = -1;

    /* encode all blits with bitmaps as necessary */
    for (i = 0; i < b; i++)
    {
        mdjvu_bitmap_t bitmap = mdjvu_image_get_blit_bitmap(image, i);
        if (!bitmap) continue;
        if (mdjvu_image_has_bitmap(image, bitmap))
            bitmap = mdjvu_image_get_substitution(image, bitmap);
        assert(bitmap);

        if (mdjvu_image_has_bitmap(image, bitmap)) // not the same check as before!
        {
            int32 bmp_i = mdjvu_bitmap_get_index(bitmap);

            if (library_table[bmp_i] == -1)
            {
                // we have not encoded the bitmap yet
                // encode it!
                if (!open_bitmap_record(image, bmp_i, /* with_blit: */ true, library_table,
                                                                       library_size,
                                                                       jb2, perr, erosion))
                {
                    free(library_table);
                    return 0;
                }
            }
            else
            {
                // we have our bitmap already
                assert(library_table[bmp_i] >= 0);
                jb2.open_record(jb2_matched_symbol_copy_to_image_without_refinement);
                jb2.matching_symbol_index.set_interval(0, d + library_size - 1);

                // encode the index
                jb2.zp.encode(d + library_table[bmp_i], jb2.matching_symbol_index);
            }
        }
        else
        {
            // the bitmap in question belongs to dictionary
            assert(d > 0);
            jb2.open_record(jb2_matched_symbol_copy_to_image_without_refinement);
            jb2.matching_symbol_index.set_interval(0, d + library_size - 1);
            assert(dictionary);
            assert(mdjvu_image_has_bitmap(dictionary, bitmap));
            assert(mdjvu_image_has_dictionary_indices(dictionary));
            jb2.zp.encode(mdjvu_image_get_dictionary_index(dictionary, bitmap),
                          jb2.matching_symbol_index);
        }

        jb2.encode_blit(image, i, mdjvu_bitmap_get_width(bitmap),
                        mdjvu_bitmap_get_height(bitmap));
        jb2.close_record();
    }

    /* end of data record */
    jb2.open_record(jb2_end_of_data);
    jb2.close_record();

    free(library_table);
    return 1;
}

MDJVU_IMPLEMENT int mdjvu_save_jb2(mdjvu_image_t image, const char *path, mdjvu_error_t *perr, int erosion)
{
    FILE *f = fopen(path, "wb");
    if (perr) *perr = NULL;
    if (!f)
    {
        if (perr) *perr = mdjvu_get_error(mdjvu_error_fopen_write);
        return 0;
    }
    int result = mdjvu_file_save_jb2(image, (mdjvu_file_t) f, perr, erosion);
    fclose(f);
    return result;
}

/* FIXME: code too close to mdjvu_save_jb2(); maybe, unite them? */
MDJVU_IMPLEMENT int mdjvu_save_jb2_dictionary(mdjvu_image_t image, const char *path, mdjvu_error_t *perr, int erosion)
{
    FILE *f = fopen(path, "wb");
    if (perr) *perr = NULL;
    if (!f)
    {
        if (perr) *perr = mdjvu_get_error(mdjvu_error_fopen_write);
        return 0;
    }
    int result = mdjvu_file_save_jb2_dictionary(image, (mdjvu_file_t) f, perr, erosion);
    fclose(f);
    return result;
}

