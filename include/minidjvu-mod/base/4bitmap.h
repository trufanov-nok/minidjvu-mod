/*
 * 4bitmap.h - routines for handling packed bitmaps
 */

typedef struct MinidjvuBitmap *mdjvu_bitmap_t;


#ifndef NDEBUG
extern int32 alive_bitmap_counter;
#endif


/* Create a bitmap.
 * Width and height must be positive, or NULL is returned.
 */
MDJVU_FUNCTION mdjvu_bitmap_t mdjvu_bitmap_create(int32 width, int32 height);

/* Destroy a bitmap. Each created bitmap must be destroyed sometime. */
MDJVU_FUNCTION void mdjvu_bitmap_destroy(mdjvu_bitmap_t);

/* Return size of a bitmap in memory in bytes. */
MDJVU_FUNCTION int mdjvu_bitmap_mem_size(mdjvu_bitmap_t bmp);

/* Get the width and height of a bitmap. */
MDJVU_FUNCTION int32 mdjvu_bitmap_get_width(mdjvu_bitmap_t);
MDJVU_FUNCTION int32 mdjvu_bitmap_get_height(mdjvu_bitmap_t);

/* Each bitmap keeps its index.
 * If this bitmap is not attached to an image, the index will be -1.
 */
MDJVU_FUNCTION int32 mdjvu_bitmap_get_index(mdjvu_bitmap_t);
MDJVU_FUNCTION void mdjvu_bitmap_set_index(mdjvu_bitmap_t, int32 new_value);

/* Returns the size of a packed row in bytes.
 * Packing stores 8 pixels to a byte,
 *     so the answer is (width + 7) / 8.
 */
MDJVU_FUNCTION int32 mdjvu_bitmap_get_packed_row_size(mdjvu_bitmap_t);


/* Get a pointer to the bitmap's packed row. Use with caution.
 * Packing is PBM-ish:
 *     most significant bit is the leftmost one,
 *     bytes go left to right.
 */
MDJVU_FUNCTION unsigned char *
    mdjvu_bitmap_access_packed_row(mdjvu_bitmap_t, int32);

MDJVU_FUNCTION unsigned char **
    mdjvu_bitmap_access_packed_data(mdjvu_bitmap_t b);

/* Fill a given row by the shape's row with the given Y coordinate.
 * The coordinate varies from 0 (top) to height-1 (bottom).
 * The memory should be enough to write <width> bytes.
 * White is 0, black is an undefined nonzero value.
 */
MDJVU_FUNCTION void
    mdjvu_bitmap_unpack_row(mdjvu_bitmap_t, unsigned char *, int32);

/* Same as mdjvu_bitmap_unpack_row, but writes exactly 1 for black. */
MDJVU_FUNCTION void
    mdjvu_bitmap_unpack_row_0_or_1(mdjvu_bitmap_t, unsigned char *, int32);

/* Fill the shape's row with the given array of <width> bytes.
 * 0 is white, nonzero is black.
 */
MDJVU_FUNCTION void
    mdjvu_bitmap_pack_row(mdjvu_bitmap_t, unsigned char *, int32 y);

/* Copy given bytes from or to the given shape.
 * The given array should contain height rows, top to bottom, by width bytes.
 */
MDJVU_FUNCTION void mdjvu_bitmap_pack_all(mdjvu_bitmap_t, unsigned char **);
MDJVU_FUNCTION void mdjvu_bitmap_unpack_all(mdjvu_bitmap_t, unsigned char **);
MDJVU_FUNCTION void mdjvu_bitmap_unpack_all_0_or_1
    (mdjvu_bitmap_t, unsigned char **);

/* Exchange and assign do NOT touch indices. */
MDJVU_FUNCTION void mdjvu_bitmap_assign(mdjvu_bitmap_t d, mdjvu_bitmap_t src);
MDJVU_FUNCTION void mdjvu_bitmap_exchange(mdjvu_bitmap_t d, mdjvu_bitmap_t src);
MDJVU_FUNCTION void mdjvu_bitmap_clear(mdjvu_bitmap_t);
MDJVU_FUNCTION mdjvu_bitmap_t mdjvu_bitmap_crop
    (mdjvu_bitmap_t b, int32 left, int32 top, int32 w, int32 h);
MDJVU_FUNCTION mdjvu_bitmap_t mdjvu_bitmap_clone(mdjvu_bitmap_t b);
MDJVU_FUNCTION int mdjvu_bitmap_match(mdjvu_bitmap_t img1, mdjvu_bitmap_t img2);

MDJVU_FUNCTION void mdjvu_bitmap_get_bounding_box
    (mdjvu_bitmap_t b, int32 *pl, int32 *pt, int32 *pw, int32 *ph);

/* Remove white edges from the bitmap.
 * The given bitmap is changed.
 * The left top corner of the bitmap's bounding box is written into *x and *y.
 * Blits that access this shape may become invalid (add x and y to them).
 */
MDJVU_FUNCTION void mdjvu_bitmap_remove_margins
    (mdjvu_bitmap_t, int32 *x, int32 *y);

/* Count the number of black pixels in the bitmap.
 * The results are not cached, so it uses O(width * height) time each call.
 */
MDJVU_FUNCTION int32 mdjvu_bitmap_get_mass(mdjvu_bitmap_t);

#ifdef MINIDJVU_WRAPPERS
    struct MinidjvuBitmap
    {
        inline static MinidjvuBitmap *create(int32 w, int32 h)
            { return mdjvu_bitmap_create(w,h); }
        inline void destroy()
            { mdjvu_bitmap_destroy(this); }

        inline int32 get_width()
            { return mdjvu_bitmap_get_width(this); }
        inline int32 get_height()
            { return mdjvu_bitmap_get_height(this); }

        inline int32 get_index()
            { return mdjvu_bitmap_get_index(this); }
        inline void set_index(int32 new_value)
            { mdjvu_bitmap_set_index(this, new_value); }

        inline int32 get_packed_row_size()
            { return mdjvu_bitmap_get_packed_row_size(this); }

        inline unsigned char *access_packed_row(int32 y)
            { return mdjvu_bitmap_access_packed_row(this, y); }

        inline void pack_row(unsigned char *buf, int32 y)
            { mdjvu_bitmap_pack_row(this, buf, y); }

        inline void unpack_row(unsigned char *buf, int32 y)
            { mdjvu_bitmap_unpack_row(this, buf, y); }

        inline void unpack_row_0_or_1(unsigned char *buf, int32 y)
            { mdjvu_bitmap_unpack_row_0_or_1(this, buf, y); }

        inline void pack_all(unsigned char **buf2d)
            { mdjvu_bitmap_pack_all(this, buf2d); }

        inline void unpack_all(unsigned char **buf2d)
            { mdjvu_bitmap_unpack_all(this, buf2d);}

        inline void unpack_all_0_or_1(unsigned char **buf2d)
            { mdjvu_bitmap_unpack_all_0_or_1(this, buf2d);}

        inline void assign(mdjvu_bitmap_t src)
            { mdjvu_bitmap_assign(this, src); }
        inline void clear()
            { mdjvu_bitmap_clear(this); }

        inline mdjvu_bitmap_t crop(int32 l, int32 t, int32 w, int32 h)
            { return mdjvu_bitmap_crop(this, l, t, w, h); }

        inline void remove_margins(int32 *x, int32 *y)
            { mdjvu_bitmap_remove_margins(this, x, y); }

        inline int32 get_mass()
            { return mdjvu_bitmap_get_mass(this); }
    };
#endif
