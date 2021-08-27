/*
 * split.c - splitting bitmaps to letters
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

/* _______________________   managing splitter options   ___________________ */

mdjvu_split_options_t mdjvu_split_options_create(void)
{
    int32 *p = (int32 *) malloc(sizeof(int32));
    mdjvu_init();
    *p = 0;
    return (mdjvu_split_options_t) p;
}

void mdjvu_split_options_set_maximum_shape_size(mdjvu_split_options_t opt, int32 size)
{
    assert(size > 0);
    * (int32 *) opt = size;
}

void mdjvu_split_options_destroy(mdjvu_split_options_t opt)
{
    free(opt);
}

/* _________________________   interpreting runs   _________________________ */

/* Runs are stored in a "map".
 * A map entry may be 0 (same-as-prev), BLACK_RUN(1) and WHITE_RUN(2).
 * Each run is started by BLACK_RUN or WHITE_RUN and continues with zeros.
 * First run is white.
 *
 * Here's an example:
 *
 *               map: 001002001000210
 *  interpreted line: 001110001111011
 */
#define BLACK_RUN 1
#define WHITE_RUN 2

/* This function does three things simultaneously:
 *    1) renders `map' into `result_line' and applies AND with line
 *    2) deletes black pixels in `line'
 *    3) clears `map' including the last (max_x + 1) column
 *
 * Indices in `map'         are min_x..max_x.
 * Indices in `line'        are min_x..max_x.
 * Indices in `result_line' are 0..max_x-min_x.
 *
 */
static void interpret_runs_in_a_line(int32 min_x, int32 max_x,
                                     unsigned char *map,
                                     unsigned char *result_line,
                                     unsigned char *line)
{
    int32 i = min_x;
    memset(result_line, 0, max_x - min_x + 1);
    map[max_x + 1] = 0;
    while (1)
    {
        /* white run */
        while (map[i] != BLACK_RUN)
        {
            if (++i > max_x) return;
        }

        map[i] = 0;

        /* black run */
        while (map[i] != WHITE_RUN)
        {
            result_line[i - min_x] = line[i]; /* copy instead of paint */
            line[i] = 0;
            if (++i > max_x) return;
        }

        map[i] = 0;
    }
}

/* This function does the same as interpret_runs_in_a_line(),
 *    but on a rectangle.
 * Returns the rendered result.
 * Also clears it in pixels[][].
 * `map' is cleared.
 */
static mdjvu_bitmap_t interpret_runs(int32 min_x, int32 max_x,
                                     int32 min_y, int32 max_y,
                                     unsigned char **map,
                                     unsigned char **pixels)
{
    int32 w = max_x - min_x + 1;
    int32 h = max_y - min_y + 1;
    unsigned char *line_buf = (unsigned char *) malloc(w);
    mdjvu_bitmap_t bmp = mdjvu_bitmap_create(w, h);
    int32 y;
    for (y = min_y; y <= max_y; y++)
    {
        interpret_runs_in_a_line(min_x, max_x, map[y], line_buf, pixels[y]);
        mdjvu_bitmap_pack_row(bmp, line_buf, y - min_y);
    }
    free(line_buf);
    return bmp;
}

/* ________________________   walking along contours   _____________________ */

/* (start_x, start_y) should be the leftmost of the topmost black pixels.
 * Walking is counterclockwise.
 * map must be cleared.
 */

static void walk_around_a_black_contour(unsigned char **pixels,
                                        unsigned char **map,
                                        int32 start_x, int32 start_y,
                                        int32 *pmin_x, int32 *pmax_x,
                                        int32 *pmin_y, int32 *pmax_y)
{
    int32 x = start_x, y = start_y;   /* current black pixel */
    int32 dx = 0, dy = 1;  /* vector along the contour*/
                   /* it's always true that abs(dx) == 1 && dy == 0
                    *                        ||
                    *                       abs(dy) == 1 && dx == 0
                    */

    /* Bounding box. start_y will always be minimal y. */
    int32 min_x = start_x;
    int32 max_x = start_x;
    int32 max_y = start_y;

    /* Starting picture:
     *
     *      | w | w |
     *      ----+---|
     *      | w | b |  <- start_y
     *            ^
     *           start_x
     *
     *  and we're heading down
     */

    /* Generic picture (===> denotes contour fragment along (dx,dy)):
     *
     *  (1) | b | ? | (3)
     *      ====>----
     *  (2) | w | ? | (4)
     *
     *  ===> vector is (dx,dy)
     *  up vector (that is, (dx,dy) rotated counterclockwise) is (dy, -dx).
     *  down vector (that is, (dx,dy) rotated clockwise) is      (-dy, dx).
     *  (the coordinate system is left)
     */

    /* MAIN IMAGE WALKING LOOP */
    while (x != start_x || y != start_y || dx != -1)
    {
        assert(pixels[y][x] /* is black */);

        /* see if this pixel enlarges bbox */
        if (x < min_x) min_x = x; else if (x > max_x) max_x = x;
        if (y > max_y) max_y = y;

        /* fill contour map if needed */
        if (dx == 0)
        {
            if (dy > 0)
                map[y][x] = BLACK_RUN;
            else
                map[y][x + 1] = WHITE_RUN; /* map has the right margin */
        }

        /* proceed with contour */
        if (pixels[y + dy][x + dx])  /* testing (3) */
        {
            int32 x4 = x + dx - dy, y4 = y + dy + dx;

            if (pixels[y4][x4]) /* testing (4) */
            {
                /* (3) and (4) are black - right turn */
                /* (-dy,dx) becomes new (dx,dy) (coordinate system is left) */
                /* pixel (4) becomes (x,y) */
                int32 old_dx = dx;
                dx = -dy;
                dy = old_dx;
                x = x4; y = y4;
            }
            else
            {
                /* go straight: pixel (3) becomes (x,y) */
                x += dx; y += dy;
            }
        }
        else
        {
            /* (3) is white - left turn */
            /* (dy,-dx) becomes new (dx,dy) (coordinate system is left) */
            int32 old_dx = dx;
            dx = dy;
            dy = -old_dx;
        }
    }
    *pmin_x = min_x;
    *pmax_x = max_x;
    *pmin_y = start_y;
    *pmax_y = max_y;
}

/* ____________________________   process row   ____________________________ */

static void add_to_image(mdjvu_image_t image,
                         mdjvu_bitmap_t bitmap,
                         int32 dpi,
                         mdjvu_split_options_t opt,
                         int32 blit_shift_x,
                         int32 blit_shift_y,
                         int big);


/* This function scans the y-th row for contours. */
/* Margins (1 pixel) from each side are required. */
static void process_row(unsigned char **pixels, unsigned char **map,
                        int32 y, int32 w, int32 h,
                        mdjvu_image_t image, int32 shift_x, int32 shift_y,
                        int32 max_shape_width,
                        int32 blit_shift_x, int32 blit_shift_y,
                        int32 dpi, mdjvu_split_options_t opt,
                        int big)
{
    unsigned char *row = pixels[y];
    int32 i;
    for (i = 0; i < w; i++)
    {
        if (row[i])
        {
            int32 min_x, max_x, min_y, max_y, shape_width;

            /* extract the contour */
            mdjvu_bitmap_t bitmap;
            walk_around_a_black_contour(pixels, map, i, y,
                                        &min_x, &max_x, &min_y, &max_y);
            bitmap = interpret_runs(min_x, max_x,
                                    min_y, max_y,
                                    map, pixels);
            shape_width = mdjvu_bitmap_get_width(bitmap);
            assert(shape_width == max_x - min_x + 1);
            assert(mdjvu_bitmap_get_height(bitmap) == max_y - min_y + 1);
            if (shape_width <= max_shape_width)
            {
                mdjvu_image_add_bitmap(image, bitmap);
                mdjvu_image_add_blit(image, shift_x + min_x + blit_shift_x,
                                            y + shift_y + blit_shift_y, bitmap);
                mdjvu_image_set_suspiciously_big_flag(image, bitmap, big);
            }
            else
            {
                /* further split the bitmap */
                int32 number_of_chunks = (shape_width + max_shape_width - 1)
                                            /
                                          max_shape_width;
                int32 j;
                int32 shape_height = mdjvu_bitmap_get_height(bitmap);
                for (j = 0; j < number_of_chunks; j++)
                {
                    int32 chunk_x = shape_width * j / number_of_chunks;
                    mdjvu_bitmap_t chunk = mdjvu_bitmap_crop(bitmap,
                      chunk_x, 0,
                      shape_width * (j+1) / number_of_chunks - chunk_x,
                      shape_height
                    );
                    /* After splitting, some white margins may be left,
                     * or the bitmap may lose connectivity.
                     * Apply the algorithm recursively to the chunk.
                     */
                    add_to_image(image, chunk, dpi, opt,
                                 shift_x + chunk_x + min_x + blit_shift_x,
                                 y + shift_y + blit_shift_y, /* big: */ 1);
                    mdjvu_bitmap_destroy(chunk);
                }
                mdjvu_bitmap_destroy(bitmap);
            } /* if (shape_width <= max_shape_width) */
        }
    }
}

/* _________________________   the main routines   _________________________ */

/* Not much job left to do here, mostly moving a window through the image. */

static void add_to_image(mdjvu_image_t image,
                         mdjvu_bitmap_t bitmap,
                         int32 dpi,
                         mdjvu_split_options_t opt,
                         int32 blit_shift_x,
                         int32 blit_shift_y,
                         int big)
{
    int32 max_shape_size = opt ? * (int32 *) opt : 0;
    int32 width = mdjvu_bitmap_get_width(bitmap);
    int32 height = mdjvu_bitmap_get_height(bitmap);
    unsigned char **buf, **window_base, **window_buf, **map;
    int32 window_shift = 0, y = 0, i;

    if (!max_shape_size)
        max_shape_size = dpi;
    int32 max_shape_height = max_shape_size;
    if (max_shape_height > height)
        max_shape_height = height;

    /* n-th line will be unpacked into buf[n % max_shape_height] + 1.
     * ( +1 is to make the left margin)
     * buf[max_shape_height] will always be blank.
     *
     * window_base[window_shift - 1]
     *      points to buf[max_shape_height] + 1 (blank line) - top margin
     * window_base[window_shift + max_shape_height]
     *      points to buf[max_shape_height] + 1 (blank line) - bottom margin
     * window_base[window_shift + i]
     *      points to buf[(window_shift + i) % max_shape_height] + 1.
     */

    /* map has the right margin of 1 */
    map = mdjvu_create_2d_array(width + 1, max_shape_height);

    /* buf has left, right and bottom margins of 1 */
    buf = mdjvu_create_2d_array(width + 2, max_shape_height + 1);
    window_buf = (unsigned char **)
        malloc(2 * (max_shape_height + 2) * sizeof(unsigned char *));
    window_base = window_buf + 1;

    /* Unpack initial portion of the bitmap; bind the window to the buffer */
    for (i = 0; i < max_shape_height; i++)
    {
        window_base[i] = window_base[max_shape_height + i] = buf[i] + 1;
        mdjvu_bitmap_unpack_row(bitmap, buf[i] + 1, i);
    }

    /* Setup top and bottom white margins */
    window_base[-1] =
        window_base[2 * max_shape_height - 1] =
            buf[max_shape_height] + 1;

    /* The "window moving" loop.
     * We're moving a (width x max_shape_height) window through the image.
     */
    while(1)
    {
        /* Extract some shapes from the window
         * (shapes touching the topmost row will be affected).
         */

        unsigned char *top_margin_save =    /* make the top margin */
            window_base[window_shift - 1];  /* (save what was there) */
        unsigned char *bot_margin_save =    /* same with the bottom margin */
            window_base[window_shift + max_shape_height];
        int32 old_window_shift;
        window_base[window_shift - 1] = buf[max_shape_height] + 1; /* clear them */
        window_base[window_shift + max_shape_height] = buf[max_shape_height] + 1;

        process_row(window_base + window_shift, map,
                    /* index of a row to process: */ 0,
                    width, max_shape_height, image, 0, y, max_shape_size,
                    blit_shift_x, blit_shift_y, dpi, opt, big);

        window_base[window_shift - 1] = top_margin_save; /* restore margins */
        window_base[window_shift + max_shape_height] = bot_margin_save;

        /* Shift the window */
        y++;
        old_window_shift = window_shift;
        window_shift = y % max_shape_height;
        if (y + max_shape_height > height)
            break;

        /* Unpack a new row into the bottom window row */
        mdjvu_bitmap_unpack_row(bitmap, buf[old_window_shift] + 1,
                                        y + max_shape_height - 1);
    }

    /* Process the last window fully */
    for (i = 0; y + i < height; i++)
    {
        process_row(window_base + window_shift, map,
                    /* index of a row to process: */ i,
                    width, max_shape_height, image, 0, y, max_shape_size,
                    blit_shift_x, blit_shift_y, dpi, opt, big);
    }

    /* Clean up */
    free(window_buf);
    mdjvu_destroy_2d_array(map);
    mdjvu_destroy_2d_array(buf);
}

mdjvu_image_t
mdjvu_split(mdjvu_bitmap_t bitmap, int32 dpi, mdjvu_split_options_t opt)
{
    int32 width = mdjvu_bitmap_get_width(bitmap);
    int32 height = mdjvu_bitmap_get_height(bitmap);
    mdjvu_image_t result = mdjvu_image_create(width, height);
    mdjvu_image_enable_suspiciously_big_flags(result);
    mdjvu_image_set_resolution(result, dpi);
    add_to_image(result, bitmap, dpi, opt, 0, 0, /* big: */ 0);
    return result;
}
