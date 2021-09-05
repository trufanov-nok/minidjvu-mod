/*
 * frames.c - extracting frameworks and calculating "importance rating"
 */

/* Frameworks are funny things...
 * The algorithm is to be commented yet.
 * Here's a picture illustrating what is a frame
 * (view with monospace font):
 *
 *  Original letter:        Its framework:
 *
  .....@@@@@@@@........ .....................
  ...@@@@@@@@@@@@...... ......@@@@@..........
  ..@@@@@@@@@@@@@@..... .....@@...@@@........
  ..@@@@@...@@@@@@@.... ....@@......@@.......
  ..@@@@.....@@@@@@.... ....@........@.......
  .@@@@@.....@@@@@@.... ....@........@.......
  .@@@@@.....@@@@@@.... ....@........@.......
  ..@@@@.....@@@@@@.... ....@........@.......
  ..........@@@@@@@.... .............@.......
  .......@@@@@@@@@@.... .............@.......
  .....@@@@@@@@@@@@.... ........@@@@@@.......
  ...@@@@@@@@@@@@@@.... ......@@@....@.......
  ..@@@@@@@..@@@@@@.... .....@@......@.......
  .@@@@@@....@@@@@@.... ...@@@.......@.......
  .@@@@@.....@@@@@@.... ...@.........@.......
  @@@@@......@@@@@@.... ..@@.........@.......
  @@@@@......@@@@@@.... ..@..........@.......
  @@@@@.....@@@@@@@.... ..@..........@.......
  @@@@@@....@@@@@@@.@@@ ..@@.........@.......
  .@@@@@@@@@@@@@@@@@@@@ ...@@.....@@@@@......
  .@@@@@@@@@@@@@@@@@@@@ ....@@..@@@...@@@@@..
  ..@@@@@@@@@.@@@@@@@@. .....@@@@............
  ....@@@@@....@@@@@... .....................

 * A letter is converted into grayshades,
 *   and a frame is the set of its purely black pixels after the transformation.
 * In the grayshade version of a letter,
 *   all pixels that were white remain absolutely white,
 *   the frame is black and the blackness falls down from it to the border.
 */

/* Offtopic: I wonder if this thing could help OCRing because frameworks
 * perfectly retain readability while becoming essentially 1-dimensional.
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

/* Stuff for not using malloc in C++
 * (made by Leon Bottou; has no use in minidjvu-mod,
 * but left here for potential DjVuLibre compatibility)
 */
#ifdef __cplusplus
# define MALLOC(Type)    new Type
# define FREE(p)         delete p
# define MALLOCV(Type,n) new Type[n]
# define FREEV(p)        delete [] p
#else
# define MALLOC(Type)    ((Type*)malloc(sizeof(Type)))
# define FREE(p)         do{if(p)free(p);}while(0)
# define MALLOCV(Type,n) ((Type*)malloc(sizeof(Type)*(n)))
# define FREEV(p)        do{if(p)free(p);}while(0)
#endif


/* This determines the gray level of the border (ratio of black).
 * Setting it to 1 will effectively eliminate grayshading.
 */
#define BORDER_FALLOFF 0.95 /* this is the main constant in all the matcher... */

typedef unsigned char byte;

static int donut_connectivity_test(byte ul, byte u, byte ur, byte l, byte r, byte dl, byte d, byte dr)/*{{{*/
{
    /*(on the pictures below 0 is white, 1 is black or gray)
     *
     * 01.
     * 1 . -> 1
     * ...
     *
     * .0.
     * 1 1 ->  1
     * .0.
     *
     * all others -> 0
     */

    int sum = u + d + l +r;

    switch(sum)
    {
        case 3:/*{{{*/
        {
            int x = 6 - (u + (l << 1) + d + (d << 1));
            switch(x)
            {
                case 0: /* l */
                    return ul && dl ? 0 : 1;
                case 1: /* d */
                    return dl && dr ? 0 : 1;
                case 2: /* r */
                    return ur && dr ? 0 : 1;
                case 3: /* u */
                    return ul && ur ? 0 : 1;
                default: assert(0); return 0;
            }
        }
        break;/*}}}*/
        case 2:/*{{{*/
        {
            int s = l + r;
            if (s & 1)
            {
                /*   A1.
                 *   1 0 - should be !A (2x2 square extermination)
                 *   .0.
                 */
                if (l)
                {
                    if (u)
                        return ul ? 0 : 1;
                    else
                        return dl ? 0 : 1;
                }
                else /* r */
                {
                    if (u)
                        return ur ? 0 : 1;
                    else
                        return dr ? 0 : 1;
                }
            }
            else
            {
                /*   .0.
                 *   1 1 - surely should be 1 to preserve connection
                 *   .0.
                 */
                return 1;
            }
        }
        break;/*}}}*/
        case 0: case 4:
            return 1;
        case 1:
            return 0;
        default: assert(0); return 0;
    }
}/*}}}*/

static byte donut_transform_pixel(byte *upper, byte *row, byte *lower)/*{{{*/
{
    /* (center pixel should be gray in order for this to work)
     * (on the pictures below 0 is white, 1 is black or gray)
     *
     * 01.
     * 1 . -> center will become 1
     * ...
     *
     * .0.
     * 1 1 -> center will become 1
     * .0.
     *
     * 00.
     * 1 0 -> center will become 1
     * .0.
     *
     * 1..
     * 1 0 -> center will become 0
     * 1..
     *
     * 11.
     * 1 0 -> center will become 0
     * .0.
     *
     * .A.
     * A A -> center will become 1
     * .A.
     */

    int sum, l, u, d, r;
    if (!*row) return 0;

    sum = (u = *upper ? 1 : 0) + (d = *lower ? 1 : 0) +
          (l = row[-1] ? 1 : 0) + (r = row[1] ? 1 : 0);

    switch(sum)
    {
        case 1: case 3:/*{{{*/
        {
            int x = u + (l << 1) + d + (d << 1);
            if (sum == 3) x = (6 - x) ^ 2;
            switch(x)
            {
                case 0: /* r */
                    return upper[1] && lower[1] ? 0 : 1;
                case 1: /* u */
                    return upper[-1] && upper[1] ? 0 : 1;
                case 2: /* l */
                    return upper[-1] && lower[-1] ? 0 : 1;
                case 3: /* d */
                    return lower[-1] && lower[1] ? 0 : 1;
                default: assert(0); return 0;
            }
        }
        break;/*}}}*/
        case 2:/*{{{*/
        {
            int s = l + r;
            if (s & 1)
            {
                /*   A1.
                 *   1 0 - should be !A (2x2 square extermination)
                 *   .0.
                 */
                if (l)
                {
                    if (u)
                        return upper[-1] ? 0 : 1;
                    else
                        return lower[-1] ? 0 : 1;
                }
                else /* r */
                {
                    if (u)
                        return upper[1] ? 0 : 1;
                    else
                        return lower[1] ? 0 : 1;
                }
            }
            else
            {
                /*   .0.
                 *   1 1 - surely should be 1 to preserve connection
                 *   .0.
                 */
                return 1;
            }
        }
        break;/*}}}*/
        case 0: case 4:
            return 1; /* lone pixels are NOT omitted */
        default: assert(0); return 0;
    }
}/*}}}*/

/* `pixels' should have a margin of 1 pixel at each side
 * returns true if the image was changed
 */
static int flay(byte **pixels, byte *buf, int w, int h, int rank, int **ranks)
{
    assert(pixels);
    assert(buf);
    assert(ranks);
    assert(w);
    assert(h);

    int result = 0;

    for (int i = 0; i < h; i++) {
        byte *up = pixels[i-1], *row = pixels[i], *dn = pixels[i+1];
        byte *buf_up = (i>0)? buf + w * (i-1): NULL;
        byte *buf_row = buf + w * i;

        for (int j = 0; j < w; j++)
        {
            buf_row[j] = donut_transform_pixel(up + j, row + j, dn + j);

            if (row[j] && !buf_row[j])
            {
                if (!donut_connectivity_test( (buf_up && j>0)?(buf_up[j-1]?1:0):0, (buf_up)?(buf_up[j]?1:0):0, (buf_up && j<w-1)?(buf_up[j+1]?1:0):0,
                                              (j>0)?(buf_row[j-1]?1:0):0, row[j+1]?1:0,
                                              dn[j-1]?1:0, dn[j]?1:0, dn[j+1]?1:0))
                {
                    ranks[i][j] = rank;
                    result = 1;
                } else {
                    buf_row[j] = 1;
                }
            }
        }
    }

    for (int i = 0; i < h; i++)
    {
        memcpy(pixels[i], buf+w*i, w);
    }

    return result;
}

/* TODO: use less temporary buffers and silly copyings */
MDJVU_IMPLEMENT void mdjvu_soften_pattern(byte **result, byte **pixels, int32 w, int32 h)/*{{{*/
{
    byte *r = MALLOCV(byte, (w + 2) * (h + 2));
    byte **pointers = MALLOCV(byte *, h + 2);
    int *ranks_buf = MALLOCV(int, w * h);
    int **ranks = MALLOCV(int *, h);

    int i, j, passes = 1;
    double level = 1, falloff;
    byte *colors;

    memset(r, 0, (w + 2) * (h + 2));
    memset(ranks_buf, 0, w * h * sizeof(int));

    for (i = 0; i < h + 2; i++)
        pointers[i] = r + (w + 2) * i + 1;

    for (i = 0; i < h; i++)
        memcpy(pointers[i+1], pixels[i], w);

    for (i = 0; i < h; i++)
        ranks[i] = ranks_buf + w * i;

    byte *buf = MALLOCV(byte, w * h);
    while(flay(pointers + 1, buf, w, h, passes, ranks)) passes++;
    FREEV(buf);

    colors = MALLOCV(byte, passes + 1);

    falloff = pow(BORDER_FALLOFF, 1./passes);

    for (i = 0; i < passes; i++)
    {
        colors[i] = (byte) (level * 255);
        level *= falloff;
    }
    /* TODO: test the next line */
    /* colors[passes - 1] = 50; pay less attention to border pixels */
    colors[passes] = 0;

    pointers++;
    for (i = 0; i < h; i++)
    {
        for (j = 0; j < w; j++)
        {
            if (pointers[i][j])
            {
                result[i][j] = 255;
            }
            else
            {
                result[i][j] = colors[passes - ranks[i][j]];
            }
        }
    }
    pointers--;

    FREEV(colors);
    FREEV(ranks);
    FREEV(ranks_buf);
    FREEV(r);
    FREEV(pointers);
}/*}}}*/
