/*
 * smooth.c - pre-filtering bitmap before splitting
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdlib.h>
#include <string.h>
#ifndef _MSC_VER
#include <endian.h>
#endif

/* all input rows must be 0-or-1 unpacked */
//static void smooth_row(unsigned char *r, /* result    */
//                       unsigned char *u, /* upper row */
//                       unsigned char *t, /* this row - must have margin 1 */
//                       unsigned char *l, /* lower row */
//                       int32 n)
//{
//    int32 i;
//    for (i = 0; i < n; i++)
//    {
//	int score = u[i] + l[i] + t[i-1] + t[i+1];

//	if (t[i] == 0)
//	{
//	    /* Only turn white into black for a good reason. */
//	    r[i] = (score == 4);
//	}
//	else if (score == 0)
//	{
//	    /* Good reason to clear the pixel. */
//	    r[i] = 0;
//	}
//	else if (score == 1)
//	{
//	    /* Check for weak horizontal or vertical linking. */
//	    if (u[i] | l[i])
//		r[i] = (u[i-1] & l[i-1]) | (u[i+1] & l[i+1]);
//	    else
//		r[i] = (u[i-1] & u[i+1]) | (l[i-1] & l[i+1]);
//	}
//	else
//	{
//	    /* Keep it black. */
//	    r[i] = 1;
//	}

//    }
//}

__inline size_t get_smooth(size_t u, size_t t, size_t d)
{
    size_t ul = u >> 1, l = t >> 1, dl = d >> 1;
    size_t ur = u << 1, r = t << 1, dr = d << 1;

    size_t res0 = l & r & u & d; // score 4 regardles t
    size_t res1 = /*t &*/ ( (l & u)  |  (u & r) |  (r & d) | (l & d) | (l & r) | (u & d) ); // score >= 2

    // score 1
    size_t res2 = /*t &*/ (u | d) & ( (ul & dl) | (ur & dr) );
    res2 |= /*t &*/ (l | r) & ( (ul & ur) | (dl & dr) );

    return res0 | (t & (res1 | res2));
}

#if __BYTE_ORDER == __BIG_ENDIAN
__inline size_t swap_t(size_t* val, unsigned int size) {
    if (size >= sizeof (size_t))
        return *val;

    size_t res = 0;
    memcpy(&res, val, size);
    return res;
}
#elif __BYTE_ORDER == __LITTLE_ENDIAN
__inline size_t swap_t(size_t* val, unsigned int size) {
    size_t res = 0;
    unsigned char * a = (unsigned char *) &res;
    memcpy(&res, val, size);
    unsigned char t;
	for (unsigned int i = 0; i < (sizeof (size_t)/2); i++) {
        t = a[i];
        a[i] = a[sizeof (size_t) - i - 1];
        a[sizeof (size_t) - i - 1] = t;
    }
    return res;
}
#endif

static void smooth_row(unsigned char *r, /* result    */
                       unsigned char *u, /* upper row */
                       unsigned char *t, /* this row - must have margin 1 */
                       unsigned char *l, /* lower row */
                       int32 n)
{
    if ( !n ) return;
    const size_t int_len_in_bits = sizeof (size_t)*8;
    const size_t len = (n + (int_len_in_bits -1) ) / int_len_in_bits;
    const size_t tail_len = (n % int_len_in_bits) ? ((n % int_len_in_bits) + 7) >> 3 : sizeof (size_t);

    size_t *r_p = (size_t *) r; /* result    */
    size_t *u_p = (size_t *) u;
    size_t *t_p = (size_t *) t;
    size_t *l_p = (size_t *) l;

    size_t u_buf = 0, t_buf = 0, l_buf = 0;
    size_t u_val = 0, t_val = 0, l_val = 0;
    size_t u_cur = 0, t_cur = 0, l_cur = 0;

    const size_t mask1 = (~(size_t)0x0) << 1; //0b11111..110
    const size_t mask2 = (size_t)0x01 << (int_len_in_bits-1); //0b100000.00
    const size_t mask3 = mask2 >> 1; //0b01000..00
    const size_t mask4 = mask3 >> 1; //0b00100..00
    const size_t mask5 = mask1 & ~mask2; //0b01111..10


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

        size_t res = get_smooth(u_val, t_val, l_val);

        size_t tail = res & mask3;
        size_t head = res & mask4;

        if (tail) {
            // for i == 0 tail is always false and this is not called
            // we access last byte instead of size_t* to not mess with int endiannes
			*(((unsigned char *)r_p)-1) |= 1;; // last bit is always 0 bcs of mask5
        }


        res = get_smooth(u_cur, t_cur, l_cur);
        res &= mask5;
        if (head) {
            res |= mask2;
        }

        res = swap_t(&res, /*i==len-1?tail_len:*/sizeof (size_t));
        memcpy(r_p++, &res, (i==len-1)?tail_len:sizeof (size_t));
    }


    if (n % 8) {
        r += ((n+7)>>3)-1;
        *r &= 0xFF << (8- n % 8);
    }
}

MDJVU_IMPLEMENT void mdjvu_smooth(mdjvu_bitmap_t b)
{
    int32 w = mdjvu_bitmap_get_width(b);
    int32 h = mdjvu_bitmap_get_height(b);
    int32 i;
    unsigned char *u = NULL, /* upper row */
                  *t = NULL, /* this row */
                  *l = NULL; /* lower row */

    if (h < 3) return;

    const int32 row_size = mdjvu_bitmap_get_packed_row_size(b);

    unsigned char *r = (unsigned char *) malloc(row_size*h); /* result */

    l = mdjvu_bitmap_access_packed_row(b, 0);
    for (i = 0; i < h; i++) {
        u = t;
        t = l;

        if (i + 1 < h)
            l = mdjvu_bitmap_access_packed_row(b, i+1);
        else
            l = NULL;

        smooth_row(r+row_size*i, u, t, l, w);
    }


    memcpy(mdjvu_bitmap_access_packed_row(b, 0), r, row_size*h);

    free(r);
}
