/*
 * split.c - splitting bitmaps to letters
 */

#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

// This is a Connected Components detection algorithm from
// DjVuLibre's cjb2 encoder.
//
// https://github.com/barak/djvulibre/blob/master/tools/cjb2.cpp


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


// --------------------------------------------------
// UTILITIES
// --------------------------------------------------

#ifdef MIN
#undef MIN
#endif

inline int
MIN(int a, int b)
{
    return ( a<b ?a :b);
}

#ifdef MAX
#undef MAX
#endif

inline int
MAX(int a, int b)
{
    return ( a>b ?a :b);
}

// --------------------------------------------------
// CONNECTED COMPONENT ANALYSIS AND CLEANING
// --------------------------------------------------

struct GRect
{
    int xmin;
    int xmax;
    int ymin;
    int ymax;
};


// -- A run of black pixels
struct Run
{
    int y;         // vertical coordinate
    short x1;      // first horizontal coordinate
    short x2;      // last horizontal coordinate
    int ccid;      // component id
};

// -- Compares runs
static inline int
comp_runs (const struct Run *a, const struct Run *b)
{
    return (a->y<b->y) || (a->y==b->y && a->x1<=b->x1);
}

// -- A component descriptor
struct CC
{
    struct GRect bb;      // bounding box
    int npix;      // number of black pixels
    int nrun;      // number of runs
    int frun;      // first run in cc ordered array of runs
};

// -- An image composed of runs
struct CCImage
{
    int height;            // Height of the image in pixels
    int width;             // Width of the image in pixels
    int dpi;

    struct Run* runs;             // array of runs
    int  runs_count, runs_allocated;

    struct CC*  ccs;              // Array of component descriptors
    int  ccs_count, ccs_allocated;

    int nregularccs;       // Number of regular ccs (set by merge_and_split_ccs)
    int largesize;         // CCs larger than that are special
    int smallsize;         // CCs smaller than that are special
    int tinysize;          // CCs smaller than that may be removed
};

// -- Adds a run to the CCImage
void
ccimage_add_single_run(struct CCImage* image, int y, int x1, int x2, int ccid)
{
    if (image->runs_count >= image->runs_allocated) {
        image->runs_allocated <<= 1;
        image->runs = (struct Run *) realloc(image->runs,
                                             image->runs_allocated * sizeof(struct Run));
    }

    struct Run* run = &image->runs[image->runs_count++];
    run->y = y;
    run->x1 = x1;
    run->x2 = x2;
    run->ccid = ccid;
}

// -- Adds runs extracted from a bitmap
void
ccimage_add_bitmap_runs(struct CCImage* image, const mdjvu_bitmap_t bm, int offx, int offy, int ccid)
{
    int w = mdjvu_bitmap_get_width(bm);
    int h = mdjvu_bitmap_get_height(bm);
    unsigned char * row = MDJVU_MALLOCV(unsigned char, w);

    // Iterate over rows
    for (unsigned int y=0; y<h; y++)
    {
        mdjvu_bitmap_unpack_row(bm, row, y);

        int x = 0;
        // Iterate over runs
        while (x < w)
        {
            while (x < w  && !row[x]) x++;
            if (x < w)
            {
                int x1 = x;
                while (x < w && row[x]) x++;
                ccimage_add_single_run(image, offy+y, offx+x1, offx+x-1, ccid);
            }
        }
    }
    MDJVU_FREEV(row);
}



void sort_array(struct Run* data, int size, int lo, int hi)
{
    while(1)
    {
        if (hi <= lo)
            return;
        if (hi >= size || lo<=0)
            return;
        // Test for insertion sort
        if (hi <= lo + 50)
        {
            for (int i=lo+1; i<=hi; i++)
            {
                int j = i;
                struct Run tmp = data[i];
                while ((--j>=lo) && !comp_runs (&data[j],&tmp))
                    data[j+1] = data[j];
                data[j+1] = tmp;
            }
            return;
        }
        // -- determine median-of-three pivot
        struct Run tmp = data[lo];
        struct Run pivot = data[(lo+hi)/2];
        if (comp_runs (&pivot,&tmp))
        { tmp = pivot; pivot=data[lo]; }
        if (comp_runs (&data[hi],&tmp))
        { pivot = tmp; }
        else if (comp_runs (&data[hi], &pivot))
        { pivot = data[hi]; }
        // -- partition set
        int h = hi;
        int l = lo;
        while (l < h)
        {
            while (! comp_runs (&pivot, &data[l])) l++;
            while (! comp_runs (&data[h], &pivot)) h--;
            if (l < h)
            {
                tmp = data[l];
                data[l] = data[h];
                data[h] = tmp;
                l = l+1;
                h = h-1;
            }
        }
        // -- recurse, small partition first
        //    tail-recursion elimination
        if (h - lo <= hi - l) {
            sort_array(data,size,lo,h);
            lo = l; // sort(l,hi)
        } else {
            sort_array(data,size,l,hi);
            hi = h; // sort(lo,h)
        }
    }
}


// -- Performs connected component analysis
void
ccimage_make_ccids_by_analysis(struct CCImage* image)
{
    // Sort runs
    sort_array(image->runs, image->runs_count, 0, image->runs_count-1);
    // Single Pass Connected Component Analysis (with unodes)
    int n;
    int p=0;

    int umap_count = 0;
    int umap_allocated = 16;
    int* umap = MDJVU_MALLOCV(int, 16);


    for (n=0; n<image->runs_count; n++)
    {
        int y = image->runs[n].y;
        int x1 = image->runs[n].x1 - 1;
        int x2 = image->runs[n].x2 + 1;
        int id = umap_count;
        // iterate over previous line runs
        for(;image->runs[p].y < y-1;p++);
        for(;(image->runs[p].y < y) && (image->runs[p].x1 <= x2);p++ )
        {
            if ( image->runs[p].x2 >= x1 )
            {
                // previous run touches current run
                int oid = image->runs[p].ccid;
                while (umap[oid] < oid)
                    oid = umap[oid];
                if (id+1 > umap_count) {
                    id = oid;
                } else if (id < oid) {
                    umap[oid] = id;
                } else {
                    umap[id] = oid;
                    id = oid;
                }
                // freshen previous run id
                image->runs[p].ccid = id;
                // stop if previous run goes past current run
                if (image->runs[p].x2 >= x2)
                    break;
            }
        }
        // create new entry in umap
        image->runs[n].ccid = id;
        if (id+1 >= umap_count) // if (id > umap.hbound())
        {
            if (id+1 >= umap_allocated) { // umap.touch(id);
                while (id+1 >= umap_allocated) umap_allocated <<= 1;
                umap = (int*) realloc(umap, umap_allocated * sizeof(int));
            }
            umap[id] = id;
            umap_count++;
        }

    }
    // Update umap and ccid
    for (n=0; n<image->runs_count; n++)
    {
        struct Run *run = &image->runs[n];
        int ccid = run->ccid;
        while (umap[ccid] < ccid)
        {
            ccid = umap[ccid];
        }
        umap[run->ccid] = ccid;
        run->ccid = ccid;
    }
    MDJVU_FREEV(umap);
}

// -- Constructs the ``ccs'' array from run's ccids.
void
ccimage_make_ccs_from_ccids(struct CCImage* image)
{
    if (!image->runs_count) return; // empty page
    int n;
    struct Run *pruns = image->runs;
    // Find maximal ccid
    int maxccid = image->nregularccs-1;
    for (n=0; n<image->runs_count; n++)
        if (pruns[n].ccid > maxccid)
            maxccid = image->runs[n].ccid;

    // Renumber ccs
    int* armap = MDJVU_MALLOCV(int, maxccid+1);
    int *rmap = armap;
    for (n=0; n<=maxccid; n++)
        armap[n] = -1;
    for (n=0; n<image->runs_count; n++)
        if (pruns[n].ccid >= 0)
            rmap[ pruns[n].ccid ] = 1;
    int nid = 0;
    for (n=0; n<=maxccid; n++)
        if (rmap[n] > 0)
            rmap[n] = nid++;

    // Adjust nregularccs (since ccs are renumbered)
    while (image->nregularccs>0 && rmap[image->nregularccs-1]<0)
        image->nregularccs -= 1;
    if (image->nregularccs>0)
        image->nregularccs = 1 + rmap[image->nregularccs-1];

    // Prepare cc descriptors
    if (image->ccs_allocated < nid) {
        image->ccs_allocated = nid; // image->ccs->resize(0,nid-1);
        image->ccs = (struct CC*) realloc(image->ccs, image->ccs_allocated * sizeof(struct CC));
    }
    image->ccs_count = nid;

    for (n=0; n<nid; n++)
        image->ccs[n].nrun = 0;

    // Relabel runs
    for (n=0; n<image->runs_count; n++)
    {
        struct Run *run = &pruns[n];
        if (run->ccid < 0) continue;  // runs with negative ccids are destroyed
        int oldccid = run->ccid;
        int newccid = rmap[oldccid];
        struct CC *cc = &image->ccs[newccid];
        run->ccid = newccid;
        cc->nrun += 1;
    }

    // Compute positions for runs of cc
    int frun = 0;
    for (n=0; n<nid; n++)
    {
        image->ccs[n].frun = rmap[n] = frun;
        frun += image->ccs[n].nrun;
    }

    // Copy runs
    struct Run* rtmp =  MDJVU_MALLOCV(struct Run, image->runs_count); //rtmp.steal(image->runs);
    const int rtmp_hbound = image->runs_count-1;
    memcpy(rtmp, image->runs, sizeof(struct Run)*image->runs_count);
    struct Run *ptmp = rtmp;

    image->runs_allocated = frun; // image->runs.resize(0,frun-1);
    image->runs_count = frun;
    image->runs = (struct Run *) realloc(image->runs, image->runs_allocated * sizeof(struct Run));
    pruns = image->runs;

    for (n=0; n<=rtmp_hbound; n++)
    {
        int id = ptmp[n].ccid;
        if (id < 0) continue;
        int pos = rmap[id]++;
        pruns[pos] = ptmp[n];
    }

    MDJVU_FREEV(rtmp);
    MDJVU_FREEV(armap);

    // Finalize ccs
    for (n=0; n<nid; n++)
    {
        struct CC *cc = &image->ccs[n];
        int npix = 0;
        sort_array(image->runs, image->runs_count, cc->frun, cc->frun+cc->nrun-1);
        struct Run *run = &image->runs[cc->frun];
        int xmin = run->x1;
        int xmax = run->x2;
        int ymin = run->y;
        int ymax = run->y;
        for (int i=0; i<cc->nrun; i++, run++)
        {
            if (run->x1 < xmin)  xmin = run->x1;
            if (run->x2 > xmax)  xmax = run->x2;
            if (run->y  < ymin)  ymin = run->y;
            if (run->y  > ymax)  ymax = run->y;
            npix += run->x2 - run->x1 + 1;
        }
        cc->npix = npix;
        cc->bb.xmin = xmin;
        cc->bb.ymin = ymin;
        cc->bb.xmax = xmax + 1;
        cc->bb.ymax = ymax + 1;
    }
}

// Removes ccs which are too small.
void
ccimage_erase_tiny_ccs(struct CCImage* image)
{
    for (int i=0; i<image->ccs_count; i++)
    {
        struct CC* cc = &image->ccs[i];
        if (cc->npix <= image->tinysize)
        {
            // Mark cc to be erased
            struct Run *r = &image->runs[cc->frun];
            int nr = cc->nrun;
            cc->nrun = 0;
            cc->npix = 0;
            while (--nr >= 0)
                (r++)->ccid = -1;
        }
    }
}


// -- Merges small ccs and split large ccs
void
ccimage_merge_and_split_ccs(struct CCImage* image)
{
    int ncc = image->ccs_count;
    int nruns = image->runs_count;
    int splitsize = image->largesize;
    if (ncc <= 0) return;
    // Grid of special components
    int gridwidth = (image->width+splitsize-1)/splitsize;
    image->nregularccs = ncc;
    // Set the correct ccids for the runs
    for (int ccid=0; ccid<ncc; ccid++)
    {
        struct CC* cc = &image->ccs[ccid];
        if (cc->nrun <= 0) continue;
        int ccheight = cc->bb.ymax - cc->bb.ymin;
        int ccwidth = cc->bb.xmax - cc->bb.xmin;
        if (ccheight<=image->smallsize && ccwidth<=image->smallsize)
        {
            int gridi = (cc->bb.ymin+cc->bb.ymax)/splitsize/2;
            int gridj = (cc->bb.xmin+cc->bb.xmax)/splitsize/2;
            int newccid = ncc + gridi*gridwidth + gridj;
            for(int runid=cc->frun; runid<cc->frun+cc->nrun; runid++)
                image->runs[runid].ccid = newccid;
        }
        else if (ccheight>=image->largesize || ccwidth>=image->largesize)
        {
            for(int runid=cc->frun; runid<cc->frun+cc->nrun; runid++)
            {
                struct Run* r = &image->runs[runid];
                int y = r->y;
                int x_start = r->x1;
                int x_end = r->x2;
                int gridi = y/splitsize;
                int gridj_start = x_start/splitsize;
                int gridj_end = x_end/splitsize;
                int gridj_span = gridj_end-gridj_start;
                int newccid = ncc + gridi*gridwidth + gridj_start;
                if (! gridj_span)
                {
                    r->ccid = newccid;
                }
                else // gridj_span>0
                {
                    // truncate the current run
                    r->ccid = newccid++;
                    int x = (gridj_start+1)*splitsize;
                    r->x2 = x-1;

                    //runs.touch(nruns+gridj_span-1);
                    if (nruns+gridj_span > image->runs_allocated) {
                        while (nruns+gridj_span > image->runs_allocated) image->runs_allocated <<= 1;
                        image->runs = (struct Run *) realloc(image->runs,
                                                             image->runs_allocated * sizeof(struct Run));
                    }
                    // append additional runs to the runs array
                    image->runs_count = nruns+gridj_span;
                    for(int gridj=gridj_start+1; gridj<gridj_end; gridj++)
                    {
                        struct Run* newrun = &image->runs[nruns++];
                        newrun->y = y;
                        newrun->x1 = x;
                        x += splitsize;
                        newrun->x2 = x-1;
                        newrun->ccid = newccid++;
                    }
                    // append last run to the run array
                    struct Run* newrun = &image->runs[nruns++];
                    newrun->y = y;
                    newrun->x1 = x;
                    newrun->x2 = x_end;
                    newrun->ccid = newccid++;
                }
            }
        }
    }
    // Recompute cc descriptors
    ccimage_make_ccs_from_ccids(image);
}


// -- Helps sorting cc
static int
top_edges_descending (const void *pa, const void *pb)
{
    if (((struct CC*) pa)->bb.ymax != ((struct CC*) pb)->bb.ymax)
        return (((struct CC*) pb)->bb.ymax - ((struct CC*) pa)->bb.ymax);
    if (((struct CC*) pa)->bb.xmin != ((struct CC*) pb)->bb.xmin)
        return (((struct CC*) pa)->bb.xmin - ((struct CC*) pb)->bb.xmin);
    return (((struct CC*) pa)->frun - ((struct CC*) pb)->frun);
}


// -- Helps sorting cc
static int
left_edges_ascending (const void *pa, const void *pb)
{
    if (((struct CC*) pa)->bb.xmin != ((struct CC*) pb)->bb.xmin)
        return (((struct CC*) pa)->bb.xmin - ((struct CC*) pb)->bb.xmin);
    if (((struct CC*) pb)->bb.ymax != ((struct CC*) pa)->bb.ymax)
        return (((struct CC*) pb)->bb.ymax - ((struct CC*) pa)->bb.ymax);
    return (((struct CC*) pa)->frun - ((struct CC*) pb)->frun);
}


// -- Helps sorting cc
static int
integer_ascending (const void *pa, const void *pb)
{
    return ( *(int*)pb - *(int*)pa );
}


// -- Sort ccs in approximate reading order
void
ccimage_sort_in_reading_order(struct CCImage* image)
{
    if (image->nregularccs<2) return;
    struct CC *ccarray = MDJVU_MALLOCV(struct CC, image->nregularccs);
    // Copy existing ccarray (but segregate special ccs)
    int ccid;
    for(ccid=0; ccid<image->nregularccs; ccid++)
        ccarray[ccid] = image->ccs[ccid];
    // Sort the ccarray list into top-to-bottom order.
    qsort (ccarray, image->nregularccs, sizeof(struct CC), top_edges_descending);
    // Subdivide the ccarray list roughly into text lines [LYB]
    // - Determine maximal top deviation
    int maxtopchange = image->width / 40;
    if (maxtopchange < 32)
        maxtopchange = 32;
    // - Loop until processing all ccs
    int ccno = 0;
    int *bottoms = MDJVU_MALLOCV(int, image->nregularccs);
    while (ccno < image->nregularccs)
    {
        // - Gather first line approximation
        int nccno;
        int sublist_top = ccarray[ccno].bb.ymax-1;
        int sublist_bottom = ccarray[ccno].bb.ymin;
        for (nccno=ccno; nccno < image->nregularccs; nccno++)
        {
            if (ccarray[nccno].bb.ymax-1 < sublist_bottom) break;
            if (ccarray[nccno].bb.ymax-1 < sublist_top - maxtopchange) break;
            int bottom = ccarray[nccno].bb.ymin;
            bottoms[nccno-ccno] = bottom;
            if (bottom < sublist_bottom)
                sublist_bottom = bottom;
        }
        // - If more than one candidate cc for the line
        if (nccno > ccno + 1)
        {
            // - Compute median bottom
            qsort(bottoms, nccno-ccno, sizeof(int), integer_ascending);
            int bottom = bottoms[ (nccno-ccno-1)/2 ];
            // - Compose final line
            for (nccno=ccno; nccno < image->nregularccs; nccno++)
                if (ccarray[nccno].bb.ymax-1 < bottom)
                    break;
            // - Sort final line
            qsort (ccarray+ccno, nccno-ccno, sizeof(struct CC), left_edges_ascending);
        }
        // - Next line
        ccno = nccno;
    }
    // Copy ccarray back and renumber the runs
    for(ccid=0; ccid<image->nregularccs; ccid++)
    {
        const struct CC* cc = &ccarray[ccid];
        image->ccs[ccid] = *cc;
        for(int r=cc->frun; r<cc->frun+cc->nrun; r++)
            image->runs[r].ccid = ccid;
    }
    // Free memory
    MDJVU_FREEV(bottoms);
    MDJVU_FREEV(ccarray);
}

// -- Creates a bitmap for a particular component
mdjvu_bitmap_t
ccimage_get_bitmap_for_cc(struct CCImage* image, const int ccid)
{
    const struct CC *cc = &image->ccs[ccid];
    const struct GRect *bb = &cc->bb;

    unsigned char ** data = mdjvu_create_2d_array(bb->xmax - bb->xmin, bb->ymax - bb->ymin);
    const struct Run *prun = & image->runs[(int)cc->frun];
    for (int i=0; i<cc->nrun; i++,prun++)
    {
        if (prun->y<bb->ymin || prun->y>=bb->ymax)
            return NULL;
        if (prun->x1<bb->xmin || prun->x2>=bb->xmax)
            return NULL;
        unsigned char *row = data[prun->y - bb->ymin];
        for (int x=prun->x1; x<=prun->x2; x++)
            row[x - bb->xmin] = 1;
    }

    mdjvu_bitmap_t bitmap = mdjvu_bitmap_create(bb->xmax - bb->xmin, bb->ymax - bb->ymin);
    mdjvu_bitmap_pack_all(bitmap, data);
    mdjvu_destroy_2d_array(data);
    return bitmap;
}


// -- Creates a JB2Image with the remaining components
mdjvu_image_t
ccimage_get_jb2image(struct CCImage* image)
{
    mdjvu_image_t result = mdjvu_image_create(image->width, image->height);
    mdjvu_image_enable_suspiciously_big_flags(result);
    mdjvu_image_enable_not_a_letter_flags(result);
    mdjvu_image_set_resolution(result, image->dpi);

    if (image->runs_count <= 0)
        return result;

    // Iterate over CCs
    for (int ccid=0; ccid<image->ccs_count; ccid++)
    {
        mdjvu_bitmap_t bitmap = ccimage_get_bitmap_for_cc(image, ccid);
        mdjvu_image_add_bitmap(result, bitmap);
        mdjvu_image_add_blit(result, image->ccs[ccid].bb.xmin,
                             image->ccs[ccid].bb.ymin, bitmap);
        mdjvu_image_set_suspiciously_big_flag(result, bitmap, ccid >= image->nregularccs);
        mdjvu_image_set_not_a_letter_flag(result, bitmap, ccid >= image->nregularccs);
    }
    // Return
    return result;
}



struct CCImage *
        ccimage_create(int width, int height, int dpi)
{
    struct CCImage * ccimage = MDJVU_MALLOC(struct CCImage);
    ccimage->height = height;
    ccimage->width = width;
    ccimage->dpi = dpi;
    ccimage->nregularccs = 0;

    ccimage->runs_allocated = 16;
    ccimage->runs = MDJVU_MALLOCV(struct Run, 16);
    ccimage->ccs_allocated = 16;
    ccimage->ccs = MDJVU_MALLOCV(struct CC, 16);
    ccimage->runs_count = ccimage->ccs_count = 0;

    dpi = MAX(200, MIN(900, dpi));
    ccimage->largesize = MIN( 500, MAX(64, dpi));
    ccimage->smallsize = MAX(2, dpi/150);
    ccimage->tinysize = MAX(0, dpi*dpi/20000 - 1);
    return ccimage;
}

void
ccimage_free(struct CCImage* image)
{
    MDJVU_FREEV(image->runs);
    MDJVU_FREEV(image->ccs);
    MDJVU_FREE(image);
}

mdjvu_image_t
mdjvu_split(mdjvu_bitmap_t bitmap, int32 dpi, mdjvu_split_options_t opt)
{
    int32 width = mdjvu_bitmap_get_width(bitmap);
    int32 height = mdjvu_bitmap_get_height(bitmap);
    struct CCImage* ccimage = ccimage_create(width, height, dpi);
    ccimage_add_bitmap_runs(ccimage, bitmap, 0, 0, 0);
    // Component analysis
    ccimage_make_ccids_by_analysis(ccimage); // obtain ccids
    ccimage_make_ccs_from_ccids(ccimage);    // compute cc descriptors

    /* we don't use cjb2's cleaning algorithm */
    //    if (opts.losslevel > 0)
    //    rimg.erase_tiny_ccs();       // clean
    ccimage_merge_and_split_ccs(ccimage);    // reorganize weird ccs
    ccimage_sort_in_reading_order(ccimage);  // sort cc descriptors

    mdjvu_image_t result = ccimage_get_jb2image(ccimage); // get ``raw'' mdjvu_image_t
    ccimage_free(ccimage);
    return result;
}
