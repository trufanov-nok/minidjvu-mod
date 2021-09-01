/*
 * classify.c - classifying patterns
 */


#include "../base/mdjvucfg.h"
#include <minidjvu-mod/minidjvu-mod.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>


/* Stuff for not using malloc in C++
 * (left here for DjVuLibre compatibility)
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


typedef struct MinidjvuClassifyOptions
{
    int classifier;
} MinidjvuClassifyOptions;

MDJVU_IMPLEMENT mdjvu_classify_options_t mdjvu_classify_options_create()
{
    mdjvu_classify_options_t opt = (mdjvu_classify_options_t)
        malloc(sizeof(struct MinidjvuClassifyOptions));
    mdjvu_init();
    opt->classifier = 1;
    return opt;
}

MDJVU_IMPLEMENT void mdjvu_classify_options_destroy(mdjvu_classify_options_t opt)
{
    free(opt);
}

MDJVU_IMPLEMENT int mdjvu_get_classifier(mdjvu_classify_options_t opt)
    {return opt->classifier;}
MDJVU_IMPLEMENT void mdjvu_set_classifier(mdjvu_classify_options_t opt, int v)
    {opt->classifier = v;}

/* Classes are single-linked lists with an additional pointer to the last node.
 * This is an class item.
 */
typedef struct ClassNode
{
    mdjvu_pattern_t ptr;
    int32 id;
    int32 pos;
    int32 dpi;
    struct ClassNode *next;        /* NULL if this node is the last one */
    struct ClassNode *global_next; /* next among all nodes to classify  */
    int32 tag;                     /* filled before the final dumping   */
} ClassNode;

/* Classes themselves are composed in double-linked list. */
typedef struct Class
{
    ClassNode *first, *last;
    struct Class *prev_class;
    struct Class *next_class;
    ClassNode * compare_start_trick;
    int32 count;
} Class;


typedef struct Classification
{
    Class *first_class;
    ClassNode *first_node, *last_node;
} Classification;

typedef struct PatternList
{
    mdjvu_pattern_t p;
    int32 id;
    int32 pos;
    int32 dpi;
    struct PatternList *prev;
    struct PatternList *next;
} PatternList;

/* Creates an empty class and links it to the list of classes. */
static Class *new_class(Classification *cl)
{
    Class *c = MALLOC(Class);
    c->first = c->last = NULL;
    c->prev_class = NULL;
    c->count = 0;
    c->next_class = cl->first_class;
    if (cl->first_class) cl->first_class->prev_class = c;
    cl->first_class = c;
    return c;
}

/* Unlinks a class and deletes it. Its nodes are not deleted. */
static void delete_class(Classification *cl, Class *c)
{
    Class *prev = c->prev_class, *next = c->next_class;

    if (prev)
        prev->next_class = next;
    else
        cl->first_class = next;

    if (next)
        next->prev_class = prev;
    FREE(c);
}

/* Creates a new node and adds it to the given class. */
static ClassNode *new_node(Classification *cl, Class *c, PatternList * pl)
{
    ClassNode *n = MALLOC(ClassNode);
    n->ptr = pl->p;
    n->id  = pl->id;
    n->pos = pl->pos;
    n->dpi = pl->dpi;
    n->next = NULL;
    if (c->last) c->last->next = n;
    c->last = n;
    if (!c->first) c->first = n;
    n->global_next = NULL;
    c->count++;

    if (cl->last_node)
        cl->last_node->global_next = n;
    else
        cl->first_node = n;

    cl->last_node = n;
    return n;
}

/* Merge two classes and delete one of them. */
static Class *merge(Classification *cl, Class *c1, Class *c2)
{
    if (!c1->first)
    {
        delete_class(cl, c1);
        return c2;
    }
    if (c2->first)
    {
        c1->last->next = c2->first;
        c1->last = c2->last;
        c1->count += c2->count;
    }
    delete_class(cl, c2);
    return c1;
}

/* Puts a tag on each node corresponding to its class. */
static unsigned put_tags(Classification *cl)
{
    int32 tag = 1;
    Class *c = cl->first_class;
    while (c)
    {
        ClassNode *n = c->first;
        while (n)
        {
            n->tag = tag;
            n = n->next;
        }
        c = c->next_class;
        tag++;
    }
    return tag - 1;
}

/* Deletes all classes; nodes are untouched. */
static void delete_all_classes(Classification *cl)
{
    Class *c = cl->first_class;
    while (c)
    {
        Class *t = c;
        c = c->next_class;
        FREE(t);
    }
}


typedef struct CachedResults
{
    unsigned char** cache;
    size_t size;
} CachedResults;

CachedResults new_cache(size_t size) {
	CachedResults c;
        c.size = size;
        c.cache = MALLOCV(unsigned char*, size);
        double mem_req = 0;
        for (int i = 0; i < size; i++) {
                mem_req += (size - i + 3) >> 2;
        }

        fprintf(stdout, "Allocating cache for %lu elements: %0.2f MiB .... ", size, mem_req / 1024 / 1024);

        unsigned char* buf = MALLOCV(unsigned char, mem_req);
        if (buf) {
            fprintf(stdout, "done\n", size);
        } else {
            fprintf(stdout, "ERROR!\n");
            exit(-1);
        }

        memset(buf, 0xFF, mem_req);

        for (int i = 0; i < size; i++) {
                c.cache[i] = buf;
                buf += (size - i + 3) >> 2;
        }

	return c;
}

unsigned char masks[4] = {0xFC, 0xF3, 0xCF, 0x3F};

static void set_cache(CachedResults* c, int a, int b, int val) {
    val++;
    if (a > b) { int t = a; a = b; b = t; }
    unsigned char* v = c->cache[a] + ((c->size - b - 1) >> 2);
    const int shift = b & 0x3;
    *v = (*v & masks[shift]) | ( val << (shift*2) );
}

static void set_cache_by_line(CachedResults* c, unsigned char* line, int a, int b, int val) {
    val++;
    if (a > b) {int t = a; a = b; b = t;}
    unsigned char* v = line + ((c->size - b - 1) >> 2);
    const int shift = b & 0x3;
    *v = (*v & masks[shift]) | ( val << (shift*2) );
}

static char get_cache_and_line(const CachedResults* c, int a, int b, unsigned char** line) {
    if (a > b) {int t = a; a = b; b = t;}
    *line = c->cache[a];
    char val = (*line)[(c->size - b - 1) >> 2];
    val >>= (b & 0x3)*2;
    return (val & 0x3) - 1;
}


static void delete_cache(CachedResults* c) {
        FREEV(c->cache[0]);
	FREEV(c->cache);
	c->cache = NULL;
	c->size = 0;
}

/* Compares p with nodes from c until a meaningful result. */
static int compare_to_class(ClassNode* o, ClassNode* start_from, mdjvu_matcher_options_t options, CachedResults* cache)
{
    int r = 0;
    ClassNode *n = start_from;
    int positive_matches = 0;

    while(n)
    {
        if (cache) {
            unsigned char* line;
            r = get_cache_and_line(cache, o->id, n->id, &line);
            if (r == 2) {
                r = mdjvu_match_patterns(o->ptr, n->ptr, n->dpi, options);
                set_cache_by_line(cache, line, o->id, n->id, r);
            }
        } else {
            r = mdjvu_match_patterns(o->ptr, n->ptr, n->dpi, options);
        }

        if (r == -1) { // definetely wrong class
            return 0;
        }

        positive_matches += (r==1);
        n = n->next;
    }

    // return 0 if comparison to all examples in class was "0 (unknown, but probably different)"
    return positive_matches ? 1 : 0;
}

static void classify(Classification *cl, PatternList * all_patterns, mdjvu_matcher_options_t options, CachedResults* cache)
{
    if (!all_patterns->p) return;

    while (all_patterns) {
        mdjvu_pattern_t cur = all_patterns->p;
        Class * c = new_class(cl);
        new_node(cl, c, all_patterns);

        PatternList * maybe_seed = all_patterns->next;
        while (maybe_seed) {
            int res = mdjvu_match_patterns(cur, maybe_seed->p, all_patterns->dpi, options);
            if (cache) {
                set_cache(cache, all_patterns->id, maybe_seed->id, res);
            }

            if (res == 1) {
                if (maybe_seed->prev) {
                    maybe_seed->prev->next = maybe_seed->next;
                }
                if (maybe_seed->next) {
                    maybe_seed->next->prev = maybe_seed->prev;
                }

                new_node(cl, c, maybe_seed);
            }
            maybe_seed = maybe_seed->next;
        }

        all_patterns = all_patterns->next;
    }

    Class * c = cl->first_class;
    int classifier_level = mdjvu_get_classifier(mdjvu_get_classify_options(options));


    while (c && c->next_class != NULL) {


        Class * nc = c->next_class;
        while (nc != NULL) {
            nc->compare_start_trick = c->first;
            nc = nc->next_class;
        }

        int changed;

        do {
            changed = 0;
            Class * next_c = c->next_class;
            Class * recheck_to_last_merged = NULL;
            Class * next_recheck_to_last_merged = NULL;
            while (next_c != recheck_to_last_merged) {
                Class * next_to_next_c = next_c->next_class; /* That's because c may be deleted in merging */
                ClassNode* n = c->count >= next_c->count ? next_c->compare_start_trick : next_c->first;
                ClassNode* n2 = c->count >= next_c->count ? next_c->first : next_c->compare_start_trick;
                char need_merge = 0;
                while (n) {
                    if (compare_to_class(n, n2, options, cache)) {
                        need_merge = 1;
                        break;
                    }
                    n = n->next;
                }

                if (classifier_level > 2 && !need_merge) {
                    n = c->count < next_c->count ? next_c->compare_start_trick : next_c->first;
                    n2 = c->count < next_c->count ? next_c->first : next_c->compare_start_trick;
                    while (n) {
                        if (compare_to_class(n, n2, options, cache)) {
                            need_merge = 1;
                            break;
                        }
                        n = n->next;
                    }
                }

                if (need_merge) {
                    next_recheck_to_last_merged = next_to_next_c;
                    merge(cl, c, next_c);
                    changed = 1;
                } else {
                    next_c->compare_start_trick = c->last;
                }

                next_c = next_to_next_c;
            }

            recheck_to_last_merged = next_recheck_to_last_merged;

        } while (classifier_level > 1 && changed);

        c = c->next_class; // next may be NULL if two last classes were merged
    }

}

static int32 get_tags_from_classification(int32 *r, int32 n, Classification *cl)
{
    int32 max_tag = put_tags(cl);
    ClassNode *node;

    delete_all_classes(cl);

    memset(r, 0, sizeof(int32) * n);
    node = cl->first_node;
    while (node)
    {
        ClassNode *t;
        r[node->pos] = node->tag;
        t = node;
        node = node->global_next;
        FREE(t);
    }

    return max_tag;
}

static void init_classification(Classification *c)
{
    c->first_class = NULL;
    c->first_node = c->last_node = NULL;
}

MDJVU_IMPLEMENT int32 mdjvu_classify_patterns
    (mdjvu_pattern_t *b, int32 *r, int32 n, int32 dpi,
     mdjvu_matcher_options_t options)
{
    if (!n) return 0;

    int32 i;
    Classification cl;
    init_classification(&cl);

    PatternList* pl = MALLOCV(PatternList, n);
    memset(pl, 0, sizeof(PatternList)*n);
    PatternList* head = pl;
    head->p = NULL;
    PatternList* tail = NULL;

    double allocated_mem_stat = 0;

    for (i = 0; i < n; i++) {
        if (b[i]) {
            head->p = b[i];
            head->id = head->pos = i;
            head->dpi = dpi;
            head->next = head+1;
            head->prev = tail;
            allocated_mem_stat += mdjvu_pattern_mem_size(head->p);
            tail = head++;
        }
    }

    fprintf(stdout, "Classifier allocated memory: %0.2f MiB\n", allocated_mem_stat / 1024 / 1024);

    if (tail) {
        tail->next = NULL;
    }

    int classifier_level = mdjvu_get_classifier(mdjvu_get_classify_options(options));
    if (classifier_level == 1) {
        classify(&cl, pl, options, NULL);
    } else {
        CachedResults cache = new_cache(n);
        classify(&cl, pl, options, &cache);
        delete_cache(&cache);
    }
    MDJVU_FREEV(pl);

    return get_tags_from_classification(r, n, &cl);
}


static void get_cheap_center(mdjvu_bitmap_t bitmap, int32 *cx, int32 *cy)
{
    *cx = mdjvu_bitmap_get_width(bitmap) / 2;
    *cy = mdjvu_bitmap_get_height(bitmap) / 2;
}


#ifndef NO_MINIDJVU

MDJVU_IMPLEMENT int32 mdjvu_classify_bitmaps
    (mdjvu_image_t image, int32 *result, mdjvu_matcher_options_t options,
        int centers_needed)
{
    int32 i, n = mdjvu_image_get_bitmap_count(image);
    int32 dpi = mdjvu_image_get_resolution(image);
    mdjvu_pattern_t *patterns = MALLOCV(mdjvu_pattern_t, n);
    int32 max_tag;

    fprintf(stdout,"Size of JB2 image in memory: %0.2f MiB\n", (double) mdjvu_image_get_bitmap_count(image) / 1024 / 1024);

    for (i = 0; i < n; i++)
    {
        mdjvu_bitmap_t bitmap = mdjvu_image_get_bitmap(image, i);
        if (mdjvu_image_get_not_a_letter_flag(image, bitmap))
            patterns[i] = NULL;
        else
            patterns[i] = mdjvu_pattern_create(options, bitmap);
    }

    max_tag = mdjvu_classify_patterns(patterns, result, n, dpi, options);

    if (centers_needed)
    {
        mdjvu_image_enable_centers(image);
        for (i = 0; i < n; i++)
        {
            int32 cx, cy;
            mdjvu_bitmap_t bitmap = mdjvu_image_get_bitmap(image, i);
            if (patterns[i])
                mdjvu_pattern_get_center(patterns[i], &cx, &cy);
            else
                get_cheap_center(bitmap, &cx, &cy);
            mdjvu_image_set_center(image, bitmap, cx, cy); 
        }
    }

    for (i = 0; i < n; i++)
        if (patterns[i]) mdjvu_pattern_destroy(patterns[i]);
    FREEV(patterns);

    return max_tag;
}


/* ____________________________   multipage stuff   ________________________ */

/* FIXME: wrong dpi handling */
MDJVU_IMPLEMENT int32 mdjvu_multipage_classify_patterns
	(int32 npages, int32 total_patterns_count, const int32 *npatterns,
     mdjvu_pattern_t **patterns, int32 *result,
	 const int32 *dpi, mdjvu_matcher_options_t options,
     void (*report)(void *, int), void *param)
{
    /* a kluge for NULL patterns */
    /* FIXME: do it decently */
    if (!total_patterns_count) return 0;
    int32 page;
    int32 max_tag;

    Classification cl;
    init_classification(&cl);

    mdjvu_pattern_t* all_patterns = MALLOCV(mdjvu_pattern_t, total_patterns_count);
    PatternList* pl = MALLOCV(PatternList, total_patterns_count);
    PatternList* head = pl;
    head->p = NULL;
    PatternList* tail = NULL;

    int32 patterns_gathered = 0;
    int32 pl_num = 0;
    double allocated_mem_stat = 0;

    for (page = 0; page < npages; page++)
    {
        int32 n = npatterns[page];
        int32 d = dpi[page];
        mdjvu_pattern_t *p = patterns[page];

        int32 i;
        for (i = 0; i < n; i++) {
            if (*p) {
                head->p = *p;
                head->id = pl_num++;
                head->pos = patterns_gathered;
                head->dpi = d;
                head->next = head+1;
                head->prev = tail;
                allocated_mem_stat += mdjvu_pattern_mem_size(head->p);
                tail = head++;

            }
            all_patterns[patterns_gathered++] = *p;
            p++;
        }
        //        report(param, page);
    }

    fprintf(stdout, "Classifier allocated memory: %0.2f MiB\n", allocated_mem_stat / 1024 / 1024);

    if (tail) {
        tail->next = NULL;
    } else {
        MDJVU_FREEV(pl);
        MDJVU_FREEV(all_patterns);
        memset(result, 0, sizeof(int32) * total_patterns_count);
        return 0;
    }

    int classifier_level = mdjvu_get_classifier(mdjvu_get_classify_options(options));
    if (classifier_level == 1) {
        classify(&cl, pl, options, NULL);
    } else {
        CachedResults cache = new_cache(pl_num);
        classify(&cl, pl, options, &cache);
        delete_cache(&cache);
    }

    MDJVU_FREEV(pl);


    max_tag = get_tags_from_classification
            (result, total_patterns_count, &cl);

    MDJVU_FREEV(all_patterns);
    return max_tag;
}


MDJVU_IMPLEMENT int32 mdjvu_multipage_classify_bitmaps
    (int32 npages, int32 total_patterns_count, mdjvu_image_t *pages,
     int32 *result, mdjvu_matcher_options_t options,
     void (*report)(void *, int), void *param, int centers_needed)
{
    int32 max_tag, k, page;
    int32 *npatterns = (int32 *) malloc(npages * sizeof(int32));
    int32 *dpi = (int32 *) malloc(npages * sizeof(int32));
    mdjvu_pattern_t *patterns = (mdjvu_pattern_t *)
        malloc(total_patterns_count * sizeof(mdjvu_pattern_t));
    mdjvu_pattern_t **pointers = (mdjvu_pattern_t **)
        malloc(npages * sizeof(mdjvu_pattern_t *));

    double images_size_in_mem = 0;
    int32 patterns_created = 0;
    for (page = 0; page < npages; page++)
    {
        mdjvu_image_t current_image = pages[page];
        images_size_in_mem += mdjvu_image_mem_size(current_image);
        int32 c = npatterns[page] = mdjvu_image_get_bitmap_count(current_image);
        int32 i;
        dpi[page] = mdjvu_image_get_resolution(current_image);

        pointers[page] = patterns + patterns_created;
        for (i = 0; i < c; i++)
        {
            if (mdjvu_image_get_not_a_letter_flag(current_image, mdjvu_image_get_bitmap(current_image, i)))
                patterns[patterns_created++] = NULL;
            else
            {
                patterns[patterns_created++] = mdjvu_pattern_create(
                    options,
                    mdjvu_image_get_bitmap(current_image, i)
                );
            }
        }
    }

    fprintf(stdout,"Size of %u JB2 images in memory: %0.2f MiB\n", npages, images_size_in_mem / 1024 / 1024);

    max_tag = mdjvu_multipage_classify_patterns
        (npages, total_patterns_count, npatterns,
         pointers, result, dpi, options, report, param);

    if (centers_needed)
    {
        int32 patterns_processed = 0;
        
        for (page = 0; page < npages; page++)
        {
            mdjvu_image_t current_image = pages[page];
            int32 n = mdjvu_image_get_bitmap_count(current_image);
            int32 i;
            mdjvu_image_enable_centers(current_image);
            for (i = 0; i < n; i++)
            {
                int32 cx, cy;
                mdjvu_bitmap_t bitmap = mdjvu_image_get_bitmap(current_image, i);
                if (patterns[patterns_processed])
                    mdjvu_pattern_get_center(patterns[patterns_processed], &cx, &cy);
                else
                    get_cheap_center(bitmap, &cx, &cy);
                patterns_processed++;
                mdjvu_image_set_center(current_image, bitmap, cx, cy); 
            }
        }
    }

    for (k = 0; k < total_patterns_count; k++)
    {
        if (patterns[k])
            mdjvu_pattern_destroy(patterns[k]);
    }
    free(patterns);
    free(pointers);
    free(npatterns);
    free(dpi);

    return max_tag;
}


MDJVU_IMPLEMENT void mdjvu_multipage_get_dictionary_flags
   (int32 n,
	const int32 *npatterns,
    int32 max_tag,
	const int32 *tags,
    unsigned char *dictionary_flags)
{
    int32 page_number;
    int32 *first_page_met = (int32 *) malloc((max_tag + 1) * sizeof(int32));
    int32 i, total_bitmaps_passed = 0;
    memset(dictionary_flags, 0, max_tag + 1);
    for (i = 0; i <= max_tag; i++) first_page_met[i] = -1;

    for (page_number = 0; page_number < n; page_number++)
    {
        int32 bitmap_count = npatterns[page_number];

        for (i = 0; i < bitmap_count; i++)
        {
            int32 tag = tags[total_bitmaps_passed++];
            if (!tag) continue; /* skip non-substitutable bitmaps */

            if (first_page_met[tag] == -1)
            {
                /* never met this tag before */
                first_page_met[tag] = page_number;
            }
            else if (first_page_met[tag] != page_number)
            {
                /* met this tag on another page */
                dictionary_flags[tag] = 1;
            }
        }
    }

    free(first_page_met);
}


#endif /* NO_MINIDJVU */
