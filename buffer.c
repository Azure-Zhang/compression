// ------------------------------------------------------------------
//   buffer.c
//   Copyright (C) 2019-2020 Divon Lan <divon@genozip.com>
//   Please see terms and conditions in the files LICENSE.non-commercial.txt and LICENSE.commercial.txt

// memory management - when running the same code by the same thread for another variant block - we reuse
// the previous variant's block memory. this way we save repetitive malloc/free cycles which might
// be very time consuming.

#include "genozip.h"
#include "profiler.h"
#include "buffer.h"
#include "vb.h"
#ifndef _MSC_VER // Microsoft compiler
#include <pthread.h>
#else
#include "compatability/visual_c_pthread.h"
#endif

#define DISPLAY_ALLOCS_AFTER 0 // display allocations, except the first X allocations. reallocs are always displayed

#define UNDERFLOW_TRAP 0x574F4C4652444E55ULL // "UNDRFLOW" - inserted at the begining of each memory block to detected underflows
#define OVERFLOW_TRAP  0x776F6C667265766FULL // "OVERFLOW" - inserted at the end of each memory block to detected overflows

static const unsigned overhead_size = 2*sizeof (uint64_t) + sizeof(uint16_t); // underflow, overflow and user counter

static pthread_mutex_t overlay_mutex; // used to thread-protect overlay counters
static uint64_t abandoned_mem_current = 0;
static uint64_t abandoned_mem_high_watermark = 0;

void buf_initialize()
{
    pthread_mutex_init (&overlay_mutex, NULL);

    vb_external_vb_initialize();
}

char *buf_human_readable_size (int64_t size, char *str /* out */)
{
    if      (size > (1LL << 40)) sprintf (str, "%3.1lf TB", ((double)size) / (double)(1LL << 40));
    else if (size > (1LL << 30)) sprintf (str, "%3.1lf GB", ((double)size) / (double)(1LL << 30));
    else if (size > (1LL << 20)) sprintf (str, "%3.1lf MB", ((double)size) / (double)(1LL << 20));
    else if (size > (1LL << 10)) sprintf (str, "%3.1lf KB", ((double)size) / (double)(1LL << 10));
    else                         sprintf (str, "%3d B"    ,     (int)size)                       ;

    return str; // for convenience so caller can use in printf directly
}

char *buf_human_readable_uint (int64_t n, char *str /* out */)
{
    unsigned len = 0, orig_len=0;

    if (n==0) 
        str[0] = '0';
        
    else {
        char rev[50] = {}; // "initialize" to avoid compiler warning
        while (n) {
            if (orig_len && orig_len % 3 == 0) rev[len++] = ',';    
            rev[len++] = '0' + n % 10;
            orig_len++;
            n /= 10;
        }
        // now reverse it
        for (int i=0; i < len; i++) str[i] = rev[len-i-1];
    }

    str[len] = '\0'; // string terminator
    return str;
}

#define POINTER_STR_LEN 19
char *buf_display_pointer (const void *p, char *str /* POINTER_STR_LEN bytes allocated by caller*/)
{
#ifdef _MSC_VER
    sprintf (str, "0x%I64x", (uint64_t)p);
#else
    sprintf (str, "0x%"PRIx64, (uint64_t)p);
#endif
    return str;
}

// get string with buffer's metadata for debug message. this function is NOT thread-safe
char *buf_display (const Buffer *buf)
{
    static char str[200]; // NOT thread-safe

    char s1[POINTER_STR_LEN], s2[POINTER_STR_LEN];
    sprintf (str, "Buffer %s (%u): size=%u len=%u data=%s memory=%s",
             buf->name, buf->param, buf->size, buf->len, buf_display_pointer(buf->data, s1), buf_display_pointer(buf->memory,s2));
    return str;    
}

static inline bool buf_has_overflowed (const Buffer *buf)
{
    return *((uint64_t*)(buf->memory + buf->size + sizeof(uint64_t))) != OVERFLOW_TRAP;
}

static inline bool buf_has_underflowed (const Buffer *buf)
{
    return *(uint64_t*)buf->memory != UNDERFLOW_TRAP;
}

static bool buf_test_overflows_do (const VariantBlock *vb, bool test_all_if_underflow);
static void buf_test_overflows_all_other_vb(const VariantBlock *caller_vb)
{
    VariantBlockPool *vb_pool = vb_get_pool();

    fprintf (stderr, "Testing all other VBs:\n");
    for (int vb_i=-1; vb_i < (int)vb_pool->num_vbs; vb_i++) {

        VariantBlock *vb = (vb_i == -1) ? evb : &vb_pool->vb[vb_i]; 
        if (vb == caller_vb) continue; // skip caller's VB

        fprintf (stderr, "Testing vb.id=%d:\n", vb->id);
        buf_test_overflows_do (vb, false);
    }
}

static bool buf_test_overflows_do (const VariantBlock *vb, bool test_all_if_underflow)
{
    if (!vb) return false;

    const Buffer *buf_list = &vb->buffer_list;

    bool corruption = false;
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
        const Buffer *buf = ((Buffer **)buf_list->data)[buf_i];

        if (!buf) continue; // buf was 'buf_destroy'd

        char s1[POINTER_STR_LEN], s2[POINTER_STR_LEN], s3[POINTER_STR_LEN];
        if (buf->memory) {
            if (buf->type < BUF_UNALLOCATED || buf->type > BUF_PARTIAL_OVERLAY) {
                fprintf (stderr, "\n\nMemory corruption: buffer=%s (buf_i=%u): Corrupt Buffer structure OR invalid buffer pointer - invalid buf->type\n", buf_display_pointer (buf, s1), buf_i);
                corruption = true;
            }
            else if (!buf->name) {
                fprintf (stderr, "\n\nMemory corruption: buffer=%s (buf_i=%u): Corrupt Buffer structure - null name\n", buf_display_pointer (buf, s1), buf_i);
                corruption = true;
            }
            else if (buf->data && buf->data != buf->memory + sizeof(uint64_t)) {
                fprintf (stderr, 
                        buf_display_pointer(buf,s1), buf_display_pointer(buf->memory,s2), buf_display_pointer(buf->memory+buf->size+2*sizeof(uint64_t)-1,s3), buf->name, buf->param,
                         "\n\nMemory corruption: data!=memory+8: vb_id=%d allocating_vb_i=%u buf_i=%u buffer=%s memory=%s func=%s:%u : Corrupt Buffer structure - expecting data+8 == memory. name=%s:%u buf->data=%s\n", 
                         vb ? vb->id : 0, buf->vb_i, buf_i, buf_display_pointer(buf,s1), buf_display_pointer(buf->memory,s2), buf->func, buf->code_line, buf->name, buf->param, buf_display_pointer(buf->data,s3));
                corruption = true;
            }
            else if (buf_has_underflowed(buf)) {
                fprintf (stderr, 
                        "\n\nMemory corruption: Underflow: buffer: %s memory: %s-%s name: %s:%u. Allocated: %s:%u vb_i=%u buf_i=%u. Found in buf_list of vb.id=%d. Fence=%c%c%c%c%c%c%c%c\n",
                        buf_display_pointer(buf,s1), buf_display_pointer(buf->memory,s2), buf_display_pointer(buf->memory+buf->size+2*sizeof(uint64_t)-1,s3), buf->name, buf->param,
                        buf->func, buf->code_line, buf->vb_i, buf_i, vb ? vb->id : -999, 
                        buf->memory[0], buf->memory[1], buf->memory[2], buf->memory[3], buf->memory[4], buf->memory[5], buf->memory[6], buf->memory[7]);

                if (test_all_if_underflow) buf_test_overflows_all_other_vb (vb);

                corruption = true;
            }
            else if (buf_has_overflowed(buf)) {
                char *of = &buf->memory[buf->size + sizeof(uint64_t)];
                fprintf (stderr,
                        "\n\nMemory corruption: Overflow: buffer: %s memory: %s-%s name: %s:%u. Allocated: %s:%u vb_i=%u buf_i=%u. Found in buf_list of vb.id=%d. Fence=%c%c%c%c%c%c%c%c\n",
                        buf_display_pointer(buf,s1), buf_display_pointer(buf->memory,s2), buf_display_pointer(buf->memory+buf->size+2*sizeof(uint64_t)-1,s3), buf->name, buf->param,
                        buf->func, buf->code_line, buf->vb_i, buf_i, vb ? vb->id : -999, 
                        of[0], of[1], of[2], of[3], of[4], of[5], of[6], of[7]);
                
                if (test_all_if_underflow) buf_test_overflows_all_other_vb (vb);

                corruption = true;
            }
        }
    }
    return corruption;
}

void buf_test_overflows (const VariantBlock *vb)
{
    ASSERT0 (!buf_test_overflows_do (vb, true), "Aborting due to memory corruption");
}

typedef struct {
    const char *name;
    uint64_t bytes; 
    unsigned buffers;
} MemStats;

static int buf_stats_sort_by_bytes(const void *a, const void *b)  
{ 
    return ((MemStats*)a)->bytes < ((MemStats*)b)->bytes ? 1 : -1;
}

void buf_display_memory_usage (bool memory_full, unsigned max_threads, unsigned used_threads)
{
    #define MAX_MEMORY_STATS 100
    static MemStats stats[MAX_MEMORY_STATS]; // must be pre-allocated, because buf_display_memory_usage is called when malloc fails, so it cannot malloc
    unsigned num_stats = 0, num_buffers = 0;

    if (memory_full)
        fprintf (stderr, "\n\nError memory is full:\n");
    else
        fprintf (stderr, "\n-------------------------------------------------------------------------------------\n");

    VariantBlockPool *vb_pool = vb_get_pool ();

    for (int vb_i=-1; vb_i < (int)vb_pool->num_vbs; vb_i++) {

        Buffer *buf_list = (vb_i == -1) ? &evb->buffer_list
                                        : &vb_pool->vb[vb_i].buffer_list; // a pointer to a buffer, which contains an array of pointers to buffers of a single vb/non-vb

        if (!buf_list->len) continue; // no buffers allocated yet for this VB

        for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
    
            ASSERT (buf_list->memory, "Error: memory of buffer_list of vb_i=%u is not allocated", vb_i); // this should never happen

            Buffer *buf = ((Buffer **)buf_list->data)[buf_i];
            
            if (!buf || !buf->memory) continue; // exclude destroyed, not-yet-allocated, overlay buffers and buffers that were src in buf_move

            bool found = false;
            for (unsigned st_i=0; st_i < num_stats && !found; st_i++) {
                MemStats *st = &stats[st_i];

                if (!strcmp (st->name, buf->name)) {
                    st->buffers++;
                    st->bytes += buf->size + overhead_size;
                    found = true;
                }
            }

            if (!found) {
                stats[num_stats].name    = buf->name;
                stats[num_stats].bytes   = buf->size + overhead_size;
                stats[num_stats].buffers = 1;
                num_stats++;
                ASSERT (num_stats < MAX_MEMORY_STATS, "# memory stats exceeded %u, consider increasing MAX_MEMORY_STATS", MAX_MEMORY_STATS);
            }

            num_buffers++;
        }
    }

    // add data_lines that are calloced
    stats[num_stats].name    = "data_lines";
    stats[num_stats].bytes   = 0;
    stats[num_stats].buffers = 0;
    num_stats++;

    for (int vb_i=0; vb_i < (int)vb_pool->num_vbs; vb_i++) 
        if (vb_pool->vb[vb_i].num_data_lines_allocated) {
            stats[num_stats-1].bytes += vb_pool->vb[vb_i].num_data_lines_allocated * (command==ZIP ? sizeof (ZipDataLine) : sizeof (PizDataLine));
            stats[num_stats-1].buffers += 1;
        }

    // sort stats by bytes
    qsort (stats, num_stats, sizeof (MemStats), buf_stats_sort_by_bytes);

    uint64_t total_bytes=0;
    for (unsigned i=0; i< num_stats; i++) total_bytes += stats[i].bytes;

    char str[30];
    buf_human_readable_size (total_bytes, str);
    fprintf (stderr, "Total bytes: %s in %u buffers in %u buffer lists:\n", str, num_buffers, vb_pool->num_vbs);
    fprintf (stderr, "Compute threads: max_permitted=%u actually_used=%u\n", max_threads, used_threads);

    for (unsigned i=0; i < num_stats; i++) {
        buf_human_readable_size (stats[i].bytes, str);
        fprintf (stderr, "%-30s: %-8s (%4.1f%%) in %u buffers\n", stats[i].name, str, 100.0 * (float)stats[i].bytes / (float)total_bytes, stats[i].buffers);
    }
}

int64_t buf_vb_memory_consumption (const VariantBlock *vb)
{
    const Buffer *buf_list   = &vb->buffer_list;

    // memory of the structure itself
    int64_t vb_memory = sizeof (*vb);

    // small BUG: This doesn't work - at the time it is called, num_data_lines_allocated is still 0
    vb_memory += vb->num_data_lines_allocated * (command == ZIP ? sizeof (ZipDataLine) : sizeof (PizDataLine));

    // memory allocated outside of Buffer (direct calloc)
    if (vb->haplotype_sections_data) vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->genotype_sections_data)  vb_memory += vb->num_sample_blocks * sizeof (Buffer);
    if (vb->phase_sections_data)     vb_memory += vb->num_sample_blocks * sizeof (Buffer);

    // memory allocated via Buffer (the vast majority, usually)
    for (unsigned buf_i=0; buf_i < buf_list->len; buf_i++) {
        Buffer *buf = ((Buffer **)buf_list->data)[buf_i]; 
        
        if (!buf || !buf->memory) continue; // exclude destroyed or not-yet-allocated buffers
        
        vb_memory += buf->size;
    }

    return vb_memory;
}

void buf_add_to_buffer_list (VariantBlock *vb, Buffer *buf)
{
#define INITIAL_MAX_MEM_NUM_BUFFERS 10000 /* for files that have ht,gt,phase,variant,and line - the factor would be about 5.5 so there will be 1 realloc per vb, but most files don't */
    Buffer *bl = &vb->buffer_list;

    buf_alloc (vb, bl, MAX (INITIAL_MAX_MEM_NUM_BUFFERS, bl->len+1) * sizeof(Buffer *), 2, "buffer_list", vb->id);

    ((Buffer **)bl->data)[bl->len++] = buf;

    if (flag_debug_memory && vb->buffer_list.len > DISPLAY_ALLOCS_AFTER) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "Init: %s:%u: size=%u buffer=%s vb->id=%d buf_i=%u\n", buf->name, buf->param, buf->size, buf_display_pointer(buf,s), vb->id, vb->buffer_list.len-1);
    }
}

static void buf_init (VariantBlock *vb, Buffer *buf, unsigned size, unsigned old_size, 
                      const char *func, unsigned code_line, const char *name, unsigned param)
{
    if (!buf->memory) {
        buf_test_overflows(vb);
#ifdef DEBUG
        buf_display_memory_usage (true, 0, 0);
#endif
        ABORT ("Error: Out of memroy. %sDetails: failed to allocate %u bytes name=%s:%u in %s:%u", 
               (command==ZIP ? "Try running with a lower variant block size using --vblock. " : ""), 
                size + overhead_size, name, param, func, code_line);
    }

    buf->data        = buf->memory + sizeof (uint64_t);
    buf->size        = size;
    buf->overlayable = false;
    buf->vb_i        = vb->variant_block_i;
    buf->func        = func;
    buf->code_line   = code_line;

    if (name) {
        buf->name = name;
        buf->param = param;
    } 
    ASSERT0 (buf->name, "Error: buffer has no name");

    *(uint64_t *)buf->memory        = UNDERFLOW_TRAP;        // underflow protection
    *(uint64_t *)(buf->data + size) = OVERFLOW_TRAP;         // overflow prortection (underflow protection was copied with realloc)
    *(uint16_t *)(buf->data + size + sizeof (uint64_t)) = 1; // counter of buffers that use of this memory (0 or 1 main buffer + any number of overlays)
}

// allocates or enlarges buffer
// if it needs to enlarge a buffer fully overlaid by an overlay buffer - it abandons its memory (leaving it to
// the overlaid buffer) and allocates new memory
unsigned buf_alloc_do (VariantBlock *vb,
                       Buffer *buf, 
                       uint32_t requested_size,
                       double grow_at_least_factor, // IF we need to allocate or reallocate physical memory, we get this much more than requested
                       const char *func, unsigned code_line,
                       const char *name, unsigned param)      
{
    START_TIMER;

    if (!requested_size) return 0; // nothing to do

    // sanity checks
    ASSERT (buf->type == BUF_REGULAR || buf->type == BUF_UNALLOCATED, "Error: cannot buf_alloc an overlayed buffer. name=%s", buf->name ? buf->name : "");
    ASSERT0 (vb, "Error: null vb");

    // case 1: we have enough memory already
    if (requested_size <= buf->size) {

        if (!buf->data) buf_init (vb, buf, buf->size, buf->size, func, code_line, name, param);
            
        //buf->data = buf->memory + sizeof (uint64_t); // allocate if not already allocated
        goto finish;
    }

    // add an epsilon to avoid floating point mutliplication ending up slightly less that the integer
    grow_at_least_factor = MAX (1, grow_at_least_factor) + 0.000000001; 

    // grow us requested - rounding up to 64 bit boundary to avoid aliasing errors with the overflow indicator
    uint32_t new_size = (uint32_t)(requested_size * grow_at_least_factor + 7) & 0xfffffff8UL; 

    ASSERT (new_size >= requested_size, "Error: allocated to little memory: requested=%u, allocated=%u", requested_size, new_size); // floating point paranoia

    // case 2: we need to allocate memory - buffer is already allocated so copy over the data
    if (buf->memory) {

        unsigned old_size = buf->size;

        // special handing if we have an overlaying buffer
        if (buf->overlayable) {
            pthread_mutex_lock (&overlay_mutex);
            uint16_t *overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

            char *old_data = buf->data;
            uint32_t old_len = buf->len;

            // if there is currently an overlay buffer on top of our buffer - abandon the memory
            // (leave it to the overlay buffer(s) that will eventually free() it), and allocate fresh memory
            if (*overlay_count > 1) {

                abandoned_mem_current += buf->size;
                abandoned_mem_high_watermark = MAX (abandoned_mem_high_watermark, abandoned_mem_current);

                (*overlay_count)--; // overlaying buffers are now on their own - no regular buffer
                buf->memory = buf->data = NULL;
                buf->size = buf->len = 0;
                buf_alloc (vb, buf, new_size, 1, name, param); // recursive call - simple alloc
                
                // copy old data
                memcpy (buf->data, old_data, old_size);
                buf->len = old_len;
            }
            else {
                // buffer is overlayable - but no current overlayers - regular realloc - however,
                // still within mutex to prevent another thread from overlaying while we're at it
                buf->memory = (char *)buf_low_level_realloc (buf->memory, new_size + overhead_size, func, code_line);
                buf_init (vb, buf, new_size, old_size, func, code_line, name, param);
            }
            buf->overlayable = true;
            pthread_mutex_unlock (&overlay_mutex);
        }

        else { // non-overlayable buffer - regular realloc without mutex
            buf->memory = (char *)buf_low_level_realloc (buf->memory, new_size + overhead_size, func, code_line);
            buf_init (vb, buf, new_size, old_size, func, code_line, name, param);
        }
    }

    // case 3: we need to allocate memory - buffer is not yet allocated, so no need to copy data
    else {
        buf->memory = (char *)malloc (new_size + overhead_size);
        buf->type  = BUF_REGULAR;

        buf_init (vb, buf, new_size, 0, func, code_line, name, param);
        buf_add_to_buffer_list(vb, buf);
    }

finish:
    if (vb) COPY_TIMER (vb->profile.buf_alloc);
    return buf->size;
}

// an overlay buffer is a buffer using some of the memory of another buffer - it doesn't have its own memory
void buf_overlay_do (VariantBlock *vb, Buffer *overlaid_buf, Buffer *regular_buf, const Buffer *copy_from /* optional */, 
                     unsigned *regular_buf_offset, const char *func, uint32_t code_line, const char *name, unsigned param)
{
    bool full_overlay = !regular_buf_offset && !copy_from;

//fprintf (stderr, "Overlaying onto buffer old_name=%s:=%u new_name=%s:%u\n", 
//         overlaid_buf->name, overlaid_buf->param, regular_buf->name, regular_buf->param);      

    // if this buffer was used by a previous VB as a regular buffer - we need to "destroy" it first
    if (overlaid_buf->type == BUF_REGULAR && overlaid_buf->data == NULL && overlaid_buf->memory) {
        buf_low_level_free (overlaid_buf->memory, func, code_line);
        overlaid_buf->type = BUF_UNALLOCATED;
    }
    
    ASSERT (overlaid_buf->type == BUF_UNALLOCATED, "Error: cannot buf_overlay to a buffer already in use. overlaid_buf->name=%s", overlaid_buf->name ? overlaid_buf->name : "");
    ASSERT (regular_buf->type == BUF_REGULAR, "Error: regular_buf in buf_overlay must be a regular buffer. regular_buf->name=%s", regular_buf->name ? regular_buf->name : "");
    ASSERT (!full_overlay || regular_buf->overlayable, "Error: buf_overlay: only overlayble buffers can be fully overlaid. regular_buf->name=%s", regular_buf->name ? regular_buf->name : "");

    overlaid_buf->size        = 0;
    overlaid_buf->len         = copy_from ? copy_from->len : 0;
    overlaid_buf->type        = full_overlay ? BUF_FULL_OVERLAY : BUF_PARTIAL_OVERLAY;
    overlaid_buf->memory      = 0;
    overlaid_buf->overlayable = false;
    overlaid_buf->vb_i        = vb->variant_block_i;

    if (name) {
        overlaid_buf->name    = name;
        overlaid_buf->param   = param;
    }
    else {
        overlaid_buf->name    = regular_buf->name;
        overlaid_buf->param   = regular_buf->param;
    }

    if (!full_overlay) {
        overlaid_buf->data = regular_buf->data + (regular_buf_offset ? *regular_buf_offset : 0);

        if (copy_from && copy_from->len) {

            ASSERT ((regular_buf_offset ? *regular_buf_offset : 0) + copy_from->len <= regular_buf->size, 
                    "Error: buf_overlay exceeds the size of the regular buffer: offset=%u size=%u regular_buf.size=%u", 
                    *regular_buf_offset, copy_from->len, regular_buf->size);
        
            memcpy (overlaid_buf->data, copy_from->data, copy_from->len);

            // new data was copied to overlaid buffer - move the offset forward and update the len of the regular buffer
            // to enable to the next buffer to be overlaid subsequently
            if (regular_buf_offset) {
                *regular_buf_offset += overlaid_buf->len;
                //regular_buf->len = *regular_buf_offset;
                
            }
        }
    }

    // full buffer overlay - copy len too and update overlay counter
    else {
        pthread_mutex_lock (&overlay_mutex);

        overlaid_buf->size = regular_buf->size;
        overlaid_buf->len  = regular_buf->len;
        overlaid_buf->data = regular_buf->data;
        uint16_t *overlay_count = (uint16_t*)(regular_buf->data + regular_buf->size + sizeof(uint64_t));
        (*overlay_count)++; // counter of users of this memory

        pthread_mutex_unlock (&overlay_mutex);
    }
}

// free buffer - without freeing memory. A future buf_alloc of this buffer will reuse the memory if possible.
void buf_free_do (Buffer *buf, const char *func, uint32_t code_line) 
{
    uint16_t *overlay_count; // number of buffers (overlay and regular) sharing buf->memory

    switch (buf->type) {

        case BUF_UNALLOCATED:
            return; // nothing to do

        case BUF_REGULAR: 

            if (buf->overlayable) {

                pthread_mutex_lock (&overlay_mutex);
                overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));

                if (*overlay_count > 1) { // current overlays exist - abandon memory - leave it to the overlaid buffer(s) which will free() this memory when they're done with it
                    (*overlay_count)--;
             
                    abandoned_mem_current += buf->size;
                    abandoned_mem_high_watermark = MAX (abandoned_mem_high_watermark, abandoned_mem_current);

                    memset (buf, 0, sizeof (Buffer)); // make this buffer UNALLOCATED
                }
                // if no overlay exists then we just keep .memory and reuse it in future allocations

                pthread_mutex_unlock (&overlay_mutex);            
            }
            
            buf->data = NULL; 
            buf->len = 0;
            buf->overlayable = false;
            buf->vb_i = 0;
            
            // name, param, memory and size are not changed

            break;

        case BUF_FULL_OVERLAY:
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (uint16_t*)(buf->data + buf->size + sizeof(uint64_t));
            (*overlay_count)--;
            pthread_mutex_unlock (&overlay_mutex);            

            // we are the last user - we can free the memory now.
            // do this outside of the mutex - free is a system call and can take some time.
            // this is safe because if we ever observe *overlay_count==0, it means that no buffer has this memory,
            // therefore there is no possibility it would be subsequently overlayed between the test and the free().
            if (! (*overlay_count)) {
                buf_low_level_free (buf->data - sizeof(uint64_t), func, code_line); // the original buf->memory
                abandoned_mem_current -= buf->size;
            }
            
            // fall through

        case BUF_PARTIAL_OVERLAY:
            memset (buf, 0, sizeof (Buffer));
            break;

        default:
            ABORT0 ("Error: invalid buf->type");
    }
} 

// remove from buffer_list of this vb
void buf_remove_from_buffer_list (VariantBlock *vb, Buffer *buf)
{
    unsigned i=0; for (; i < vb->buffer_list.len; i++) 
        if (((Buffer **)vb->buffer_list.data)[i] == buf) {
            ((Buffer **)vb->buffer_list.data)[i] = NULL;
    
            if (flag_debug_memory) {
                char s[POINTER_STR_LEN];
                fprintf (stderr, "Destroy:%s:%u: buffer=%s vb->id=%d buf_i=%u\n", buf->name, buf->param, buf_display_pointer(buf,s), vb->id, i);
            }
            break;
        }

    // note: it is possible that the buffer is not found in the list if it is never allocated. that's fine.
}

void buf_destroy_do (VariantBlock *vb, Buffer *buf, const char *func, uint32_t code_line)
{
    if (!buf) return; // nothing to do

    buf_remove_from_buffer_list (vb, buf); // remove buffer regardless of memory - since the buffer itself might be free()d following the destroy

    if (buf->memory) {
    
        uint16_t overlay_count = 1;
        if (buf->overlayable) {
            pthread_mutex_lock (&overlay_mutex);
            overlay_count = (*(uint16_t*)(buf->data + buf->size + sizeof(uint64_t)));
            pthread_mutex_unlock (&overlay_mutex);            
        }

        ASSERT (overlay_count==1, "Error: cannot destroy buffer %s because it is currently overlaid", buf->name);

        buf_low_level_free (buf->memory, func, code_line);
    }

    memset (buf, 0, sizeof (Buffer)); // reset to factory defaults
}

void buf_copy (VariantBlock *vb, Buffer *dst, const Buffer *src, 
               unsigned bytes_per_entry, // how many bytes are counted by a unit of .len
               unsigned src_start_entry, unsigned max_entries,  // if 0 copies the entire buffer 
               const char *name, unsigned param)
{
    ASSERT0 (src->data, "Error in buf_copy: src->data is NULL");
    
    ASSERT (!max_entries || src_start_entry < src->len, 
            "Error buf_copy of name=%s:%u: src_start_entry=%u is larger than src->len=%u", src->name, src->param, src_start_entry, src->len);

    unsigned num_entries = max_entries ? MIN (max_entries, src->len - src_start_entry) : src->len - src_start_entry;
    if (!bytes_per_entry) bytes_per_entry=1;
    
    buf_alloc(vb, dst, num_entries * bytes_per_entry, 1, 
              name ? name : src->name, name ? param : src->param); // use realloc rather than malloc to allocate exact size

    memcpy (dst->data, &src->data[src_start_entry * bytes_per_entry], num_entries * bytes_per_entry);

    dst->len = num_entries;  
}   

// moves all the data from one buffer to another, leaving the source buffer unallocated
void buf_move (VariantBlock *vb, Buffer *dst, Buffer *src)
{
    ASSERT (dst->type == BUF_UNALLOCATED, "Error: attempt to move to an already-allocated src: name=%s:%u dst: name=%s:%u",
            src->name, src->param, dst->name, dst->param);

    buf_remove_from_buffer_list (vb, src); // must be before the memset - otherwise memory will be 0 and it won't work

    memcpy (dst, src, sizeof(Buffer));
    memset (src, 0, sizeof(Buffer));
    
    if (dst->memory) buf_add_to_buffer_list (vb, dst);
}

void buf_add_string (VariantBlock *vb, Buffer *buf, const char *str) 
{ 
    unsigned len = strlen (str); 
    buf_alloc (vb, buf, MAX (1000, buf->len + len + 1), 2, "string_buf", 0);
    buf_add (buf, str, len);
    buf->data[buf->len] = '\0'; // string terminator without increasing buf->len
}

void buf_print (Buffer *buf, bool add_newline)
{
    for (unsigned i=0; i < buf->len; i++) 
        putchar (buf->data[i]);  // safer than printf %.*s ?

    if (add_newline) putchar ('\n');
}

void buf_low_level_free (void *p, const char *func, uint32_t code_line)
{
    if (flag_debug_memory) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "Memory freed by free(): %s %s:%u\n", buf_display_pointer (p, s), func, code_line);
    }

    free (p);
}

void *buf_low_level_realloc (void *p, size_t size, const char *func, uint32_t code_line)
{
    void *new = realloc (p, size);

    if (flag_debug_memory && new != p) {
        char s[POINTER_STR_LEN];
        fprintf (stderr, "Memory freed by realloc(): %s %s:%u\n", buf_display_pointer (p, s), func, code_line);
    }

    return new;
}
