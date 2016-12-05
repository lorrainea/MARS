/******************************************************************************
*                                                                             *
* Programmer : Guillermo Peris                                                *
* Version : Febrero, 2001                                                     *
* Use : This procedure computes Levenshtein distance restricted to            *
*       some limits, in order to use it for the Maes algorithm cycle.         *
******************************************************************************/

#ifndef _HEAP_H_INCLUDED

#define HEAP_TYPE Range           /* Ojo con los formatos en printf y scanf */
#define COST_TYPE float           

typedef struct Range_ {
        int left  ;             /* Left  limit of range                   */
        int right ;             /* Right limit of range                   */
        } Range ;

typedef struct HeapElmt_ {
        HEAP_TYPE  element;     /* element of the graph                   */
        COST_TYPE  priority;    /* Priority in queue                      */
        } HeapElmt ;

typedef struct Heap_ {
        int      size ;         /* Size of the heap                       */
        HeapElmt *root ;        /* Root node of the heap                  */
        HeapElmt *array ;       /* Array representing the heap            */
        } Heap ;

/* struct definition for range in branch-and-bound                        */


extern void HeapInit (Heap *heap, int elements) ;
extern void HeapInsert(Heap *heap, HEAP_TYPE element, COST_TYPE priority);
extern COST_TYPE HeapMin(Heap *heap) ;
extern HEAP_TYPE HeapExtract(Heap *heap, COST_TYPE *priority) ;
extern void HeapDestroy(Heap *heap) ;
void HeapUp(Heap *heap, int i) ;
void HeapDown(Heap *heap, int i) ;

#define HeapSize(_pqueue_) ( (_pqueue_)->size )  

#define _HEAP_H_INCLUDED

#endif
