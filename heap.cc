/******************************************************************************
*                                                                             *
* Programmer : Guillermo Peris                                                *
* Version : Febrero, 2001                                                     *
* Use : This procedure computes Levenshtein distance restricted to            *
*       some limits, in order to use it for the Maes algorithm cycle.         *
******************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<heap.h>

#define Queue_Priority(_i_) ( ( (heap->array[(_i_)]).priority )  ) 
#define Queue_State(_i_) ( ((heap->array[(_i_)]).element ) ) 
#define Parent(_i_) ( ( ((_i_) -1)/2)  )
#define LeftChild(_i_) ( (2*(_i_) + 1)  ) 
#define RightChild(_i_) ( (2*(_i_) + 2)  ) 


void HeapInit (Heap *heap, int elements) 
{
	heap->size = 0 ;
	heap->root = NULL ;
	heap->array = (HeapElmt *) malloc(elements * sizeof(heap->array[0]));
}


void HeapInsert(Heap *heap, HEAP_TYPE element, COST_TYPE priority)
{
        Queue_Priority(heap->size) = priority ;
        Queue_State(heap->size)    = element ;
        heap->size++ ;
        HeapUp(heap, heap->size - 1) ; 
}
	

COST_TYPE HeapMin(Heap *heap) 
{
	return heap->root->priority ;
}

HEAP_TYPE HeapExtract(Heap *heap, COST_TYPE *priority)
{
	HEAP_TYPE rootnode ;

	rootnode = Queue_State(0) ;
	*priority = Queue_Priority(0) ;
	heap->array[0] = heap->array[heap->size - 1] ;
	(heap->size)-- ;
	HeapDown(heap, 0) ;
return rootnode ;
}


void HeapUp(Heap *heap, int i)
{
	HeapElmt t ;

	while( ( i != 0 ) && (Queue_Priority(Parent(i)) > Queue_Priority(i))  )
  	{
	  t = heap->array[Parent(i)] ;
	  heap->array[Parent(i)] = heap->array[i] ;
	  heap->array[i] = t ;
	  i = Parent(i) ;
  	}

	if(i == 0) 
		heap->root = heap->array ;
}

void HeapDown(Heap *heap, int i)
{
	HeapElmt t ;
	int max_son ;
	int left, right ;

	left  = LeftChild(i) ;
	right = RightChild(i) ;

	if(left >= heap->size) 
		max_son = -1 ;
	else if(right >= heap->size)
		max_son = left ;
	else
		max_son = (Queue_Priority(left) < Queue_Priority(right) ) ? left : right ;

	if(  (max_son > 0) && (Queue_Priority(i) > Queue_Priority(max_son) ) )
	{
		t = heap->array[i] ;
		heap->array[i] = heap->array[max_son] ;
		heap->array[max_son] = t ; 
		HeapDown(heap, max_son) ;
	}
	if(i == 0) 
		heap->root = heap->array ;
}

void HeapDestroy(Heap *heap) 
{
	heap->size = 0 ;
	heap->root = NULL ;
	free(heap->array) ;
}
