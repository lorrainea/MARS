#include <stdio.h>
#include <stdlib.h> 
#include <iostream>
#include "heap.h"
#include "RestrictedLevenshtein.h"

#define min(_x_,_y_) ( (_x_) < (_y_) ? (_x_) : (_y_) )
#define max(_x_,_y_) ( (_x_) > (_y_) ? (_x_) : (_y_) )

using namespace std;

float bb(float distance2, int length1, int length2, unsigned char *pattern1, unsigned char *pattern2, Path *BestPath, Path *path, int bb_type, float *gamma, int *cyc, unsigned int * rotation, unsigned int * distance, int * SlopeToMinimumPath)
{
	int rotation2 , cc_bound, c_bound , r_bound, l_bound, i, lll, rrr, ccc;
	float null, min_distance, RotatedDistance[length1 + 1], lower_bound ;
	Heap partition;
	Range range;
	Limits Limit ;
	int row_length = 2*length1+1 ;
	float external_bound=INFINITE ;

	/* Memory allocation */

	min_distance = external_bound ;
	AllocateMemoryPath(&(Limit.Left), row_length) ;
	AllocateMemoryPath(&(Limit.Right), row_length) ;

	/* Initialize RotatedDistance */

	RotatedDistance[0] = distance2 ;
	RotatedDistance[length1] = distance2 ; 

	/*  Translate BestPath[0] */

	for( i = 0; i< length1; ++i) 
	{
		BestPath[length1].Minimum[i] = 0 ;
		BestPath[length1].Maximum[i] = 0 ;
		BestPath[length1].Minimum[i+length1] = BestPath[0].Minimum[i]  ;
		BestPath[length1].Maximum[i+length1] = BestPath[0].Maximum[i]  ;
	}
	BestPath[length1].Minimum[2*length1] = BestPath[0].Minimum[length1]  ;
	BestPath[length1].Maximum[2*length1] = BestPath[0].Maximum[length1]  ;

	/* Heap initialization. Insert initial state in heap */

	range.left  = 0 ;
	range.right = length1   ;        
	if(min_distance > distance2) min_distance = distance2 ;
		rotation2 = 0 ; 

	lower_bound = BoundFunction(0, length1, RotatedDistance[0], RotatedDistance[length1], pattern1, pattern2, length1, length2, bb_type, gamma) ;

	HeapInit(&partition, length1+1) ; 
	HeapInsert(&partition, range, lower_bound ); 

	/* Begin b&b cycle */

	while ( (HeapSize(&partition) > 0) && (min_distance > HeapMin(&partition) ) ) 
	{
		range = HeapExtract(&partition, &null)  ; 
		l_bound = range.left ;
		r_bound = range.right ;
		c_bound = l_bound + (r_bound - l_bound) / 2 ; 

		CopyPath(BestPath + l_bound, &(Limit.Left), 2*length1 + 1 ) ;
		CopyPath(BestPath + r_bound, &(Limit.Right), 2*length1 + 1 ) ; 
		distance2 = RestrictedLevenshtein(c_bound, pattern1, pattern2, length1, length2, Limit, path, gamma, SlopeToMinimumPath) ;

		RotatedDistance[c_bound] = distance2 ;
		CopyPath(path, BestPath + c_bound, 2*length1 + 1 ) ;
		if(distance2 < min_distance) 
		{         
			min_distance = distance2 ;
			rotation2 = c_bound; 
		}

		/* Analisis de particion izquierda */
		lower_bound = BoundFunction(l_bound, c_bound, RotatedDistance[l_bound], RotatedDistance[c_bound], pattern1, pattern2, length1, length2, bb_type, gamma) ;
     		if(lower_bound < min_distance ) 
		{
			lll = l_bound ;
			ccc = c_bound ;
			
			if( ccc > lll +2 ) 
			{
				range.left = lll ;
				range.right = ccc ;
				HeapInsert(&partition, range , lower_bound);
			}
	     		else if( ccc == lll +2 ) 
			{
				CopyPath(BestPath + lll, &(Limit.Left), 2*length1 + 1 ) ;
				CopyPath(BestPath + ccc, &(Limit.Right), 2*length1 + 1 ) ; 
				cc_bound = lll + 1;
				distance2 = RestrictedLevenshtein(cc_bound, pattern1, pattern2, length1, length2, Limit, path, gamma, SlopeToMinimumPath) ;
			
				RotatedDistance[cc_bound] = distance2 ;
				CopyPath(path, BestPath + cc_bound, 2*length1 + 1 ) ;
				if(distance2 < min_distance) 
				{
					min_distance = distance2 ;
					rotation2 = cc_bound;
				}
			}
		}

		/* Analisis de particion derecha */
     		lower_bound = BoundFunction(c_bound, r_bound, RotatedDistance[c_bound], RotatedDistance[r_bound], pattern1, pattern2, length1, length2, bb_type, gamma) ;
		if(lower_bound < min_distance ) 
		{
			ccc = c_bound ;
			rrr = r_bound ;
			if( rrr > ccc +2 ) 
			{
				range.left  = ccc ;
				range.right = rrr ;
				HeapInsert(&partition, range , lower_bound);
	                }
	     		else if( rrr == ccc +2 ) 
			{
				CopyPath(BestPath + ccc, &(Limit.Left), 2*length1 + 1 ) ;
				CopyPath(BestPath + rrr, &(Limit.Right), 2*length1 + 1 ) ;
				cc_bound = ccc + 1;
				distance2 = RestrictedLevenshtein(cc_bound, pattern1, pattern2, length1, length2, Limit, path, gamma, SlopeToMinimumPath) ;

				RotatedDistance[cc_bound] = distance2 ;
				CopyPath(path, BestPath + cc_bound, 2*length1 + 1 ) ;

				if(distance2 < min_distance) 
				{
	                            min_distance = distance2 ;
	                            rotation2 = cc_bound;
				}
			}
		}
	}

	CopyPath( BestPath + rotation2, path, 2*length1 + 1 ) ;
	*cyc = rotation2 ;

	*rotation = rotation2;
	*distance = min_distance;

	/* Free memory */

	FreeMemoryPath(&(Limit.Left)) ;
	FreeMemoryPath(&(Limit.Right)) ;
	HeapDestroy(&partition) ;
	if(min_distance == external_bound)
		min_distance = -1.0;

return min_distance ;
}


float BoundFunction(int left, int right, float left_cost, float right_cost, unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, int bound_type, float *gamma)
{
	int i ;
	float distance_approx, distance, suma ;

	suma = gamma[1] + gamma[0];

	/* Calculo de bb1 */
	if(bound_type == 1) 
	{
      		distance_approx = (left_cost + right_cost)/2.0 + (float) (left - right)*suma/2.0 ;
      		if (distance_approx < 0.0) distance_approx = 0.0;  
			return distance_approx ;
	}
}

