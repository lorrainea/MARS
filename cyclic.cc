/******************************************************************************
*                                                                             *
*  -------------------------- cyclic.c / cyclic.o -------------------------   *
*                                                                             *
******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "RestrictedLevenshtein.h"
#include "mars.h"
#include <sys/time.h>
#include <time.h>
#include <math.h> 

#define min(_x_,_y_) ( (_x_) < (_y_) ? (_x_) : (_y_)  )
#define max(_x_,_y_) ( (_x_) > (_y_) ? (_x_) : (_y_)  )
#define CompareCharacters( _x_ , _y_ ) ( (_x_) != (_y_) )

void cyclic(unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, float *gamma, int time_rep, int norm, unsigned int * rotation, unsigned int * distance)
{

	int i, cyc, row_length;
	float distance2, thresh , dummy_dist;
	int *SlopeToMinimumPath ;
	Limits Limit ;
	Path path , *BestPath , *ptr;

	thresh = 0.01*min(gamma[0], min(gamma[1], gamma[2]) ) ;
	row_length = 2*length1 + 1;

	/* Memory allocation */ 

	AllocateMemoryPath(&(Limit.Left), row_length) ;
	AllocateMemoryPath(&(Limit.Right), row_length) ;
	SlopeToMinimumPath = ( int * ) calloc( (length2+1)*(2*length1+1) , sizeof(int) ) ;

	BestPath = (Path *) malloc( (length1+1) * sizeof(BestPath[0])) ;
	for( i = 0, ptr = BestPath; i <= length1; ++i, ptr++)
		AllocateMemoryPath((ptr), row_length)  ; 
	AllocateMemoryPath(&path, row_length) ;

	/* Initializing Limits for standar Lev. distance */

	LimitInitialize(&Limit, length1, length2) ;

	/* Run bb method */
		
	int bb_type = 1 ;
	distance2 = RestrictedLevenshtein(0, pattern1, pattern2, length1, length2, Limit, BestPath, gamma, SlopeToMinimumPath) ;
	dummy_dist = bb(distance2, length1, length2, pattern1, pattern2, BestPath, &path,  bb_type , gamma, &cyc, rotation, distance, SlopeToMinimumPath ) ;

	/* Free memory */

	for( i = 0; i <= length1; ++i) 
		FreeMemoryPath(BestPath+i) ; 
	free(BestPath) ;
	FreeMemoryPath(&(Limit.Left)) ;
	FreeMemoryPath(&(Limit.Right)) ;
	free(SlopeToMinimumPath) ;
	FreeMemoryPath(&path) ;  

}

void LimitInitialize(Limits *Limit, int length1, int length2) 
{
	/*  We define the paths defining the limits for computing Levensh. distance   */

	Limit->Left.Minimum[0] = 0 ;
	Limit->Left.Maximum[0] = length2 ;
	Limit->Right.Minimum[length1] = 0 ;
	Limit->Right.Maximum[length1] = length2 ;
	for (int i = 1; i <= length1 ; ++i)
	{
		Limit->Left.Minimum[i] = length2 ;
		Limit->Left.Maximum[i] = length2 ;
		Limit->Right.Minimum[i-1] = 0 ;
		Limit->Right.Maximum[i-1] = 0 ; 
	}
}


void CopyPath(Path *Original, Path *Copy, int length)
{
	/* This procedure just copies paths from one pointer to another */
	for(int i =0; i< length; ++i)
	{
		Copy->Minimum[i] = Original->Minimum[i] ;
		Copy->Maximum[i] = Original->Maximum[i] ;
	}
}

void AllocateMemoryPath(Path *path, int size)
{
	path->Minimum  =   (int *) malloc(size * sizeof(int) ) ;
	path->Maximum  =   (int *) malloc(size * sizeof(int) ) ;
	return;
}

void FreeMemoryPath(Path *path)
{
	free(path->Minimum) ;
	free(path->Maximum) ;
}
