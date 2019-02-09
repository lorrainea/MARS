/******************************************************************************
*                                                                             *
* Programmer : Guillermo Peris                                                *
* Version : Febrero, 2001                                                     *
* Use : This procedure computes Levenshtein distance restricted to            *
*       some limits, in order to use it for the Maes algorithm cycle.         *
******************************************************************************/

#include<string.h>
#include<stdlib.h>
#include "RestrictedLevenshtein.h"
#include "mars.h"
#define CompareCharacters( _x_ , _y_ ) ( (_x_) != (_y_) )

float RestrictedLevenshtein(int FirstCharacter, unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, Limits Limit , Path *path, int *SlopeToMinimumPath)
{
	long long column, row, min_row, char1 ;
	unsigned long long ShortestPathWeight[length2+1][2] ;
	unsigned long long Delete, Substitute, Insert, min , path_dist;

	/* Initializing matrices */

	SlopeToMinimumPath[Limit.Right.Minimum[FirstCharacter]*(2*length1+1)+FirstCharacter]=-1;
	for(row = Limit.Right.Minimum[FirstCharacter] + 1 ; row <= Limit.Left.Maximum[FirstCharacter]; ++row) 
	{
		unsigned long long in= row*(2*length1+1)+FirstCharacter;
                if( in < (length2+1)*(2*length1+1) )
                        SlopeToMinimumPath[in] =  2 ;
        }


	/* Boundary conditions */

	for (row = 0 ; row <= Limit.Left.Maximum[FirstCharacter] ; ++row ) 
        	ShortestPathWeight[row][(FirstCharacter)%2] = row * INS ;


	/* Cycle for path under left limit path and over right path */

	for(column = FirstCharacter+1; column <= FirstCharacter+length1 ; ++column ) 
	{
	        if( Limit.Right.Minimum[column] == 0 ) 
		{
		 	ShortestPathWeight[0][column%2]  = (column - FirstCharacter) * INS ;
		      	SlopeToMinimumPath[column] = 0 ;
		      	min_row = 1;
        	}
        	else  min_row =  Limit.Right.Minimum[column];
        
		for ( row = min_row ; row <= Limit.Left.Maximum[column] ; ++row )
                {
                	min = INFINITE ;
                	
			/* Insertion */
                	if( row > Limit.Right.Minimum[column] )  
			{
                     		Insert = ShortestPathWeight[row - 1][column%2] + DEL ; 
                     		if(Insert < min) 
				{
                          		min = Insert ;
					unsigned long long in2= row*(2*length1+1)+column;

                                        if( in2 < (length2+1)*(2*length1+1) )
                                                SlopeToMinimumPath[in2] =  2 ;
                                     
                          	}
                     	}

                	/* Substitution  */
                	if( (row <= Limit.Left.Minimum[column]) && (  row > Limit.Right.Maximum[column-1]) ) 
			{
          	        	char1 = column ;
          	                if(char1 > length1) 
                                char1 = char1  -length1 ;
                         	Substitute = ShortestPathWeight[row - 1][(column-1)%2] ;
                          
				if( CompareCharacters(pattern1[char1-1],pattern2[row-1]) != 0 )
                                	Substitute = Substitute + SUB ;
                          	else
                                	Substitute = Substitute + MAT ;
                          	if(Substitute < min) 
				{
              		                min = Substitute ;
              		                unsigned long long in3= row*(2*length1+1)+column;
                                        if( in3 < (length2+1)*(2*length1+1) )
                                                SlopeToMinimumPath[in3] =  1 ;
                                }
			}
		
		        /*   Deletion   */
		        if( row <= Limit.Left.Maximum[column-1]) 
			{
		                Delete = ShortestPathWeight[row][(column-1)%2] + INS ;
		                if(Delete < min) 
				{
		                	min =Delete ;
		                      	unsigned long long in4 = row*(2*length1+1)+column;

                                        if( in4 < (length2+1)*(2*length1+1) )
                                                SlopeToMinimumPath[in4] =  0 ;
		                }
		        }
		        ShortestPathWeight[row][column%2]  = min ;
		}
	}

	/* Deriving Best Path */

	row--; 
	column--;
	
	while(column >= FirstCharacter) 
	{
		path->Maximum[column] = row ;
		while(SlopeToMinimumPath[row*(2*length1+1)+column] == 2) 
			row-- ;
		path->Minimum[column] = row ;
		if(SlopeToMinimumPath[row*(2*length1+1)+column] == 1) 
		{
		      row-- ;
		      column-- ;
		}
		else column-- ;
	}

	for(column=0; column < FirstCharacter; ++column) 
	{
		path->Minimum[column] = 0 ;
		path->Maximum[column] = 0 ;
	}

	for(column=FirstCharacter+length1+1; column <= 2*length1; ++column) 
	{
		path->Minimum[column] = length2;
		path->Maximum[column] = length2;
	}

return ShortestPathWeight[length2][(FirstCharacter+length1)%2] ; 
}
                        
