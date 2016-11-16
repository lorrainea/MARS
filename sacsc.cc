/**
    CSC: Circular Sequence Comparison
    Copyright (C) 2015 Solon P. Pissis, Ahmad Retha, Fatima Vayani 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
**/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <assert.h>
#include <sys/time.h>

#include "mars.h"
#include "sacsc.h"
#include <divsufsort.h>                                           // include header for suffix sort
#include <sdsl/bit_vectors.hpp>					  // include header for bit vectors

using namespace sdsl;
using namespace std;

unsigned int LCParray ( unsigned char *text, INT n, INT * SA, INT * ISA, INT * LCP )
{										
	INT i=0, j=0;

	LCP[0] = 0;
	for ( i = 0; i < n; i++ ) // compute LCP[ISA[i]]
		if ( ISA[i] != 0 ) 
		{
			if ( i == 0) j = 0;
			else j = (LCP[ISA[i-1]] >= 2) ? LCP[ISA[i-1]]-1 : 0;
			while ( text[i+j] == text[SA[ISA[i]-1]+j] )
				j++;
			LCP[ISA[i]] = j;
		}

	return ( 1 );
}

unsigned int circular_sequence_comparison (  unsigned char ** seq, struct TSwitch  sw, TPOcc ** D_original, unsigned int num_seqs )
{

	int total_length = 0;

	for(int j = 0; j<num_seqs; j++ )
	{
		total_length = total_length + ( 2 * strlen( ( char * ) seq[j] ) );
	}
		
	unsigned char * all_seqs = ( unsigned char * ) calloc( ( total_length+1 ) , sizeof( unsigned char * ) );

	for(int j = 0; j<num_seqs; j++ )
	{
		strcat ( ( char * ) all_seqs, ( char * ) seq[j] );
		strcat ( ( char * ) all_seqs, ( char * ) seq[j] );
	}
	all_seqs[total_length] = '\0';

	INT * SA;
	INT * LCP;
	INT * invSA;

        /* Compute the suffix array */
        SA = ( INT * ) malloc( ( total_length ) * sizeof( INT ) );
        if( ( SA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for SA.\n" );
                return ( 0 );
        }

        if( divsufsort( all_seqs, SA,  total_length ) != 0 )
        {
                fprintf(stderr, " Error: SA computation failed.\n" );
                exit( EXIT_FAILURE );
        }

        /* Compute the inverse SA array */
        invSA = ( INT * ) calloc( total_length , sizeof( INT ) );
        if( ( invSA == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for invSA.\n" );
                return ( 0 );
        }

        for ( INT i = 0; i < total_length; i ++ )
        {
                invSA [SA[i]] = i;
        }

	LCP = ( INT * ) calloc  ( total_length, sizeof( INT ) );
        if( ( LCP == NULL) )
        {
                fprintf(stderr, " Error: Cannot allocate memory for LCP.\n" );
                return ( 0 );
        }

        /* Compute the LCP array */
        if( LCParray( all_seqs, total_length, SA, invSA, LCP ) != 1 )
        {
                fprintf(stderr, " Error: LCP computation failed.\n" );
                exit( EXIT_FAILURE );
        }


	/* Ranking of q-grams and creation of x' and y' */

	int q = sw . q;
	int sigma = 0;

	INT mmnn = total_length - q + 1; 
	INT * rank; // 
	rank = ( INT * ) calloc( ( mmnn ) , sizeof( INT ) );


	/* Here we rank the first q-gram in the suffix array */	
	if( SA[0] <= total_length - q )
		rank[SA[0]] = sigma;


	/* Loop through the LCP array to rank the rest q-grams in the suffix array */	

	for ( INT i = 1; i < total_length; i++ )
	{
		INT lcp = LCP[i];
		INT ii = SA[i];

		if( ii <= total_length - q )
		{
			if ( lcp < q ) 
			{
				sigma++;
			}

			rank[ii] = sigma;
		}	
		
	}

	int pos_xx = 0;
	int pos_y = 0;

	for(int i=0; i<num_seqs; i++ )
	{
		int m = strlen( ( char * ) seq[i] );
	
		INT mm = m + m - q + 1;
		
		INT b;
		if( sw . l == 0 )
			b = (int) ( m / sqrt(m) );
		else b = (int) ( m / sw . l );

		INT * xx = ( INT * ) calloc( ( mm + 1 ) , sizeof( INT ) );


		memcpy( &xx[0], &rank[pos_xx], mm * sizeof( int ) ); 
		xx[mm] = '\0';

		int pos_y = 0;
		for(int j = 0; j<num_seqs; j++)
		{

			if( j == i )
			{
				pos_y =  pos_y + ( 2 * m ); // n = m
				continue;
			}
			
			
			int n = strlen( ( char * ) seq[j] );

			INT nn = n - q + 1; 

			INT * y = ( INT * ) calloc( ( nn + 1 ) , sizeof( INT ) );

			memcpy( &y[0], &rank[pos_y], nn * sizeof( int ));
			y[nn] = '\0';

			pos_y =  pos_y + ( 2 * n );
			
			/* Partitioning x' and y' as evenly as possible */
			INT * xind; 					//this is the starting position of the fragment
			INT * xmf; 					//this is the number of q-grams in the fragment
			xind = ( INT * ) calloc ( b, sizeof ( INT ) );
			xmf = ( INT * ) calloc ( b, sizeof ( INT ) );

			for ( INT j = 0; j < b; j++ )	
				partitioning ( 0, j, b, m - q + 1, xmf, xind );

			INT * yind; 					//this is the starting position of the fragment
			INT * ymf; 					//this is the number of q-grams in the fragment
			yind = ( INT * ) calloc ( b, sizeof ( INT ) );
			ymf = ( INT * ) calloc ( b, sizeof ( INT ) );

			for ( INT j = 0; j < b; j++ )	
				partitioning ( 0, j, b, nn, ymf, yind );
		
			#if 0
			for ( INT i = 0; i < b; i++ )
			{
				fprintf ( stderr, "(%d %d) ", xind[i], xmf[i] );
			}
			fprintf ( stderr, "\n" );
			#endif

			/* Allocate the diff vector */
			INT ** diff;
			diff = ( INT ** ) calloc ( b, sizeof ( INT * ) );
			for ( INT i = 0; i < b; i++ )	diff[i] = ( INT * ) calloc ( sigma + 1, sizeof ( INT ) );

			/* Step 1: Create diff, pvy, and D_0 */
			INT * D;
			D = ( INT * ) calloc ( b, sizeof ( INT * ) );
			for ( INT i = 0; i < b; i++ )
			{
				for ( INT j = yind[i]; j < yind[i] + ymf[i]; j++ )
				{
					diff[i][y[j]]++;	
					D[i]++;
				}
			}	
			//fprintf ( stderr, "D0 = %d\n", D[0] ); getchar();

			/* Step 2: Compute the distances for position 0 */
			int min_dist = 0;
			for ( INT i = 0; i < b; i++ )	//first window
			{	
				for ( INT j = xind[i]; j < xind[i] + xmf[i]; j++ )
				{
					diff[i][xx[j]]--;
					if ( diff[i][xx[j]] >= 0 )
					{
						D[i]--;
					}
					else
					{
						D[i]++;
					}
				}
				min_dist += D[i];
			}
			//fprintf ( stderr, "D0 = %d\n", D[0] ); getchar();

			/* Step 3: Compute the rest of the distances */
			int rot = 0;
			for ( INT i = 1; i < m; i++ )	//all the rest windows
			{
				int dist = 0;
				for ( INT j = 0; j < b; j++ )
				{
					diff[j][xx[i - 1 + xind[j]]]++; //letter out

					//For the letter we take out
					if ( diff[j][xx[i - 1 + xind[j]]] <= 0 )	
					{
						D[j]--;
					}
					else
					{
						D[j]++;
					}

					diff[j][xx[i - 1 + xind[j] + xmf[j]]]--; //letter in

					//For the letter we add in
					if ( diff[j][xx[i - 1 + xind[j] + xmf[j]]] < 0 )	
					{
						D[j]++;
					}
					else
					{							
						D[j]--;
					}
					dist += D[j];
				}
				//fprintf ( stderr, "dist = %d\n", dist );
				if ( dist < min_dist )
				{
					rot = i;
					min_dist = dist;
				}
			}
	
			D_original[i][j].err = ( unsigned int ) min_dist; 
			D_original[i][j].rot = ( unsigned int ) rot;

			/* De-allocate the memory */	
			free ( D );
			free ( y );	
			free ( xind );
			free ( xmf );	
			free ( yind );
			free ( ymf );	
			for ( INT i = 0; i < b; i++ )	
				free ( diff[i] );
			free ( diff );

		}
		pos_xx = pos_xx + ( 2 * m );

		free ( xx );
	}

	free ( invSA );
	free ( SA );
	free ( rank );
	free ( LCP );
	free ( all_seqs );	
		
	return ( 1 );
}

void partitioning ( INT i, INT j, INT f, INT m, INT * mf, INT * ind )
{
    	INT modulo = m % f;
    	double nf = ( double ) ( m ) / f;
    	INT first;
    	INT last;
	if ( j < modulo )
	{
		first = j * ( ceil( nf ) );
		last = first + ceil( nf ) - 1;
	}
	else
	{
		first = j * ( floor( nf ) ) + modulo;
		last = first + floor( nf ) - 1;
	}
	ind[j + i * f] = first;
	mf[j + i * f] = last - first + 1;
}
