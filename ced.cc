/**
    hCED: Heuristic Cyclic Edit Distance
    Copyright (C) 2016 Solon P. Pissis, Lorraine A. K. Ayad 

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

#include <seqan/align.h>
#include <stdio.h>
#include <stdlib.h>
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "mars.h"
#include "sacsc.h"
#include "ced.h"


#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))


using namespace std;
using namespace seqan;

int delta ( char a, char b, struct TSwitch sw )
 {
  	if ( a == DEL || b == DEL ) 
    	{
		return 0;
    	}
   	else
	{
		int matching_score = ( sw . matrix ? pro_delta( a, b ) : nuc_delta( a, b ) ) ;
		if ( matching_score == ERR )
			return 0;
		else return matching_score;
    	}
 }


unsigned int nw_ag_allocation( unsigned int m, unsigned int n, int * &I, int * &D, int ** &T )
{
	int i;

	if ( ( T = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	
	I = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) );
	
	D = ( int * ) calloc ( ( m + 1 ) , sizeof( int ) );

	return EXIT_SUCCESS;
}


unsigned int nw_ag ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch  sw, int * score, int * &I, int * &D, int ** &T)
{
	int i, j;
	int g = sw . O;
	int h = sw . E;
	int u, v, w;
        
   	for ( i = 0; i < m + 1; i++ )
	{
		D[i] = m * -1;
	}

	for ( j = 0; j < n + 1; j++ )
	{
		I[j] = n * -1;
	}

	T[0][0] = 0;
	if ( m > 0 )
		T[1][0] = g;
	for ( i = 2; i < m + 1; i++ )
	T[i][0] = T[i - 1][0] + h;
	if ( n > 0 )
	T[0][1] = g;
	for ( j = 2; j < n + 1; j++ )
	T[0][j] = T[0][j - 1] + h;

	for( i = 1; i < m + 1; i++ )
	{
        	for( j = 1; j < n + 1; j++ )
        	{
			D[i] = MAX2 ( D[i - 1] + h, T[i - 1][j] + g );
			u = D[i];

			I[j] = MAX2 ( I[j - 1] + h, T[i][j - 1] + g );
			v = I[j];

			w = T[i - 1][j - 1] + delta ( t[j - 1], p[i - 1], sw );

			T[i][j] = MAX3 ( w, u, v );
        	}
    	}

	( * score ) = T[m][n];
	
	return EXIT_SUCCESS;
}


unsigned int nw_allocation( unsigned int m, unsigned int n, int ** &T )
{

	if ( ( T = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: T could not be allocated!\n");
                return ( 0 );
        }
        for ( int i = 0; i < m + 1; i ++ )
        {
                if ( ( T[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: T could not be allocated!\n");
                        return ( 0 );
                }
        }
	
	return EXIT_SUCCESS;
}


unsigned int nw ( unsigned char * p, unsigned int m, unsigned char * t, unsigned int n, struct TSwitch  sw, int * score, int ** &T )
{
	init_substitution_score_tables ();
	int ins = sw . O;
	int del = sw . E;
	int i, j;
	
   	for ( i = 0; i < m + 1; i++ )
	{
		T[i][0] = del * i;
	}
	for ( j = 1; j < n + 1; j++ )
	{
		T[0][j] = ins * j;
	}

        for ( i = 1; i < m + 1; i++ )
	{
		for ( j = 1; j < n + 1; j++ )
			T[i][j] = MAX3(T[i][j-1] + ins, T[i-1][j] + del, T[i-1][j-1] + delta ( t[j - 1], p[i - 1], sw ));
	}

    	( * score ) = T[m][n];

	return EXIT_SUCCESS;
}

unsigned int sacsc_refinement ( unsigned char * x, unsigned char * xr, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance ) 
{
	unsigned int rot = *rotation;

	int * I;
	int * D;
	int ** T;

	unsigned int m = strlen ( ( char * ) x );
	unsigned int n = strlen ( ( char * ) y );

	unsigned char * X;
	unsigned char * Y;

	unsigned int sl = sw . P * ( sw . l ); //section length
	sl = MIN3 ( sl, m/2, n/2 );

	X = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );
	Y = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );

	memcpy ( &X[0], &xr[0], sl );
	for ( int i = 0; i < sl; i++ )
		X[sl + i] = DEL;
	memcpy ( &X[sl + sl], &xr[m - sl], sl );
	X[3 * sl] = '\0';
	
	memcpy ( &Y[0], &y[0], sl );
	for ( int i = 0; i < sl; i++ )
		Y[sl + i] = DEL;
	memcpy ( &Y[sl + sl], &y[n - sl], sl );
	Y[3 * sl] = '\0';

	unsigned int mm = sl + sl + sl;
	unsigned int nn = sl + sl + sl;

	int score = -INT_MAX;
	int max_score = score;
	unsigned int rrot = 0;
	unsigned char * Xr = ( unsigned char * ) calloc( ( 3 * sl + 1 ) , sizeof( unsigned char ) );

	if( sw . O != sw . E )
		nw_ag_allocation ( m, n, I, D, T );
	else 
		nw_allocation( m, n, T );

	for ( int i = 0; i < mm; i++ )
	{
		if ( i >= sl && i < 2 * sl )
			continue;

		memmove ( &Xr[0], &X[i], ( 3 * sl ) - i );
    		memmove ( &Xr[( 3 * sl ) - i], &X[0], i );
    		Xr[3*sl] = '\0';

		if( sw . O != sw . E )
			nw_ag ( Xr, mm , Y, nn, sw, &score, I, D, T );
		else
			nw ( Xr, mm, Y, nn, sw, &score, T );

		if ( score > max_score )
		{
			max_score = score;
			rrot = i;
		}	 
	}

	for ( int j = 0; j < m + 1; j ++ )
	{
		free ( T[j] );
	}

	free ( Xr );

	if( sw . O != sw . E )
	{
		free ( I );
		free ( D );
	}
	free ( T );

	int final_rot;
        if ( rrot < sl )
        {
                final_rot = rot + rrot;
        }
        else
        {
                final_rot = rot - ( 3 * sl - rrot );
        }

        if ( final_rot > ( int ) m )
        {
                ( * rotation ) = final_rot % m;
        }
        else if ( final_rot < 0 )
        {
                ( * rotation ) = m + final_rot;
        }
        else
                ( * rotation ) = final_rot;

	unsigned char * x_final_rotation = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );
	
	create_rotation( x, ( *rotation ), x_final_rotation );


	int sub = sw . S;
	int ins = sw . I; 
	int del = sw . D;

	if (  ins == 1 && del == 1 && sub == 1 )
		editDistanceMyers( x_final_rotation, y, m, n, distance );
	else
		editDistance( x_final_rotation, y, m, n, distance, sub, ins, del ); 

	free ( X );
	free ( Y );
	free( x_final_rotation );

	return EXIT_SUCCESS;
}

/*
Myers Bit-Vector algorithm implemented using SeqAn Library
www.seqan.de
*/
int editDistanceMyers( unsigned char * xInput, unsigned char * yInput, int mInput, int nInput, unsigned int * distance )
{
	typedef String<char> TSequence;

	TSequence seq1 = xInput;
	TSequence seq2 = yInput;

	int score = globalAlignmentScore( seq1, seq2, MyersBitVector() )/-1;

	( * distance ) = score;

	return EXIT_SUCCESS;
}

unsigned int editDistance(unsigned char * xInput, unsigned char * yInput, int mInput, int nInput, unsigned int * distance, int sub, int ins, int del )
{
	unsigned int x, y, lastdiag, olddiag;
	unsigned int match = 0;
	unsigned int * column;
	column = ( unsigned int * ) calloc ( mInput + 1 , sizeof(unsigned int));
	
	for (y = 1; y <= mInput; y++)
		column[y] = y;

	for (x = 1; x <= nInput; x++) 
	{
        	column[0] = x;
        	for (y = 1, lastdiag = x-1; y <= mInput; y++) 
		{
			olddiag = column[y];
            		column[y] = MIN3(column[y] + ins, column[y-1] + del, lastdiag + (xInput[y-1] == yInput[x-1] ? match : sub));
            		lastdiag = olddiag;
        	}
    	}
	
	( * distance ) = column[mInput];
	free ( column );
	return EXIT_SUCCESS;
}

