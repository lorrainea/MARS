/**
    MARS: Multiple circular sequence Alignment using Refined Sequences
    Copyright (C) 2016 Lorraine A.K. Ayad, Solon P. Pissis

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
#include <stdint.h>
#include <float.h>

#include <seqan/basic.h>
#include <seqan/score.h> 
#include "mars.h"

#include "EDNAFULL.h"
#include "EBLOSUM62.h"

unsigned int EDNA[90];
unsigned int BLOSUM[91];

using namespace seqan;

void init_substitution_score_tables ( void )
{
    unsigned int i;
    char edna[] = DNA;
    for ( i = 0; i < NUC_SCORING_MATRIX_SIZE; i ++ ) {
	EDNA[(int)edna[i]] = i;
    }
    EDNA[(int)'U'] = 1; //Setting RNA U = T
    char blosum[] = PROT;
    for ( i = 0; i < PRO_SCORING_MATRIX_SIZE; i ++ ) {
	BLOSUM[(int)blosum[i]] = i;
    }
}

/* Returns the score for matching character a and b based on EDNAFULL matrix */
int nuc_delta ( char a, char b )
 {
   unsigned int index_a = nuc_char_to_index ( a );
   unsigned int index_b = nuc_char_to_index ( b );

   if ( ( index_a < NUC_SCORING_MATRIX_SIZE ) && ( index_b < NUC_SCORING_MATRIX_SIZE ) )
     return ( EDNAFULL_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the score for matching character a and b based on EBLOSUM62 matrix */
int pro_delta ( char a, char b )
 {
   unsigned int index_a = pro_char_to_index( a );
   unsigned int index_b = pro_char_to_index( b );

   if ( ( index_a < PRO_SCORING_MATRIX_SIZE ) && ( index_b < PRO_SCORING_MATRIX_SIZE ) )
     return ( EBLOSUM62_matrix[ index_a ][ index_b ] );
   else //Error
     return ( ERR );
 }

/* Returns the index of char a in EDNAFULL matrix */
unsigned int nuc_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'T':
        index = 1; break;

      case 'G':
        index = 2; break;

      case 'C':
        index = 3; break;

      case 'S':
        index = 4; break;

      case 'W':
        index = 5; break;

      case 'R':
        index = 6; break;

      case 'Y':
        index = 7; break;

      case 'K':
        index = 8; break;

      case 'M':
        index = 9; break;

      case 'B':
        index = 10; break;

      case 'V':
        index = 11; break;

      case 'H':
        index = 12; break;

      case 'D':
        index = 13; break;

      case 'N':
        index = 14; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the nucleotide sequences!!!\n" );
        index = ERR; break;
    }
   
   return ( index );
 }

/* Returns the index of char a in EBLOSUM62 matrix */
unsigned int pro_char_to_index ( char a )
 {
   unsigned int index; 

   switch ( a )
    {
      case 'A':
        index = 0; break;

      case 'R':
        index = 1; break;

      case 'N':
        index = 2; break;

      case 'D':
        index = 3; break;

      case 'C':
        index = 4; break;

      case 'Q':
        index = 5; break;

      case 'E':
        index = 6; break;

      case 'G':
        index = 7; break;

      case 'H':
        index = 8; break;

      case 'I':
        index = 9; break;

      case 'L':
        index = 10; break;

      case 'K':
        index = 11; break;

      case 'M':
        index = 12; break;

      case 'F':
        index = 13; break;

      case 'P':
        index = 14; break;

      case 'S':
        index = 15; break;

      case 'T':
        index = 16; break;

      case 'W':
        index = 17; break;

      case 'Y':
        index = 18; break;

      case 'V':
        index = 19; break;

      case 'B':
        index = 20; break;

      case 'Z':
        index = 21; break;

      case 'X':
        index = 22; break;

      case '*':
        index = 23; break;

      default:
        fprintf ( stderr, "Error: unrecognizable character in one of the protein sequences!!!\n" );
        index = ERR; break;
    }
   return ( index );
 }

/*
namespace seqan 
{
	struct DNAMatrix {};

	template <>
	struct ScoringMatrixData_<int, Dna5, DNAMatrix>
	{
    		enum
    		{
        		VALUE_SIZE = ValueSize<Dna5>::VALUE,
        		TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    		};

    		static inline int const * getData()
    		{
        		static int const _data[TAB_SIZE] =
        		{
            			5, -4, -4, -4, -2,
            			-4, 5, -4, -5, -2,
            			-4, -4, 5, -4, -2,
            			-4, -4, -4, 5, -2,
            			-2, -2, -2, -2, -1
        		};
        	return _data;
    		}

	};
}*/
