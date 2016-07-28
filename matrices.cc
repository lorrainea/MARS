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
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "mars.h"

unsigned int EDNA[90];
unsigned int BLOSUM[91];

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
