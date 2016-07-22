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

#define INITIAL_SC		-100000
#define ALLOC_SIZE               104857
#define DEL                     '$'
#define NA			'N'
#define DEL_STR                 "$"
#define GAP 			'-'
#define ERR                      24
#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*"   //Proteins alphabet
#define DNA                     "ATGCSWRYKMBVHDN"            //IUPAC alphabet
#define RNA                     "AUGCN"                      //RNA alphabet
#define NUC_SCORING_MATRIX_SIZE 15
#define PRO_SCORING_MATRIX_SIZE 24
#define WORD_LEN 		64

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)  
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))

#include <tuple>

using namespace std;
typedef unsigned long int WORD;

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   unsigned int         matrix;
   double	        P;
   int 			O, E, S, I, D;
   int 			f, g;
   unsigned int         l;
   unsigned int         q;
 };

struct TPOcc
 {
   double	        err;
   unsigned int		rot;
 };


typedef int32_t INT;

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
TPOcc * unique ( TPOcc * first, TPOcc * last );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
double gettime ( void );
void usage ( void );
int nuc_delta ( char a, char b );
int pro_delta ( char a, char b );
unsigned int nuc_char_to_index ( char a );
unsigned int pro_char_to_index ( char a );
