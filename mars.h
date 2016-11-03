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
#define DL                      '$'
#define NA			'N'
#define DEL_STR                 "$"
#define GAP 			'-'
#define ERR                      24
#define PROT                    "ARNDCQEGHILKMFPSTWYVBZX*"   //Proteins alphabet
#define DNA                     "ATGCSWRYKMBVHDN"            //IUPAC alphabet
#define INS			1
#define DEL			1
#define SUB			1
#define MAT			0

#define NUC_SCORING_MATRIX_SIZE 15
#define PRO_SCORING_MATRIX_SIZE 24

#define MAX2(a,b) ((a) > (b)) ? (a) : (b)
#define MIN2(a,b) ((a) < (b)) ? (a) : (b)  
#define MAX3(a, b, c) ((a) > (b) ? ((a) > (c) ? (a) : (c)) : ((b) > (c) ? (b) : (c)))
#define nuc_delta(a,b) EDNAFULL_matrix[ EDNA[(int)(a)] ][ EDNA[(int)(b)] ]
#define pro_delta(a,b) EBLOSUM62_matrix[ BLOSUM[(int)(a)] ][ BLOSUM[(int)(b)] ]

using namespace std;
typedef unsigned long int WORD;

struct TSwitch
 {
   char               * alphabet;
   char               * input_filename;
   char               * output_filename;
   unsigned int         matrix;
   double	        P;
   int 			O, E, U, V, S, I, D, T;
   unsigned int         l, q, m;
 };

struct TPOcc
 {
   double	        err;
   unsigned int		rot;
 };


typedef int32_t INT;

extern unsigned int EDNA[];
extern unsigned int BLOSUM[];

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation );
TPOcc * unique ( TPOcc * first, TPOcc * last );
int decode_switches ( int argc, char * argv [], struct TSwitch * sw );
double gettime ( void );
void usage ( void );
void init_substitution_score_tables ( void );

void cyclic(unsigned char *pattern1, unsigned char *pattern2, int length1, int length2, int time_rep, int norm, unsigned int * rotation, unsigned int * distance);
extern void centra(char *string0, int a, char *sep, char *dummy ) ;
