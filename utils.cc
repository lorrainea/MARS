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
#include <sys/time.h>
#include <getopt.h>
#include <assert.h>
#include "mars.h"

static struct option long_options[] =
 {
   { "alphabet",                required_argument, NULL, 'a' },
   { "seqs-file",               required_argument, NULL, 'i' },
   { "output-file",             required_argument, NULL, 'o' },
   { "block-length",            required_argument, NULL, 'l' },
   { "q-length",                required_argument, NULL, 'q' },
   { "gap-open",                optional_argument, NULL, 'O' },
   { "gap-extend",              optional_argument, NULL, 'E' },
   { "refine-blocks",           optional_argument, NULL, 'P' },
   { "cost-substitution",	optional_argument, NULL, 'S' },
   { "cost-insertion",		optional_argument, NULL, 'I' },
   { "cost-deletion",		optional_argument, NULL, 'D' },
   { "score-insertion",     	optional_argument, NULL, 'f' },
   { "score-deletion",     	optional_argument, NULL, 'g' },
   { "help",                    no_argument,       NULL, 'h' },
   { NULL,                      0,                 NULL,  0  }
 };


/* 
Decode the input switches 
*/
int decode_switches ( int argc, char * argv [], struct TSwitch * sw )
 {
   int          oi;
   int          opt;
   double       val;
   char       * ep;
   int          args;

   /* initialisation */
   sw -> alphabet                       = NULL;
   sw -> input_filename                 = NULL;
   sw -> output_filename                = NULL;
   sw -> O                              = -10;
   sw -> E                              = -2;
   sw -> P                              = 0.0;
   sw -> S				= 1;
   sw -> I 				= 1;
   sw -> D				= 1;
   sw -> f				= -4;
   sw -> g				= -4;
   sw -> l                              = 50;
   sw -> q                              = 5;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:l:q:S:I:D:f:g:O:E:T:P:h", long_options, &oi ) ) != -1 ) 
    {

      switch ( opt )
       {
         case 'a':
           sw -> alphabet = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> alphabet, optarg );
           args ++;
           break;

         case 'i':
           sw -> input_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> input_filename, optarg );
           args ++;
           break;

         case 'o':
           sw -> output_filename = ( char * ) malloc ( ( strlen ( optarg ) + 1 ) * sizeof ( char ) );
           strcpy ( sw -> output_filename, optarg );
           args ++;
          break;

          case 'l':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> l = val;
	   args++;
           break;

          case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
	   args++;
           break;

	case 'S':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> S = val;
           break;

	case 'I':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> I = val;
           break;

	case 'D':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> D = val;
           break;

	case 'f':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> f = val;
           break;

	case 'g':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> g = val;
           break;

         case 'O':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> O = val;
           break;

         case 'P':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> P = val;
           break;

         case 'E':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> E = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 4)
     {
       usage ();
       exit ( 1 );
     }
   else
     return ( optind );
 }

/* 
Usage of the tool 
*/
void usage ( void )
 {
   fprintf ( stdout, " Usage: bear <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet              <str>     `DNA' for nucleotide  sequences  or `PROT' for\n"
                     "                                        protein  sequences. \n" );
   fprintf ( stdout, "  -i, --input-file            <str>     MultiFASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file           <str>     Output filename with rotated sequences.\n" );
   fprintf ( stdout, "  -q, --q-length              <int>     The q-gram length.\n");
   fprintf ( stdout, "  -l, --block-length          <int>     The length of each block.\n\n");
   fprintf ( stdout, " Optional:\n" ); 
   fprintf ( stdout, "  -f, --score-insertion       <int>     Score of inserting a character in alignment. Default: -4.\n" );
   fprintf ( stdout, "  -g, --score-deletion        <int>     Score of deleting a character in alignment. Default: -4.\n" );
   fprintf ( stdout, "  -O, --gap-open              <int>     Affine gap open penalty. Default: -10.\n" );
   fprintf ( stdout, "  -E, --gap-extend            <int>     Affine gap extension penalty. Default: -2.\n" );
   fprintf ( stdout, "  -P, --refine-blocks         <dbl>     Refine the alignments by checking a block\n"
		     "                                        percentage of the ends.\n" );
   fprintf ( stdout, "  -S, --cost-substitution     <int>     Cost of substitution for edit distance. Default: 1.\n" );
   fprintf ( stdout, "  -I, --cost-insertion        <int>     Cost of insertion for edit distance. Default: 1.\n" );
   fprintf ( stdout, "  -D, --cost-deletion         <int>     Cost of deletion for edit distance. Default: 1.\n");

 }

double gettime( void )
{
    struct timeval ttime;
    gettimeofday( &ttime , 0 );
    return ttime.tv_sec + ttime.tv_usec * 0.000001;
}

void create_rotation ( unsigned char * x, unsigned int offset, unsigned char * rotation )
{
	unsigned int m = strlen ( ( char * ) x );
	memmove ( &rotation[0], &x[offset], m - offset );
	memmove ( &rotation[m - offset], &x[0], offset );
	rotation[m] = '\0';
}

