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
   { "gap-open-seq",            optional_argument, NULL, 'O' },
   { "gap-extend-seq",          optional_argument, NULL, 'E' },
   { "gap-open-pro",            optional_argument, NULL, 'U' },
   { "gap-extend-pro",          optional_argument, NULL, 'V' },
   { "refine-blocks",           required_argument, NULL, 'P' },
   { "method",			required_argument, NULL, 'm' },
   { "threads", 		optional_argument, NULL, 'T' },
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
   sw -> E                              = -1;
   sw -> U                              = -10;
   sw -> V                              = -1;
   sw -> P                              = 1.0;
   sw -> l                              = 25;
   sw -> q                              = 5;
   sw -> m				= 0;
   sw -> T                              = 1;
   args = 0;

   while ( ( opt = getopt_long ( argc, argv, "a:i:o:l:q:m:U:V:O:E:T:P:h", long_options, &oi ) ) != -1 ) 
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
           break;

          case 'q':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> q = val;
           break;

	case 'V':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> V = val;
           break;

	case 'U':
           val = strtol ( optarg, &ep, 10 );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> U = val;
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

	 case 'T':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> T = val;
           break;

	 case 'm':
           val = strtod ( optarg, &ep );
           if ( optarg == ep )
            {
              return ( 0 );
            }
           sw -> m = val;
           break;

         case 'h':
           return ( 0 );
       }
    }

   if ( args < 3 )
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
   fprintf ( stdout, " Usage: mars <options>\n" );
   fprintf ( stdout, " Standard (Mandatory):\n" );
   fprintf ( stdout, "  -a, --alphabet              <str>     'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.\n" );
   fprintf ( stdout, "  -i, --input-file            <str>     MultiFASTA input filename.\n" );
   fprintf ( stdout, "  -o, --output-file           <str>     Output filename with rotated sequences.\n" );    
   fprintf ( stdout, " Optional:\n" );
   fprintf ( stdout, " Cyclic edit distance computation between pairs of sequences.\n" );   
   fprintf ( stdout, "  -m, --method                <int>     0 for heuristic cyclic edit distance hCED - Faster but less accurate. \n"
   		     "                                        1 for branch and bound method - Slower but more accurate. Default: 0.\n" );
   fprintf ( stdout, " Computation of pairwise cyclic edit distance for method hCED. \n" );
   fprintf ( stdout, "  -q, --q-length              <int>     The q-gram length. Default: 5.\n" );
   fprintf ( stdout, "  -l, --block-length          <int>     The length of each block. Default: 25.\n" );   
   fprintf ( stdout, " Refining pairwise rotations.\n" );
   fprintf ( stdout, "  -P, --refine-blocks         <dbl>     Refine the rotations by aligning P blocks of the ends. Default: 1.\n" );
   fprintf ( stdout, "  -O, --gap-open-seq          <int>     Affine gap open penalty in pairwise block alignment. Default: -10.\n" );
   fprintf ( stdout, "  -E, --gap-extend-seq        <int>     Affine gap extension penalty in pairwise block alignment. Default: -1.\n" ); 
   fprintf ( stdout, " Progressive alignment of profiles.\n" ); 
   fprintf ( stdout, "  -U, --gap-open-pro          <int>     Affine gap open penalty in alignment of profiles. Default: -10.\n" );
   fprintf ( stdout, "  -V, --gap-extend-pro        <int>     Affine gap extension penalty in alignment of profiles. Default: -1.\n" );
   fprintf ( stdout, " Number of threads.\n" ); 
   fprintf ( stdout, "  -T, --threads               <int>     Number of threads to use. Default: 1. \n" );
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

