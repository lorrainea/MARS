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

#include <iostream>
#include <cstdlib>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <float.h>
#include <limits.h>
#include <omp.h>
#include "mars.h"
#include "sacsc.h"
#include "ced.h"
#include "nj.h"

int main(int argc, char **argv)
{
	struct TSwitch  sw;

	FILE *          in_fd;                  // the input file descriptor
	FILE *          out_fd;                 // the output file descriptor
        char *          input_filename;         // the input file name
        char *          output_filename;        // the output file name

        unsigned char ** seq    = NULL;         // the sequence(s) in memory
        unsigned char ** seq_id = NULL;         // the sequence(s) id in memory

	char *          alphabet;               // the alphabet
	unsigned int    i, j;
	unsigned int    q, l;
	unsigned int    total_length = 0;

	/* Decodes the arguments */
        i = decode_switches ( argc, argv, &sw );

	/* Check the arguments */
        if ( i < 10 )
        {
                usage ();
                return ( 1 );
        }
        else
        {

                if      ( ! strcmp ( "DNA", sw . alphabet ) )   { alphabet = ( char * ) DNA;  sw . matrix = 0; }
                else if ( ! strcmp ( "PROT", sw . alphabet ) )  { alphabet = ( char * ) PROT; sw . matrix = 1; }
                else
                {
                        fprintf ( stderr, " Error: alphabet argument a should be `DNA' for nucleotide sequences or `PROT' for protein sequences!\n" );
                        return ( 1 );
                }

		if ( sw . P < 0 || sw . P >= 5.1 )
		{
			fprintf ( stderr, " Error: The number of refinement blocks to check at each end should be in the range 0 to 5.\n" );
			return ( 1 );
		}
		if ( sw . q < 2 )
		{
			fprintf ( stderr, " Error: The q-gram length is too small.\n" );
			return ( 1 );	
		}
		if( sw . q >= sw . l )
		{
			fprintf ( stderr, " Error: The length of the q-gram must be smaller than the block length.\n" );
			return ( 1 );
		}
	
                input_filename       = sw . input_filename;
		if ( input_filename == NULL )
		{
			fprintf ( stderr, " Error: Cannot open file for input!\n" );
			return ( 1 );
		}
		output_filename         = sw . output_filename;

        }


	double start = gettime();

	/* Read the (Multi)FASTA file in memory */
	fprintf ( stderr, " Reading the (Multi)FASTA input file: %s\n", input_filename );
	if ( ! ( in_fd = fopen ( input_filename, "r") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", input_filename );
		return ( 1 );
	}

	char c;
        unsigned int num_seqs = 0;           // the total number of sequences considered
	unsigned int max_alloc_seq_id = 0;
	unsigned int max_alloc_seq = 0;
	c = fgetc( in_fd );
	do
	{
		if ( c != '>' )
		{
			fprintf ( stderr, " Error: input file %s is not in FASTA format!\n", input_filename );
			return ( 1 );
		}
		else
		{
			if ( num_seqs >= max_alloc_seq_id )
			{
				seq_id = ( unsigned char ** ) realloc ( seq_id,   ( max_alloc_seq_id + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
				max_alloc_seq_id += ALLOC_SIZE;
			}

			unsigned int max_alloc_seq_id_len = 0;
			unsigned int seq_id_len = 0;

			seq_id[ num_seqs ] = NULL;

			while ( ( c = fgetc( in_fd ) ) != EOF && c != '\n' )
			{
				if ( seq_id_len >= max_alloc_seq_id_len )
				{
					seq_id[ num_seqs ] = ( unsigned char * ) realloc ( seq_id[ num_seqs ],   ( max_alloc_seq_id_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
					max_alloc_seq_id_len += ALLOC_SIZE;
				}
				seq_id[ num_seqs ][ seq_id_len++ ] = c;
			}
			seq_id[ num_seqs ][ seq_id_len ] = '\0';
			
		}

		if ( num_seqs >= max_alloc_seq )
		{
			seq = ( unsigned char ** ) realloc ( seq,   ( max_alloc_seq + ALLOC_SIZE ) * sizeof ( unsigned char * ) );
			max_alloc_seq += ALLOC_SIZE;
		}

		unsigned int seq_len = 0;
		unsigned int max_alloc_seq_len = 0;

		seq[ num_seqs ] = NULL;

		while ( ( c = fgetc( in_fd ) ) != EOF && c != '>' )
		{
			if( seq_len == 0 && c == '\n' )
			{
				fprintf ( stderr, " Omitting empty sequence in file %s!\n", input_filename );
				c = fgetc( in_fd );
				break;
			}
			if( c == '\n' || c == ' ' ) continue;

			c = toupper( c );

			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) );
				max_alloc_seq_len += ALLOC_SIZE;
			}

			if( strchr ( alphabet, c ) )
			{
				seq[ num_seqs ][ seq_len++ ] = c;
			}
			else
			{
				fprintf ( stderr, " Error: input file %s contains an unexpected character %c!\n", input_filename, c );
				return ( 1 );
			}

		}

		if( seq_len != 0 )
		{
			if ( seq_len >= max_alloc_seq_len )
			{
				seq[ num_seqs ] = ( unsigned char * ) realloc ( seq[ num_seqs ],   ( max_alloc_seq_len + ALLOC_SIZE ) * sizeof ( unsigned char ) ); 
				max_alloc_seq_len += ALLOC_SIZE;
			}
			seq[ num_seqs ][ seq_len ] = '\0';
			total_length += seq_len;
			num_seqs++;
		}
		
	} while( c != EOF );

	seq[ num_seqs ] = NULL;

	if ( fclose ( in_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}
	

	fprintf ( stderr, " Computing cyclic edit distance for all sequence pairs\n" );
		

	TPOcc ** D;

	int * Rot = ( int * ) calloc ( num_seqs , sizeof ( int ) );

	if ( ( D = ( TPOcc ** ) calloc ( ( num_seqs ) , sizeof( TPOcc * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: Cannot allocate memory!\n" );
		exit( EXIT_FAILURE );
	}

	for ( int i = 0; i < num_seqs; i++ )
	{
		unsigned int m = strlen ( ( char * ) seq[i] );

		if( sw . P * sw . l * 3 > m )
		{
			fprintf( stderr, " Error: P is too large!\n" );
			exit( EXIT_FAILURE );
		}
		
		if ( ( D[i] = ( TPOcc * ) calloc ( ( num_seqs + 1 ) , sizeof( TPOcc ) ) ) == NULL )
		{
			fprintf( stderr, " Error: Cannot allocate memory!\n" );
			exit( EXIT_FAILURE );
		}
			
		for ( int j = 0; j < num_seqs; j ++ )
		{
			if ( i == j ) 
				continue;

			unsigned int n = strlen ( ( char * ) seq[j] );	
			unsigned char * xr = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );

			if ( sw . l > m - sw . q + 1  || sw . l > n - sw . q + 1 )
			{
				fprintf( stderr, " Error: Illegal block length.\n" );
				exit ( 1 );
			}

			unsigned int distance = ( int ) DBL_MAX;
			unsigned int rotation = 0;

			/*Finds an approximate rotation for every pair of sequences in the data sets*/ 
			circular_sequence_comparison ( seq[i], seq[j], sw, &rotation, &distance );

			create_rotation ( seq[i], rotation, xr );

			/*Produces more accurate rotations using refined sequences*/
			sacsc_refinement(seq[i], xr, seq[j], sw, &rotation, &distance);

			D[i][j] . err = distance;
			D[i][j] . rot = rotation;

			free( xr );

		}
	}

	fprintf ( stderr, " Creating the guide tree\n" );

	/*Creates the guide tree*/
	nj ( D, num_seqs, seq, sw, Rot);

	fprintf ( stderr, " Preparing the output\n" );

	if ( ! ( out_fd = fopen ( output_filename, "w") ) )
	{
		fprintf ( stderr, " Error: Cannot open file %s!\n", output_filename );
		return ( 1 );
	}

	unsigned int succ = 0;
	unsigned int unsucc = 0;
	for ( int i = 0; i < num_seqs; i ++ )
	{
		if ( Rot[i] >= 0 )
		{
			Rot[i] = Rot[i] % strlen ( ( char * ) seq[i] );
			unsigned char * rotation;
			unsigned int m = strlen ( ( char * ) seq[i] );
			rotation = ( unsigned char * ) calloc ( m + 1, sizeof ( unsigned char ) );
			create_rotation ( seq[i], Rot[i], rotation );
			fprintf( out_fd, ">%s\n", seq_id[i] );
			fprintf( out_fd, "%s\n", rotation );
			free ( rotation );
			succ++;
		}		
	}
		
	if ( fclose ( out_fd ) )
	{
		fprintf( stderr, " Error: file close error!\n");
		return ( 1 );
	}

	free ( Rot );
	for ( i = 0; i < num_seqs; i ++ )
	{
		free ( D[i] );
	}
	free ( D );

	double end = gettime();

        fprintf( stderr, "Elapsed time for processing %d sequence(s): %lf secs.\n", num_seqs, ( end - start ) );
	
	/* Deallocate */
	for ( i = 0; i < num_seqs; i ++ )
	{
		free ( seq[i] );
		free( seq_id[i] );
	}	
	free ( seq );
	free ( seq_id );
        free ( sw . input_filename );
        free ( sw . output_filename );
        free ( sw . alphabet );

	return ( 0 );
}
