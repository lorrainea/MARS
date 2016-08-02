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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <time.h>
#include <seqan/graph_msa.h>
#include <seqan/align.h>
#include <vector>
#include "EDNAFULL.h"
#include "EBLOSUM62.h"
#include "mars.h"
#include "sacsc.h"
#include "nj.h"

using namespace seqan;
using namespace std;

unsigned int progAlignment(TPOcc ** D, unsigned char ** seq, TGraph njTree, struct TSwitch  sw, int * Rot, vector<array<int, 2>> * branchingOrder, unsigned int num_seqs )
{
	
	init_substitution_score_tables ();
	int * R = ( int * ) calloc ( num_seqs , sizeof ( int ) );

	unsigned char ** sequences;
	if ( ( sequences = ( unsigned char  ** ) calloc ( ( num_seqs ) , sizeof( unsigned char * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: Sequences could not be allocated!\n");
                return ( 0 );
        }
        for ( int i = 0; i < num_seqs; i ++ )
        {
                if ( ( sequences[i] = ( unsigned char * ) calloc ( ( ALLOC_SIZE ) , sizeof( unsigned char ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: Sequences could not be allocated!\n");
                        return ( 0 );
                }
        }

	for(int i=0; i<num_seqs; i++)
	{
		memcpy( sequences[i], seq[i],  strlen ( (  char * ) seq[i] )*sizeof( unsigned char )  ) ;
	}

	
	for( int i=0; i<branchingOrder->size(); i++ ) //traverse through all nodes in tree
	{
		array<int , 2> children;
		children = branchingOrder->at(i); // child at each node

		if( isLeaf( njTree, children[0] ) == 1 && isLeaf( njTree, children[1] ) == 1 ) // two sequences
		{
			int m = strlen( ( char * ) seq[ children[0] ] );
	
			R[ children[0] ] = D[ children[0] ][ children[1] ] . rot; // obtain first rotation for two sequences (already refined)

			Rot[ children[0] ] =  D[ children[0] ][ children[1] ] . rot; // Accumulated rotation array

			unsigned char * rotatedSeq = ( unsigned char * ) calloc( ( m + 1 ) , sizeof( unsigned char ) );

			create_rotation( seq[ children[0] ] , R[ children[0] ], rotatedSeq );

			alignPairs ( rotatedSeq, seq[ children[1] ], sequences, sw, children[0], children[1] ); // find new alignments

			free( rotatedSeq );

		}
		else 
		{

			vector<unsigned char * > * profileA = new vector<unsigned char *>();  // holds all alignments in profile A
			vector<unsigned char * > * profileB = new vector<unsigned char *>();  // holds all alignments in profile B
			vector<char> * characters = new vector<char>(); // hold characters in profile A

			vector<int> * profileAPos = new vector<int>(); // holds all leaves of profile A
			vector<int> * profileBPos = new vector<int>(); // holds all leaves of profile B

			String<int> leaves0;
			collectLeaves(njTree, children[0] , leaves0); // find all leaves of child 0

	    		typedef Iterator<String<int>>::Type TIterator;
	    		for (TIterator it = begin(leaves0); it != end(leaves0); goNext(it))
	    		{
	       			profileAPos->push_back( value(it) );
	   		}

			String<int> leaves1;
			collectLeaves(njTree, children[1] , leaves1); // find all leaves of child 1

	    		typedef Iterator<String<int>>::Type TIterator;
	    		for (TIterator it2 = begin(leaves1); it2 != end(leaves1); goNext(it2))
	    		{
	       			profileBPos->push_back( value(it2) );
	   		}			
	   			
	   		if( profileAPos->size() < profileBPos->size() )
			{
				(*profileBPos).swap(*profileAPos); //Profile A will always have largest number of leaves
			}

			int m = strlen( ( char * ) sequences[ profileAPos->at(0) ] );
			int n = strlen( ( char * ) sequences[ profileBPos->at(0) ] );

			int rotValue = 0;
			int bValue = 0;
			int aValue = 0;

			/* Identify the most suitable initial approximate rotation using array R and initial distance matrix D */
	   		for(int j=0; j<profileAPos->size(); j++)
	   		{
	   			for(int i =0; i<profileBPos->size(); i++)
	   			{

					if( R[ profileBPos->at(i) ] == 0 && R[ profileAPos->at(j) ] == 0 )
					{
						rotValue = D[ profileBPos->at(i) ] [ profileAPos->at(j) ] . rot;
						bValue = i;
						aValue = j;
					}
				}

	   		}

			int rot = rotValue;
			int rs = sw . l * sw . P;

			int gapCountA = 0;
			int gapCountB = 0;
			int largestNoGapsA = 0;
			int largestNoGapsB = 0;

			if( n == strlen( ( char * ) seq[ bValue ] ) )
			{
				rot = rot;
			}
			else
			{
				for(int i=0; i<=rotValue; i++ )
				{
					if( sequences[ profileBPos->at(bValue) ][i] == GAP )
					{
						largestNoGapsB++;
						rotValue ++;
					}
				}
					
				rot =  rot + largestNoGapsB;
			}

			unsigned char ** initial_rotation = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );
			for(int i=0; i<profileBPos->size(); i++ )
			{
				initial_rotation[i] = ( unsigned char * ) calloc( ( n + 1 ) , sizeof( unsigned char ) );
			}
			
			for(int j=0; j<profileBPos->size(); j++)
			{
				create_rotation( sequences[ profileBPos->at( j ) ] , rot, initial_rotation[j] );
			}	
		   
			unsigned char ** profB = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );
			for(int i=0; i<profileBPos->size(); i++)
				profB[i] = ( unsigned char * ) calloc( ( 3 * rs + 1 ) , sizeof( unsigned char ) );
	   		
	   		for(int j=0; j<profileBPos->size(); j++)
	   		{
				memcpy ( &profB[j][0], &initial_rotation[j][0], rs );
				for ( int i = 0; i < rs; i++ )
					profB[j][rs + i] = DEL;
					
				memcpy ( &profB[j][rs + rs], &initial_rotation[j][n - rs], rs );
				profB[j][ 3* rs] = '\0';
	
				profileB->push_back( profB[j] );  
			}

			unsigned char ** profA = ( unsigned char ** ) calloc( ( profileAPos->size() ) , sizeof( unsigned char * ) );
			for(int j=0; j<profileAPos->size(); j++)
				profA[j] = ( unsigned char * ) calloc( ( 3*rs + 1 ) , sizeof( unsigned char ) );	
		
			for(int i=0; i<profileAPos->size(); i++)
			{	
				memcpy ( &profA[i][0], &sequences[ profileAPos->at(i)][0], rs );
				for ( int k = 0; k < rs; k++ )
					profA[i][rs + k] = DEL;
						
				memcpy ( &profA[i][rs + rs], &sequences[ profileAPos->at(i)][m - rs], rs );
				profA[i][3*rs] = '\0';
			
				profileA->push_back( profA[i] );

			}


			// begin progressive alignment for every rotation of sequences in profileB
			double score = INITIAL_SC;
			int rotation = 0;
	
			if( rs > 0 )
			{	
				int ** TB; 
				double ** SM, ** PM, * IM, * DM;
			
				if( sw . U != sw . V )
					alignAllocation_ag( PM, SM, IM, DM, TB, characters, profileA, profileB, sw );	
				else alignAllocation( PM, SM, TB, characters, profileA, profileB, sw );
			
				profileB->clear(); // clear profileB so can be re-inserted with refined sequences

				/* Find the best rotation from the refined sequence using DP score */
				unsigned char ** rotatedSeq = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );

				for(int i=0; i<3*rs; i++ )
				{
					if ( i >=rs && i < 2 * rs )
						continue;

					for(int j = 0; j< profileBPos->size(); j++)
						rotatedSeq[j] = ( unsigned char * ) calloc( ( 3 * rs + 1 ) , sizeof( unsigned char  ) );
					

					for(int j=0; j<profileBPos->size(); j++)
					{
						create_rotation( profB[j], i , rotatedSeq[j] );
						profileB->push_back( rotatedSeq[j] );
					}

					if( sw . U != sw . V )
						alignmentScore_ag( profileA, profileB, &score, sw, i, &rotation, sequences, TB, SM, PM, IM, DM, characters, 0);
					else alignmentScore( profileA, profileB, &score, sw, i, &rotation, sequences, TB, SM, PM, characters, 0);
					
					for(int k = 0; k < profileBPos->size(); k++ )
						free( rotatedSeq[k] );

					profileB->clear();
				}
				free( rotatedSeq );			
			

				for(int i=0; i<characters->size(); i++)
					free( PM[i] );
		
				for(int j=0; j<3 * rs + 1; j++)
				{
					free( SM[j] );
					free( TB[j] );
				}
				free( TB );
				if( sw . U != sw . V )
				{
					free( IM );
					free( DM );
				}
				free( SM );
				free( PM );
			}
			else profileB->clear();

			for(int i=0; i<profileBPos->size(); i++)
			{
				free( profB[i] );
				free( initial_rotation[i] );
			}	

			for(int i=0; i<profileAPos->size(); i++)
			{
				free( profA[i] );
			}	
			
			free( profA );
			free( profB );
			free( initial_rotation );

			characters->clear();
			
			//calculate the final rotation using initial rotation and new rotation calculated
			int final_rot = 0;
			if ( rotation <= rs )
			{
				final_rot = rot + rotation;
			}
			else
			{
		 		final_rot = rot - ( 3*rs - rotation );
			}

			if ( final_rot > ( int ) n )
			{
		 		final_rot = final_rot % n;
			}
			else if ( final_rot < 0 )
			{
				final_rot = final_rot + n;
			}
			else final_rot = final_rot;
				
			/* Add the rotation value to the rotation arrays */
			for(int i =0; i<profileBPos->size(); i ++)
			{
				Rot[ profileBPos->at(i) ] = Rot[ profileBPos->at(i) ] + final_rot;				
				for(int j=0; j<=final_rot; j++)
				{
					if( sequences[ profileBPos->at(i)][j] == '-' )
					{
						Rot[ profileBPos->at(i) ] = Rot[ profileBPos->at(i) ] - 1;
					}
				}
				R[ profileBPos->at(i) ] = ( R[ profileBPos->at(i) ] + final_rot ) % n;
			}

			/* Create final rotation for sequences in profile B and place back into vector B */
			unsigned char ** final_rotation = ( unsigned char ** ) calloc( ( profileBPos->size() ) , sizeof( unsigned char * ) );

			for(int i =0; i<profileBPos->size(); i++)
				final_rotation[i] = ( unsigned char * ) calloc( ( n + 1 ) , sizeof( unsigned char  ) );

			for(int j=0; j<profileBPos->size(); j++)
			{
				create_rotation( sequences[ profileBPos->at(j ) ], final_rot, final_rotation[j] );
				profileB->push_back( final_rotation[j] );
			}
			

			profileA->clear();
			
			/* Place original sequences in profile A back into vector A */
			for(int i=0; i<profileAPos->size(); i++)
			{
				profileA->push_back( sequences[ profileAPos->at(i) ] );			
			}
	
			int ** TBl;
			double ** SMl, ** PMl, * IMl, * DMl;

			if( sw . U != sw . V )
				alignAllocation_ag( PMl, SMl, IMl, DMl, TBl, characters, profileA, profileB, sw );	
			else alignAllocation( PMl, SMl, TBl, characters, profileA, profileB, sw );

			/* Calculate traceback matrix */
			if( sw . U != sw . V )
				alignmentScore_ag( profileA, profileB, &score, sw, 0, &rotation, sequences, TBl, SMl, PMl, IMl, DMl, characters, 1);
			else alignmentScore( profileA, profileB, &score, sw, 0, &rotation, sequences, TBl, SMl, PMl, characters, 1);

			for(int i=0; i<characters->size(); i++)
				free( PMl[i] );
			delete( characters );

			for(int i=0; i<m+ 1; i++)
				free( SMl[i] );

			if( sw . U != sw . V )
			{
				free( IMl );
				free( DMl );
			}
			free( PMl );
			free( SMl );
			
			alignSequences( profileA, profileB, profileAPos, profileBPos, sequences , TBl); // find the best alignment for sequences

			//free arrays 
			for(int i=0; i<profileBPos->size(); i++)
				free( final_rotation[i] );
	
			free( final_rotation );
			delete( profileA );
			delete( profileB );
			delete( profileAPos );
			delete( profileBPos );
			for(int i=0; i<m+ 1; i++)
				free( TBl[i] );
			
			free( TBl );
			
			
			
		}
	}

	for ( int i = 0; i < num_seqs; i ++ )
	{
		free ( sequences[i] );
	}	
	free ( sequences );
	free( R );


return 0;
}

unsigned int alignSequences(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, vector<int> * profileAPos, vector<int> * profileBPos, unsigned char ** sequences , int ** &TB)
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	int seqP = 0; // sizes of new sequences for profile A and profile B
	int seqS = 0; // m and n will increase when gaps are added
	
	unsigned char ** ASequences = ( unsigned char ** ) calloc( ( profileA->size() ) , sizeof( unsigned char * ) );
	for(int i=0; i< profileA->size(); i++)
		 ASequences[i] = ( unsigned char * ) calloc ( ( ALLOC_SIZE ) , sizeof( unsigned char ));
		
	unsigned char ** BSequences = ( unsigned char ** ) calloc( ( profileB->size() ) , sizeof( unsigned char * ) );
	for(int j=0; j< profileB->size(); j++)
		 BSequences[j] = ( unsigned char * ) calloc ( ( ALLOC_SIZE ) , sizeof( unsigned char ));

	int dirI = m; // direction traceback takes for i
	int dirJ = n; // direction traceback takes for j
	
	int i = m; // position of element in profile A
	int j = n; // position of element in profile B

	while ( dirI != 0 || dirJ != 0 )
	{
		if( TB[dirI][dirJ] == 0 )
		{
			for(int k=0; k<profileAPos->size(); k++)
			{
				ASequences[ k ][ seqP]  = profileA->at(k )[ i-1 ];
			}
			for(int l=0; l<profileBPos->size(); l++ )
			{
				BSequences[ l ][ seqS]  = profileB->at( l )[ j-1 ];
			}
			
			seqP++; seqS++;	
			i--;  j--;
			dirI = dirI-1; dirJ=dirJ-1;
		}
		else if( TB[dirI][dirJ] == 1 ) // place gap in sequence
		{
			for(int k=0; k<profileAPos->size(); k++)
			{
				ASequences[ k ][ seqP]  = profileA->at(k )[ i-1 ];
			}
			
			for(int l =0; l<profileB->size(); l++)
			{			
				BSequences[ l ][ seqS ] = '-';
			}

			seqP++; seqS++;
			i--;
			dirI = dirI-1; dirJ = dirJ;
		}
		else if( TB[dirI][dirJ]  == -1 ) // place gap in profile
		{
			for(int k =0; k<profileAPos->size(); k++)
			{
				ASequences[ k ][ seqP ] = '-'; 
			}
			
			for(int l=0; l<profileBPos->size(); l++ )
			{
				BSequences[ l ][ seqS]  = profileB->at( l )[ j-1 ];
			}

			seqP++; seqS++;
			j--;
			dirI = dirI; dirJ=dirJ-1;
			
		}

	}
	
	for(int a=0; a<profileA->size(); a++)
	{
		int k = seqP-1;
		for(int i=0; i<seqP; i++)
		{
			sequences[profileAPos->at(a)][i] = ASequences[a][k];
			k --;
			 			
		}	

		sequences[ profileAPos->at(a) ][ seqP ] = '\0';
	}
		
	for(int b=0; b<profileB->size(); b++)
	{
		int l = seqS-1;
		for(int j=0; j<seqS; j++)
		{
			sequences[ profileBPos->at(b) ][j] = BSequences[b][l];
			l--;
		}
	
		sequences[ profileBPos->at(b)][ seqS] = '\0'; 
	}

	
	for(int i=0; i<profileA->size(); i++)
		free( ASequences[i] );
			
	for(int j=0; j<profileB->size(); j++)
		free( BSequences[j] );
			
	free( ASequences );
	free( BSequences );

return 0;
}


unsigned int alignAllocation( double ** &PM, double ** &SM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB, struct TSwitch sw)
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	// add all characters in profileA into vector
	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			if(find(characters->begin(), characters->end(), profileA->at(i)[j]) != characters->end())  
			{
	    			continue;
			}	
			else characters->push_back( profileA->at(i)[j] );
		}
	}

	
	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
		  	return ( 0 );
		}
	}
			 
	//probability matrix  
	if ( ( PM = ( double ** ) calloc ( ( characters->size() ) , sizeof( double * ) ) ) == NULL )
   	{
               fprintf( stderr, " Error: PM could not be allocated!\n");
                return ( 0 );
        }

        for ( int i = 0; i < characters->size(); i ++ )
        {
		
                if ( ( PM[i] = ( double * ) calloc ( ( m ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: PM could not be allocated!\n");
                        return ( 0 );
                }
        }

   	//scoring matrix
	if ( ( SM = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: SM could not be allocated!\n");
                return ( 0 );
        }

        for ( int i = 0; i < m + 1; i ++ )
        {
                if ( ( SM[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: SM could not be allocated!\n");
                        return ( 0 );
                }
        }


	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
         	
         for(int j=1; j<n+1; j++)
		TB[0][j] = -1;

	SM[0][0] = 0;
	for ( int i = 1; i < m +1 ; i++ )
	SM[i][0] = i * sw . g;

	for ( int j = 1; j < n + 1 ; j++ )
	SM[0][j] = j * sw . f;

	double prob = 1.0/profileA->size();

	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			int pos = find(characters->begin(), characters->end(), profileA->at(i)[j]) - characters->begin() ;
			PM[ pos ][ j ] = PM[ pos ][ j ] + prob;
		}
	}
}

unsigned int alignmentScore(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, unsigned char ** sequences, int ** &TB, double ** &SM, double ** &PM, vector<char> * characters, unsigned int calculate_TB)
{

	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );


	double u;
	double v;
	double w;
	for(int i=1; i<m+1; i++) 
	{
		for(int j=1; j<n+1; j++)
		{
			u = SM[i-1][j-1] + probScore( characters, i, j, PM, profileA, profileB, sw);
			v = SM[i-1][j] + sw . g; // gap in sequence
			w = SM[i][j-1] + sw . f; // gap in profile

			SM[i][j] = MAX3 ( u, v, w );
			
			if( calculate_TB == 1 )
			{
			
				if( SM[i][j] == u)
				{
					TB[i][j] = 0;
				}
				else if(SM[i][j] == w )
				{
					TB[i][j] = -1; 
				}
				else if( SM[i][j] == v)
				{
					TB[i][j] = 1;
				}
			}

		}
	}

	if( SM[m][n] > (*score) )
	{	
		( * rotation ) = i;	
		( *score ) = SM[m][n];
	}

return 0;
}

double probScore( vector<char> * characters, int i, int j, double ** PM, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, struct TSwitch sw )
{
	double score = 0;
	for(int k=0; k<profileB->size(); k++)
	{
		for( int l=0; l<characters->size(); l++)
		{
			score = score + ( similarity(  profileB->at(k)[j-1], characters->at(l) , sw ) * PM[l][i-1] ) ; 
		}
	}

return score;
}

int similarity( unsigned char x, unsigned char y, struct TSwitch sw)
{
	int sim = 0;
	
	if ( x == DEL || y == DEL ) 
   	{
		sim = 0;
   	}
	else if ( x == GAP && y == GAP)
	{
		int matching_score = ( sw . matrix ? pro_delta( NA, NA ) : nuc_delta( NA, NA ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
	else if( x == GAP )
	{
		int matching_score = ( sw . matrix ? pro_delta( NA, y ) : nuc_delta( NA, y ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
	else if( y == GAP )
	{
		int matching_score = ( sw . matrix ? pro_delta( x, NA ) : nuc_delta( x, NA ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
	}
    	else
    	{
		int matching_score = ( sw . matrix ? pro_delta( x, y ) : nuc_delta( x, y ) ) ;
		if ( matching_score == ERR )
			sim = 0;
		else sim = matching_score;
    	}	

return sim;
}

unsigned int alignPairs(unsigned char * x, unsigned char * y, unsigned char ** sequences, struct TSwitch sw, int posX, int posY)
{

	typedef String<char> TSequence;
	typedef Align<TSequence, ArrayGaps> TAlign;
	typedef Row<TAlign>::Type TRow;
	typedef Iterator<TRow>::Type TRowIterator;

	TAlign align;
    	resize(rows(align), 2);
    	assignSource(row(align, 0), x);
    	assignSource(row(align, 1), y);


	typedef int TValue;
    	typedef Score<TValue, ScoreMatrix<Dna5, Default> > DNAMatrix;


	if( sw . O != sw . E )
		int score = globalAlignment(align, Score<int, Simple>(0, sw.g, sw. E, sw.O), AffineGaps() );
	else
		int score = globalAlignment(align, Score<int, Simple>(0, sw . f, sw . g));

	unsigned xLength = length(row(align, 0)); //sequence x in alignment
	unsigned yLength = length(row(align, 1)); //sequence y in alignment

	//iterate through alignment to find new sequence x
	int pos = 0;
      	TRowIterator it = iter(row(align, 0), 0);
        TRowIterator itEnd = iter(row(align, 0), xLength);
       	while (it != itEnd)
        {
           	sequences[posX][pos] = ( ( unsigned char ) value(it) );
        	++it;
		++pos;
    	}
	sequences[posX][pos] = '\0'; // add new x to sequences

	//iterate through alignment to find new sequence y
	int pos2 = 0;
      	TRowIterator it2 = iter(row(align, 1), 0);
        TRowIterator itEnd2 = iter(row(align, 1), yLength);
       	while (it2 != itEnd2)
        {
         	sequences[posY][pos2] = ( ( unsigned char ) value(it2) );
           	++it2;
		++pos2;
        }
	sequences[posY][pos2] = '\0'; // add new y to sequences


return 0;
}

unsigned int alignAllocation_ag( double ** &PM, double ** &SM, double * &IM, double * &DM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB, struct TSwitch sw )
{
	
	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	// add all characters in profileA into vector
	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			if(find(characters->begin(), characters->end(), profileA->at(i)[j]) != characters->end())  
			{
	    			continue;
			}	
			else characters->push_back( profileA->at(i)[j] );
		}
	}

	if ( ( TB = ( int ** ) calloc ( ( m + 1 ) , sizeof( int * ) ) ) == NULL )
	{
		fprintf( stderr, " Error: TB could not be allocated!\n");
		return ( 0 );
	}
	for ( int i = 0; i < m + 1; i ++ )
	{
		if ( ( TB[i] = ( int * ) calloc ( ( n + 1 ) , sizeof( int ) ) ) == NULL )
		{
			fprintf( stderr, " Error: TB could not be allocated!\n");
			return ( 0 );
		}
	}
			 
	//probability matrix
	if ( ( PM = ( double ** ) calloc ( ( characters->size() ) , sizeof( double * ) ) ) == NULL )
   	{
               fprintf( stderr, " Error: PM could not be allocated!\n");
               return ( 0 );
        }


        for ( int i = 0; i < characters->size(); i ++ )
        {
		
                if ( ( PM[i] = ( double * ) calloc ( ( m ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: PM could not be allocated!\n");
                        return ( 0 );
                }
        }
  
     
	if ( ( SM = ( double ** ) calloc ( ( m + 1 ) , sizeof( double * ) ) ) == NULL )
        {
                fprintf( stderr, " Error: SM could not be allocated!\n");
                return ( 0 );
        }

        for ( int i = 0; i < m + 1; i ++ )
        {
                if ( ( SM[i] = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) ) ) == NULL )
                {
                        fprintf( stderr, " Error: SM could not be allocated!\n");
                        return ( 0 );
                }
        }
        
	IM = ( double * ) calloc ( ( n + 1 ) , sizeof( double ) );
	DM = ( double * ) calloc ( ( m + 1 ) , sizeof( double ) );

	for ( int i = 0; i < m + 1; i++ )
	{
		DM[i] =  m * sw . g;
	}
	for ( int j = 0; j < n + 1; j++ )
	{
		IM[j] = n * sw . f;
	}

	TB[0][0] = 0;
	for(int i=1; i<m+1; i++)
		TB[i][0] = 1;
         	
         		for(int j=1; j<n+1; j++)
				TB[0][j] = -1;

	SM[1][0] = sw . U;
	for (int  i = 2; i < m + 1; i++ )
		SM[i][0] = SM[i - 1][0] + sw . E;
	
	SM[0][1] = sw . U;
	for ( int j = 2; j < n + 1; j++ )
		SM[0][j] = SM[0][j - 1] + sw . E;

	double prob = 1.0/profileA->size();

	for(int i=0; i<profileA->size(); i++)
	{
		for(int j=0; j<m; j++)
		{
			int pos = find(characters->begin(), characters->end(), profileA->at(i)[j]) - characters->begin() ;
			PM[ pos ][ j ] = PM[ pos ][ j ] + prob;
		}
	}

return 0;
}

unsigned int alignmentScore_ag(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, unsigned char ** sequences, int ** &TB,  double ** &SM, double ** &PM, double * &IM, double * &DM, vector<char> * characters, unsigned int calculate_TB)
{

	int m = strlen( ( char * ) profileA->at(0) );
	int n = strlen( ( char * ) profileB->at(0) );

	double u;
	double v;
	double w;

	for(int i=1; i<m+1; i++) 
	{
		
		for(int j=1; j<n+1; j++)
		{
			u = SM[i-1][j-1] + probScore( characters, i, j, PM, profileA, profileB , sw );

			DM[i] = MAX2 ( DM[i - 1] + sw . V, SM[i - 1][j] + sw . U );
			v = DM[i];

			IM[j] = MAX2 ( IM[j - 1] + sw . V, SM[i][j - 1] + sw . U );
			w = IM[j];

			SM[i][j] = MAX3 ( u, v, w );
			
			if( calculate_TB == 1 )
			{
				if( SM[i][j] == u)
				{
					TB[i][j] = 0;
				}
				else if(SM[i][j] == w )
				{
					TB[i][j] = -1;
				}
				else if( SM[i][j] == v)
				{
					TB[i][j] = 1;
				}
			}

		}
	}

	if( SM[m][n] > (*score) )
	{	
		( * rotation ) = i;
		( *score ) = SM[m][n];
	}

return 0;
}
