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

#include <seqan/graph_msa.h>
#include <seqan/align.h>

using namespace seqan;
typedef Graph<Tree<double> > TGraph;

unsigned int nj(TPOcc ** D, unsigned int n,  unsigned char ** seq, struct TSwitch  sw, int * Rot);

unsigned int progAlignment(TPOcc ** D, unsigned char ** seq, TGraph njTree, struct TSwitch  sw, int * Rot, vector<array<int, 2>> * branchingOrder, unsigned int num_seqs );

unsigned int alignPairs(unsigned char * x, unsigned char * y, struct TSwitch sw, int posX, int posY);

unsigned int alignmentScore(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, int ** &TB, double ** &SM, double ** &PM, vector<char> * characters, unsigned int calculate_TB);

unsigned int alignmentScore_ag(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, double * score , struct TSwitch sw, int i, int * rotation, int ** &TB,  double ** &SM, double ** &PM, double ** &IM, double ** &DM, vector<char> * characters, unsigned int calculate_TB);

unsigned int alignSequences(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, vector<int> * profileAPos, vector<int> * profileBPos, unsigned char ** sequences, int ** &TB);

double probScore( vector<char> * characters,  int i, int j, double ** PM, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB, struct TSwitch sw );

unsigned int alignAllocation_ag( double ** &PM, double ** &SM, double ** &IM, double ** &DM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB, struct TSwitch sw );

unsigned int alignAllocation( double ** &PM, double ** &SM, int ** &TB, vector<char> * characters, vector<unsigned char*> * profileA, vector<unsigned char*> * profileB,  struct TSwitch sw);

int similarity( unsigned char x, unsigned char y, struct TSwitch sw);

unsigned int alignPairs_ag(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw, int ** &TB,  double ** &SM,  double ** &IM,  double ** &DM);

unsigned int alignPairs(vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw, int ** &TB,  double ** &SM );

unsigned int pairAllocation_ag( double ** &SM, int ** &TB,  double ** &IM, double ** &DM, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw);

unsigned int pairAllocation( double ** &SM, int ** &TB, vector<unsigned char *> * profileA, vector<unsigned char *> * profileB , struct TSwitch sw);
