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
#include <seqan/basic.h>
#include <vector>
#include "mars.h"
#include "nj.h"

using namespace seqan;
using namespace std;

unsigned int nj(TPOcc ** D, unsigned int n, unsigned char ** seq, struct TSwitch  sw, int * Rot )
{
	vector<array<int, 2>> * branchingOrder = new vector<array<int, 2>>();

	TGraph njTreeOut;
	int a =0;
	String<double> mat;
	resize(mat, n*n);

	for(int i=0; i<n; i++)
	{
		for(int j=0; j<n; j++)
		{
			assignValue(mat, a, D[i][j]. err);
			a++;
		}
	}

	njTree(mat, njTreeOut);

	typedef Iterator<TGraph, EdgeIterator>::Type TEdgeIter;
	TEdgeIter edIt(njTreeOut);
	array<int, 2> children;
	for(;!atEnd(edIt);goNext(edIt)) {
		children[0] = (int) targetVertex(edIt);
		goNext(edIt);
		children[1] = (int) targetVertex(edIt);
		branchingOrder->push_back( children );
	}

	fprintf ( stderr, " Starting progressive alignment\n" );

	/*Progressively aligns sequences using refined sequences*/
	progAlignment( D, seq , njTreeOut, sw , Rot , branchingOrder, n );
	
	delete( branchingOrder );
	
	return 0;

}
