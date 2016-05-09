/**
    hCED: Heuristic Cyclic Edit Distance
    Copyright (C) 2016 Solon P. Pissis, Lorraine A. K. Ayad 

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

int editDistanceMyers(unsigned char *xInput, unsigned char *yInput, int mInput, int nInput, unsigned int * distance);

unsigned int editDistance(unsigned char * xInput, unsigned char * yInput, int mInput, int nInput, unsigned int * distance, int sub, int ins, int del);

unsigned int sacsc_refinement ( unsigned char * x, unsigned char * xr, unsigned char * y, struct TSwitch  sw, unsigned int * rotation, unsigned int * distance );
