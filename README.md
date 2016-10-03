MARS: Multiple circular sequence Alignment using Refined Sequences
===

GNU GPLv3 License; Copyright (C) 2016 Lorraine A.K. Ayad and Solon P. Pissis.

Introduction: A fundamental assumption of all widely-used multiple sequence alignment techniques is that the left- and right-most positions of the input sequences are relevant to the alignment. However, the position where a sequence starts or ends can be totally arbitrary due to a number of reasons: arbitrariness in the linearisation (sequencing) of a circular molecular structure; or inconsistencies introduced into sequence databases due to different linearisation standards. These scenarios are relevant, for instance, in the process of multiple sequence alignment of mitochondrial DNA, viroid, viral or other genomes, which have a circular molecular structure. A solution for these inconsistencies would be to identify a suitable rotation (cyclic shift) for each sequence; these refined sequences may in turn lead to improved multiple sequence alignments using the preferred multiple sequence alignment program.

MARS is a program, which can be used in conjunction with any multiple sequence alignment program, to address this problem effectively and efficiently.

To compile MARS, please follow the instructions given in file INSTALL.

Input: A set of sequences in FASTA format. The input file is specified using the -i option. 

Output: The set of refined (via cyclic shifts) sequences with no gaps added in FASTA format. The output file is specified using the -o option. This output file can then be used as input of the preferred MSA program to obtain the final alignment.

 
