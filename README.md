MARS: Multiple circular sequence Alignment using Refined Sequences
===


<b>Introduction</b>: A fundamental assumption of all widely-used multiple sequence alignment techniques is that the left- and right-most positions of the input sequences are relevant to the alignment. However, the position where a sequence starts or ends can be totally arbitrary due to a number of reasons: arbitrariness in the linearisation (sequencing) of a circular molecular structure; or inconsistencies introduced into sequence databases due to different linearisation standards. These scenarios are relevant, for instance, in the process of multiple sequence alignment of mitochondrial DNA, viroid, viral or other genomes, which have a circular molecular structure. 

MARS is a program, which can be used in conjunction with any multiple sequence alignment program, to address this problem effectively and efficiently.

<b>Installation</b>: To compile MARS, please follow the instructions given in file INSTALL.

<b>INPUT</b>: A set of sequences in FASTA format. The input file is specified using the <b>-i</b> option. Typical input sequences include mitochondrial DNA, viroid, viral or other genomes, which have a circular molecular structure. 

<b>OUTPUT</b>: The set of refined (cyclically shifted) sequences with no gaps added in FASTA format. The output file is specified using the <b>-o</b> option. This output file can then be used as input to the preferred MSA program to obtain the final alignment.

<b>Options</b>: Type <em>./mars -h</em>

<b>Example</b>: For a typical run, see file EXAMPLES.

<b>License</b>: GNU GPLv3 License; Copyright (C) 2016 Lorraine A.K. Ayad and Solon P. Pissis.

