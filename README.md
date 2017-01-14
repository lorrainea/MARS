MARS: Multiple circular sequence Alignment using Refined Sequences
===


<b>Description</b>: A fundamental assumption of all widely-used multiple sequence alignment techniques is that the left- and right-most positions of the input sequences are relevant to the alignment. However, the position where a sequence starts or ends can be totally arbitrary due to a number of reasons: arbitrariness in the linearisation (sequencing) of a circular molecular structure; or inconsistencies introduced into sequence databases due to different linearisation standards. These scenarios are relevant, for instance, in the process of multiple sequence alignment of mitochondrial DNA, viroid, viral or other genomes, which have a circular molecular structure. 

MARS is a program, which can be used in conjunction with any multiple sequence alignment program, to address this problem effectively and efficiently.

<b>Installation</b>: To compile MARS, please follow the instructions given in file INSTALL.

<b>INPUT</b>: A set of sequences in FASTA format. The input file is specified using the <b>-i</b> option. Typical input sequences include mitochondrial DNA, viroid, viral or other circular genomes. 

<b>OUTPUT</b>: The set of refined (cyclically shifted) sequences with no gaps added in FASTA format. The output file is specified using the <b>-o</b> option. This output file can then be used as input to the preferred MSA program to obtain the final alignment.

```
 Usage: mars <options>
 Standard (Mandatory):
  -a, --alphabet              <str>     'DNA' for nucleotide  sequences  or 'PROT' for protein  sequences.
  -i, --input-file            <str>     MultiFASTA input filename.
  -o, --output-file           <str>     Output filename with rotated sequences.
 Optional:
 Cyclic Edit Distance Computation.
  -m, --method                <int>     0 for heuristic Cyclic Edit Distance (hCED) - Faster but less accurate. 
                                        1 for branch and bound method - Slower but exact. Default: 0.
  -q, --q-length              <int>     The q-gram length for method hCED. Default: 5.
 Refinement Parameters. 
  -l, --block-length          <int>     The length of each block. Default: sqrt(seq_length).
  -P, --refine-blocks         <dbl>     Refine the rotations by aligning P blocks of the ends. Default: 1.
 Gap Penalties.
  -O, --gap-open-seq          <int>     Gap open penalty in pairwise block alignment. Default: -10.
  -E, --gap-extend-seq        <int>     Gap extension penalty in pairwise block alignment. Default: -1.
  -U, --gap-open-pro          <int>     Gap open penalty in alignment of profiles. Default: -10.
  -V, --gap-extend-pro        <int>     Gap extension penalty in alignment of profiles. Default: -1.
 Number of threads.
  -T, --threads               <int>     Number of threads to use. Default: 1. 
```

<b>Example</b>: For a typical run, see file EXAMPLES.

<b>Citation</b>
```
L. A. K. Ayad, S. P. Pissis
MARS: improving multiple circular sequence alignment using refined sequences
BMC Genomics, vol. 18, no. 1, 2017, pp. 86.
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2016 Lorraine A.K. Ayad and Solon P. Pissis.
