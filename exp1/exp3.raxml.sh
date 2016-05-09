#!/bin/bash

 ./raxmlHPC-PTHREADS -m GTRCAT -n 50.2500.35.true.raxml -p 12345 -s 50.2500.35.true
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS -m GTRCAT -n 50.2500.35.mars.raxml -p 12345 -s 50.2500.35.mars
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS -m GTRCAT -n 50.2500.35.orig.raxml -p 12345 -s 50.2500.35.orig
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS -m GTRCAT -n 50.2500.35.bear.raxml -p 12345 -s 50.2500.35.bear
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*
 ./raxmlHPC-PTHREADS -m GTRCAT -n 50.2500.35.cyc.raxml -p 12345 -s 50.2500.35.cyc
 rm RAxML_i*
 rm RAxML_l*
 rm RAxML_p*
 rm RAxML_r*

 cat RAxML_bestTree.50.2500.35.true.raxml RAxML_bestTree.50.2500.35.mars.raxml RAxML_bestTree.50.2500.35.orig.raxml RAxML_bestTree.50.2500.35.bear.raxml RAxML_bestTree.50.2500.35.cyc.raxml > 50.2500.35.allTrees
 ./raxmlHPC-PTHREADS -f r -z 50.2500.35.allTrees -m GTRCAT -n 50.2500.35
 rm RAxML_i*
