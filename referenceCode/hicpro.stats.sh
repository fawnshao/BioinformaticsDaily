#!/bin/bash
workdir=/home1/04935/shaojf/scratch/Xiaoyu.Seq/NovoGene.2-3
# ln -s ../HiCPro.output/NG3_HiC*/hic_results/stats/NG3_HiC* .
grep -w "valid_interaction" $workdir/HiCPro.output/*/hic_results/stats/*/*_allValidPairs.mergestat 
grep -w "valid_interaction_rmdup" $workdir/HiCPro.output/*/hic_results/stats/*/*_allValidPairs.mergestat 