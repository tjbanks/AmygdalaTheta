#!/bin/bash

RUN_NUM="${RUN_NUM:-1}"
echo $RUN_NUM

OUTPUT_DIR=./cases_paper/case1/run{$RUN_NUM} sbatch h1_batch_run_hom.sh
OUTPUT_DIR=./cases_paper/case2/run{$RUN_NUM} sbatch h2_batch_vpsi_run_hom.sh
OUTPUT_DIR=./cases_paper/case3/run{$RUN_NUM} sbatch h3_batch_vpsi_run_hom_ach_high.sh
OUTPUT_DIR=./cases_paper/case4/run{$RUN_NUM} sbatch h4_batch_vpsi_run_hom_ach_med.sh
OUTPUT_DIR=./cases_paper/case5/run{$RUN_NUM} sbatch h5_batch_vpsi_run_hom_ach_high_no_vpsi.sh
OUTPUT_DIR=./cases_paper/case6/run{$RUN_NUM} sbatch h6_batch_vpsi_run_hom_ach_med_no_vpsi.sh
