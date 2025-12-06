import sys,glob,os
import sys,glob
import os
import random

import argparse
import pandas as pd
#import Bio.SeqIO
#import Bio.Seq
#import collections
import gzip
import numpy as np
import os
#import pysam
import re

###################
cwd = os.getcwd()
threads = 4
memory = 32
per_cpu_memory = 16
rannum = '{:05}'.format(random.randrange(1, 10**5))
#/weka/labs/kzhang/kzhang-lab-scratch/hjeong/scimet_mKidney/processing/merge_bam/bismark_pseudobulk_bam/celltype/mKidney__c0.sorted.bam
qsub_list = f"wgbstools_bam2pat_{rannum}.sh"
conda_env = "~/qsub_conda_env/scimet_qsub"

systemstr   = '"cd '+cwd
systemstr  += ' && . ~/.bash_qsub'
systemstr  += ' && conda activate '+conda_env

sbatch_option = f' --ntasks 1 --cpus-per-task {threads} --mem-per-cpu={per_cpu_memory}G'
#####################
#sbatch --wrap=". ~/.bash_qsub && conda activate ~/qsub_conda_env/scimet_qsub && cd /home/hjeong/Projects/restricted/ && touch test2.txt && bismark --help" --ntasks 1 --cpus-per-task 8 --mem-per-cpu=8G

#reference_fasta = '/home/hjeong/Projects/ref/genome.fa'
#reference_dir = '/'.join(reference_fasta.split('/')[0:-1])
#reference_fasta_fai = reference_fasta+'.fai'
#cytosine_SNP_file = f"/home/hjeong/Projects/ref/dbSnp155Common.snp.bed"
#if not os.path.exists(reference_fasta_fai):
#    os.system(f"samtools faidx {reference_fasta}")
#BISMARK_REF = "/".join(reference_fasta.split('/')[0:-1])+'/' #"/mnt/nfs/home/Projects/sciMet/ref/"
#min_map_score = 10

#context_file_path = "bismark_context_by_sample"
#cytosine_report_dir = "cytosine_report_by_sample"
#coverage_file_dir = "bismark_cov_by_sample"
#blacklist_bed = "/weka/hjeong/Projects/ref/Blacklist/lists/mm39.excluderanges.bed"
blacklist_bed = "/weka/hjeong/Projects/ref/Blacklist/lists/hg38-blacklist.v2.bed"
pat_out = "bam2pat"
os.system(f"mkdir -p {pat_out}")
#flist =glob.glob("/weka/labs/kzhang/kzhang-lab-scratch/hjeong/scimet_mKidney/processing/merge_bam/bismark_pseudobulk_bam/celltype/mKidney__c*.sorted.bam")
flist =glob.glob("/weka/labs/kzhang/kzhang-lab-scratch/hjeong/scimet/processing/celltype_annotation_freeze_06282024/bismark_pseudobulk_bam/celltype/*__*.sorted.bam")
os.system(f"mkdir -p tmp_sh")
with open(qsub_list,'w') as fout:
    for sfile in flist:
        file_idx = sfile.split('/')[-1].split('.')[0]
        command = ''
        #command += f" && ~/Programs/wgbs_tools/wgbstools bam2pat --genome hg38 -f --out_dir {pat_out} --threads {threads} --blacklist {blacklist_bed} {sfile}"
        command += f" && ~/Programs/wgbs_tools/wgbstools bam2pat --genome hg38 -f --out_dir {pat_out} --threads {threads} {sfile}"
        run_systemstr = 'sbatch --wrap=' + systemstr + command
        run_systemstr += '" '+sbatch_option+f' -e {cwd}/err_bam2pat__{file_idx}.log -o {cwd}/out_bam2pat__{file_idx}.log'
        fout.write(run_systemstr+'\n')

with open(qsub_list,'r') as fp:
    for sJob in fp:
        print (sJob)
        #os.system(sJob)

