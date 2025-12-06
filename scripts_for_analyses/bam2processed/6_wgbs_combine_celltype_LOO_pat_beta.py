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
threads = 8
memory = 32
per_cpu_memory = 16
rannum = '{:05}'.format(random.randrange(1, 10**5))
#/weka/labs/kzhang/kzhang-lab-scratch/hjeong/scimet_mKidney/processing/merge_bam/bismark_pseudobulk_bam/celltype/mKidney__c0.sorted.bam
qsub_list = f"wgbstools_combine_pat_LOO_{rannum}.sh"
conda_env = "~/qsub_conda_env/scimet_qsub"

systemstr   = '"cd '+cwd
systemstr  += ' && . ~/.bash_qsub'
systemstr  += ' && conda activate '+conda_env

sbatch_option = f' --ntasks 1 --cpus-per-task {threads} --mem-per-cpu={per_cpu_memory}G'
#####################
#sbatch --wrap=". ~/.bash_qsub && conda activate ~/qsub_conda_env/scimet_qsub && cd /home/hjeong/Projects/restricted/ && touch test2.txt && bismark --help" --ntasks 1 --cpus-per-task 8 --mem-per-cpu=8G

celltype_pat_out = "LOO_celltype_pat"
celltype_beta_out = "LOO_celltype_beta"
os.system(f"mkdir -p {celltype_pat_out}")
os.system(f"mkdir -p {celltype_beta_out}")

flist =glob.glob("bam2pat/*__*.sorted.pat.gz")
celltypes = list(set([x.split('/')[-1].split('__')[1].split('.')[0] for x in flist if "__nan" not in x]))
samples = list(set([x.split('/')[-1].split('__')[0] for x in flist]))

fout = open(qsub_list,'w')
for celltype in celltypes:
    flist =glob.glob(f"bam2pat/*__{celltype}.sorted.pat.gz")
    for sample in samples:
        pat_out = f"{celltype_pat_out}/{celltype}__{sample}__LOO"
        #beta_out = f"{celltype_beta_out}/{celltype}__{sample}__LOO"
        flist_sample_except = [x for x in flist if sample not in x]
        command = f' && ~/Programs/wgbs_tools/wgbstools merge --genome hg38 -f -p {pat_out} {" ".join(flist_sample_except)}'
        command += f" && ~/Programs/wgbs_tools/wgbstools pat2beta -f --threads {threads} --genome hg38 -o {celltype_beta_out} {pat_out}.pat.gz"
        run_systemstr = 'sbatch --wrap=' + systemstr + command
        run_systemstr += '" '+sbatch_option+f' -e {cwd}/err_combine_pat__{celltype}__{sample}__LOO.log -o {cwd}/out_combine_pat__{celltype}__{sample}__LOO.log'
        fout.write(run_systemstr+'\n')
fout.close()
with open(qsub_list,'r') as fp:
    for sJob in fp:
        print (sJob)
        #os.system(sJob)

