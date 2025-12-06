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
per_cpu_memory = 8
rannum = '{:05}'.format(random.randrange(1, 10**5))
#/weka/labs/kzhang/kzhang-lab-scratch/hjeong/scimet_mKidney/processing/merge_bam/bismark_pseudobulk_bam/celltype/mKidney__c0.sorted.bam
qsub_list = f"wgbstools_combine_major_pat_{rannum}.sh"
conda_env = "~/qsub_conda_env/scimet_qsub"

systemstr   = '"cd '+cwd
systemstr  += ' && . ~/.bash_qsub'
systemstr  += ' && conda activate '+conda_env

sbatch_option = f' --ntasks 1 --cpus-per-task {threads} --mem-per-cpu={per_cpu_memory}G'
#####################
#sbatch --wrap=". ~/.bash_qsub && conda activate ~/qsub_conda_env/scimet_qsub && cd /home/hjeong/Projects/restricted/ && touch test2.txt && bismark --help" --ntasks 1 --cpus-per-task 8 --mem-per-cpu=8G

celltype_pat_out = "major_celltype_pat"
celltype_beta_out = "major_celltype_beta"
os.system(f"mkdir -p {celltype_pat_out}")
os.system(f"mkdir -p {celltype_beta_out}")

flist = glob.glob("celltype_pat/*.pat.gz")#TAL__SDKZ0052__LOO.pat.gz")
celltypes = list(set([x.split('/')[-1].split('.')[0] for x in flist if "nan." not in x]))
dd_celltype = {}
for celltype in celltypes:
    if "-" in celltype:
        celltype_base = celltype.split('-')[0]
        dd_celltype.setdefault(celltype_base,[]).append(celltype)
    elif celltype == "dPC":
        dd_celltype.setdefault("PC",[]).append(celltype)
    else:
        dd_celltype.setdefault(celltype,[]).append(celltype)
            
fout = open(qsub_list,'w')
for celltype in dd_celltype:
    celltype_pat_files = [f"celltype_pat/{x}.pat.gz" for x in dd_celltype[celltype]]
    pat_out = f"{celltype_pat_out}/{celltype}"
    #beta_out = f"{celltype_beta_out}/{celltype}__{sample}__LOO"
    if len(celltype_pat_files) == 1:
        command = f' && cp {" ".join(celltype_pat_files)} {pat_out}.pat.gz'
    else:
        command = f' && ~/Programs/wgbs_tools/wgbstools merge --genome hg38 -f -p {pat_out} {" ".join(celltype_pat_files)}'
    command += f" && ~/Programs/wgbs_tools/wgbstools pat2beta -f --threads {threads} --genome hg38 -o {celltype_beta_out} {pat_out}.pat.gz"
    run_systemstr = 'sbatch --wrap=' + systemstr + command
    run_systemstr += '" '+sbatch_option+f' -e {cwd}/err_combine_major_pat__{celltype}.log -o {cwd}/out_combine_major_pat__{celltype}.log'
    fout.write(run_systemstr+'\n')
fout.close()

with open(qsub_list,'r') as fp:
    for sJob in fp:
        print (sJob)
        #os.system(sJob)

