import sys,os,glob
import queue,threading
import argparse
import pandas as pd
#import Bio.SeqIO
#import Bio.Seq
import re

class ThreadRBH(threading.Thread):
    def __init__(self,queue):
        threading.Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            qcommand = self.queue.get()
            #print (qcommand)
            os.system(qcommand)
            self.queue.task_done()

exe_list = []
#beta_dir = "LOO_celltype_beta"
types = ["LOO_major_celltype_beta","LOO_celltype_beta"]
cwd = os.getcwd()
for beta_dir in types:
    #beta_dir = "LOO_major_celltype_beta"# celltype_beta"
    if "major" in beta_dir: # celltype_beta":
        block_path = f"{cwd}/hKidney_sciMET_segment.major_celltype.block.bed"
    else:
        block_path = f"{cwd}/hKidney_sciMET_segment.celltype.block.bed"
    rep_cell = "PT-S1"
    if "major" in beta_dir:
        rep_cell = "PT"
    thread_per_job = 2
    coverage_per_cpg = 3
    LOO_samples = list(set([x.split('/')[-1].split('__')[1].split('_')[0] for x in glob.glob(f"{beta_dir}/{rep_cell}__*__LOO.beta")]))
    cwd = os.getcwd()
    marker_path = cwd+f"/DMR_{beta_dir}/" #celltype_combined_LOO/"
    os.system("mkdir -p "+marker_path)
    for LOO_sample in LOO_samples:
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/")
        group_file = f"{marker_path}/{LOO_sample}/group.csv"
        with open(group_file,"w") as f:
            f.write("name,group\n")
            betas_glob = f"{cwd}/{beta_dir}/*__"+LOO_sample+"__LOO.beta"
            for bed in glob.glob(betas_glob): #./wgbstools_input_bed_celltype_combined_subtypes
                group = bed.split("/")[-1].split("__")[0]
                name = bed.split("/")[-1].split(".")[0]
                f.write(name+","+group+"\n")
        cpg_cov_list = range(4,30)
        for cpg_cov in cpg_cov_list:
            if len(glob.glob(f"{marker_path}/{LOO_sample}/cpg{cpg_cov}/Markers.*.bed"))>0:
                #DMR_LOO_major_celltype_beta/SDKZ0030/cpg29/Markers.B.bed
                continue
            os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg{cpg_cov}")
            systemstr = f"cd {marker_path}/{LOO_sample}/cpg{cpg_cov} && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.15 --delta_means 0.25 --tg_quant 0.01 --bg_quant 0.01 --na_rate_tg 0 --na_rate_bg 0.15 -c {cpg_cov*coverage_per_cpg} --min_cpg {cpg_cov} --max_cpg {cpg_cov} --min_bp 50 --max_bp 2000 --pval 1"
            exe_list.append(systemstr)
        if len(glob.glob(f"{marker_path}/{LOO_sample}/cpg30above/Markers.*.bed"))>0:
            continue
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg30above")
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg30above && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.15 --delta_means 0.25 --tg_quant 0.01 --bg_quant 0.01 --na_rate_tg 0 --na_rate_bg 0.15 -c 150 --min_cpg 31 --min_bp 50 --max_bp 2000 --pval 1"
        exe_list.append(systemstr)

queue = queue.Queue()
for i in range(30):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

tmp_flist = glob.glob(f'DMR_*/*/*/Markers.*.bed')

for beta_dir in types:
    #LOO_samples = list(set([x.split('/')[-1].split('__')[1].split('_')[0] for x in glob.glob(f"{beta_dir}/{rep_cell}__*__LOO.beta")]))
    all_DMR_files = glob.glob(f"DMR_{beta_dir}/*/*/Markers.*.bed")
    LOO_samples = list(set([x.split('/')[-3] for x in all_DMR_files]))
    celltypes = list(set([x.split('/')[-1].split(".")[1] for x in all_DMR_files]))
    LOO_sample_cnt = len(LOO_samples)
    print (celltypes)
    print (LOO_samples)

    dd = {}
    for LOO_sample in LOO_samples:
        LOO_sample_files = [x for x in all_DMR_files if LOO_sample in x]
        for DMR_file in LOO_sample_files:
            celltype = DMR_file.split("/")[-1].split(".")[1]
            dd.setdefault(celltype,{})
            with open(DMR_file,'r') as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    line_temp = line.strip().split('\t')
                    ttest = float(line_temp[-2])
                    if ttest>0.0005:
                        continue
                    region = line_temp[6]
                    dd[celltype].setdefault(region,{})
                    dd[celltype][region][LOO_sample] = line_temp
    for LOO_sample in LOO_samples:
        for celltype in dd:
            LOO_sample_file = f"consensus_DMR_{beta_dir}/{LOO_sample}/{celltype}.bed"
            os.system(f"mkdir -p consensus_DMR_{beta_dir}/{LOO_sample}/")
            with open(LOO_sample_file,"w") as f:
                f.write("chr\tstart\tend\tstartCpG\tendCpG\ttarget\tregion\tlenCpG\tbp\ttg_mean\tbg_mean\tdelta_means\tdelta_quants\tdelta_maxmin\tttest\tdirection\n")
                for region in dd[celltype]:
                    sample_cnt = len(dd[celltype][region])
                    if sample_cnt<LOO_sample_cnt:
                        continue
                    f.write("\t".join(dd[celltype][region][LOO_sample])+"\n")


        



'''
head DMR_LOO_major_celltype_beta/SDKZ0030/cpg10/Markers.CNT.bed 
#chr    start   end     startCpG        endCpG  target  region  lenCpG  bp      tg_mean bg_mean delta_means     delta_quants    delta_maxmin    ttest   direction
chr6    167295542       167295624       10763880        10763890        CNT     chr6:167295542-167295624        10CpGs  82bp    0.143   0.728   0.585   0.245   0.224  8.29e-08 U
chr9    130111930       130112095       14844643        14844653        CNT     chr9:130111930-130112095        10CpGs  165bp   0.0623  0.732   0.669   0.277   0.276  1.34e-06 U
chr15   31040094        31040223        20974303        20974313        CNT     chr15:31040094-31040223 10CpGs  129bp   0.293   0.901   0.608   0.204   0.183   8.38e-09U
'''




