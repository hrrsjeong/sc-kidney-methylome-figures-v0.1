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
            print (qcommand)
            os.system(qcommand)
            self.queue.task_done()


exe_list = []
#beta_dir = "LOO_celltype_beta"
cwd = os.getcwd()
types = ["LOO_major_celltype_beta","LOO_celltype_beta"]
types = ["LOO_major_celltype_beta"]#,"LOO_celltype_beta"]
for beta_dir in types:
    #beta_dir = "LOO_major_celltype_beta"# celltype_beta"
    if "major" in beta_dir: # celltype_beta":
        block_path = f"{cwd}/hKidney_sciMET_segment.major_celltype.block.bed"
    else:
        block_path = f"{cwd}/hKidney_sciMET_segment.celltype.block.bed"
    rep_cell = "PT-S1"
    if "major" in beta_dir:
        rep_cell = "PT"
 
    #block_path = "~/Projects/data/Loyfer_2022/data/GSE186458_blocks.s207.hg38.bed.gz"
    thread_per_job = 4
    coverage_per_cpg = 3
    LOO_samples = list(set([x.split('/')[-1].split('__')[1].split('_')[0] for x in glob.glob(f"{beta_dir}/{rep_cell}__*__LOO.beta")]))
    cwd = os.getcwd()
    marker_path = cwd+f"/markers_{beta_dir}/" #celltype_combined_LOO/"
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

        cpg_cov_list = range(5,30)
        for cpg_cov in cpg_cov_list:
            if len(glob.glob(f"{marker_path}/{LOO_sample}/cpg{cpg_cov}/Markers.*.bed"))>0:
                #DMR_LOO_major_celltype_beta/SDKZ0030/cpg29/Markers.B.bed
                continue
            os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg{cpg_cov}")
            systemstr = f"cd {marker_path}/{LOO_sample}/cpg{cpg_cov} && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --tg_quant 0.01 --bg_quant 0.01 --na_rate_tg 0 --na_rate_bg 0.15 -c {cpg_cov*coverage_per_cpg} --min_cpg {cpg_cov} --max_cpg {cpg_cov} --min_bp 50 --max_bp 2000 --pval 0.0001"
            exe_list.append(systemstr)
        if len(glob.glob(f"{marker_path}/{LOO_sample}/cpg30above/Markers.*.bed"))>0:
            continue
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg30above")
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg30above && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --tg_quant 0.01 --bg_quant 0.01 --na_rate_tg 0 --na_rate_bg 0.15 -c 150 --min_cpg 31 --min_bp 50 --max_bp 2000 --pval 0.0001"
        exe_list.append(systemstr)

        '''

        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg5_cov15")
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg8_cov25")
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg11_cov35")
        os.system(f"mkdir -p {marker_path}/{LOO_sample}/cpg15_cov45")
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg5_cov15 && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --na_rate_tg 0 --na_rate_bg 0.15 -c 15 --min_cpg 4 --max_cpg 6 --min_bp 50 --max_bp 1500 --pval 0.0001"
        exe_list.append(systemstr)
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg8_cov25 && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --na_rate_tg 0 --na_rate_bg 0.15 -c 25 --min_cpg 7 --max_cpg 9 --min_bp 50 --max_bp 1500 --pval 0.0001"
        exe_list.append(systemstr)
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg11_cov35 && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --na_rate_tg 0 --na_rate_bg 0.15 -c 35 --min_cpg 10 --max_cpg 12 --min_bp 50 --max_bp 1500 --pval 0.0001"
        exe_list.append(systemstr)
        systemstr = f"cd {marker_path}/{LOO_sample}/cpg15_cov45 && ~/Programs/wgbs_tools/wgbstools find_markers --blocks_path {block_path} --threads {thread_per_job} --groups_file {group_file} --betas {betas_glob} --delta_quants 0.35 --delta_means 0.3 --na_rate_tg 0 --na_rate_bg 0.15 -c 45 --min_cpg 13 --min_bp 50 --max_bp 1500 --pval 0.0001"
        exe_list.append(systemstr)
        '''

        '''
        os.system("mkdir -p group_major_celltype_LOO/"+LOO_sample)
        with open("group_major_celltype_LOO/"+LOO_sample+"/group.csv","w") as f:
            f.write("name,group\n")
            for bed in glob.glob("./LOO_major_celltype_beta/*__"+LOO_sample+"__LOO.beta"):
                group = bed.split("/")[-1].split("__")[0]
                name = bed.split("/")[-1].split(".")[0]
                f.write(name+","+group+"\n")
        os.system("mkdir -p markers_major_celltype_combined_LOO/"+LOO_sample)
        os.system("wgbstools find_markers --blocks_path ~/Projects/data/Loyfer_2022/data/GSE186458_blocks.s207.hg38.bed.gz --threads 48 --groups_file group_major_celltype_LOO/"+LOO_sample+"/group.csv --betas LOO_major_celltype_beta/*__"+LOO_sample+"__LOO.beta --delta_quants 0.35 --delta_means 0.3 --threads 48 --only_hypo --na_rate_tg 0 --na_rate_bg 0.15 -c 10 --min_cpg 4 --min_bp 50 --max_bp 1500 --pval 0.0001")
        os.system("mv Markers.*.bed ./markers_major_celltype_combined_LOO/"+LOO_sample+"/")
        '''
queue = queue.Queue()
for i in range(16):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

#cat <(echo "name,group") <(ls wgbstools_input_bed_celltype/*.bed | grep -v tmp | cut -f 1 -d "." | awk '{print $0","$0}') >| group_celltype/group.csv

#~/Programs/wgbs_tools/wgbstools bed2beta --genome hg38 wgbstools_input_bed_celltype/*.bed
#mkdir -p beta_celltype
#mv *.beta beta_celltype/
#~/Programs/wgbs_tools/wgbstools segment --threads 48 --betas beta_celltype/*.beta -o hKidney_celltype_segment.block.bed
#~/Programs/wgbs_tools/wgbstools beta_to_table --betas beta_celltype/*.beta -c 3 --threads 32 hKidney_celltype_segment.block.bed >| hKidney_celltype_segment.methyl.txt

#mkdir -p markers_scimet_segment
#~/Programs/wgbs_tools/wgbstools find_markers --blocks_path hKidney_celltype_segment.block.bed --groups_file group_celltype/group.csv --betas beta_celltype/*beta --delta_quants .2 --threads 32 --na_rate_tg 0 --na_rate_bg 0.15 -c 4 --min_cpg 3 --min_bp 50 --max_bp 20000 --pval 1
#mv Markers.*.bed ./markers_celltype/

#~/data/net/home/Program/wgbs_tools/wgbstools find_markers --blocks_path sciMETv2_celltype_segment.block.bed --groups_file group/group.csv --betas beta/*beta --delta_quants .25 --threads 32 --na_rate_tg 0 --na_rate_bg 0.1 -c 4 --min_cpg 3 --min_bp 100 --max_bp 20000 --pval 0.005 --only_hyper --top 100
#mv Markers.*.bed ./markers_scimet_segment_hyper_top100/

#~/data/net/home/Program/wgbs_tools/wgbstools find_markers --blocks_path sciMETv2_celltype_segment.block.bed --groups_file group/group.csv --betas beta/*beta --delta_quants .25 --threads 32 --na_rate_tg 0 --na_rate_bg 0.1 -c 4 --min_cpg 3 --min_bp 100 --max_bp 20000 --pval 0.005 --only_hypo --top 200
#mkdir -p markers_scimet_segment_hypo_top200
#mv Markers.*.bed ./markers_scimet_segment_hypo_top200/

#~/data/net/home/Program/wgbs_tools/wgbstools find_markers --blocks_path ../../../data/Loyfer_2022/GSE186458_blocks.s207.hg38.autosome.index.bed.gz --groups_file group/group.csv --betas beta/*beta --delta_quants .25 --threads 32 --na_rate_tg 0 --na_rate_bg 0.2 -c 10 --min_cpg 3 --min_bp 50 --pval 0.05
