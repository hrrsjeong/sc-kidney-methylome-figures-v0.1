import sys,glob,os
import sys,glob
import os
import random

import queue,threading

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

input_base = "LOO_major_celltype"
input_idx = f"{input_base}_beta"
marker_input_idx = "markers_" + input_idx

pat_out = "ONT_bam2pat"

os.system(f"mkdir -p {pat_out}")
threads = 8
flist = glob.glob("/weka/labs/kzhang/kzhang-lab-scratch/hjeong/ONT_sciMET/output_minimap2_len1000/*.bam")
print (flist)
#sys.exit()
#/weka/labs/kzhang/kzhang-lab-scratch/hjeong/ONT_sciMET/output_minimap2/hKid_SDKZ0032_ONT.bam
os.system(f"mkdir -p deconv_LOO_ONT")
os.system(f"mkdir -p tmp_deconv")
os.system(f"mkdir -p tmp_sh")
for sfile in flist:
    file_idx = sfile.split('/')[-1].split('.')[0]
    if "hKid" in file_idx:
        sample_idx = file_idx.split('hKid_')[-1].split('_ONT')[0]
        genome = "hg38"
    elif "IMR90" in file_idx:
        sample_idx = "SDKZ0032"
        genome = "hg38"
    elif "mKid" in file_idx:
        sample_idx = "mKidney"
        genome = "mm39"
        continue
    os.system(f"mkdir -p tmp_deconv/{file_idx}")
    command = "echo ''"
    #command += f" && ~/Programs/wgbs_tools/wgbstools bam2pat --genome {genome} -f --out_dir {pat_out} --nanopore --min_cpg 2 --np_thresh 0.8 --threads {threads} {sfile}"
    #sample_idx_list = ["SDKZ0037","SDKZ0057","SDKZ0036"]
    #for sample_idx in sample_idx_list:
    #command += f" && ~/Programs/UXM_deconv/uxm deconv --atlas Atlas_{marker_input_idx}/{sample_idx}.top250.tsv --output deconv_LOO_ONT/{file_idx}.{sample_idx}.deconv.csv --tmp_dir tmp_deconv/{file_idx}/ -v -f --threads 8 -l 4 {pat_out}/{file_idx}.pat.gz"
    command += f" && ~/Programs/UXM_deconv/uxm deconv --atlas Atlas_{marker_input_idx}/SDKZ0052.top250.tsv --output deconv_LOO_ONT/{file_idx}.deconv.csv --tmp_dir tmp_deconv/{file_idx}/ -v -f --threads 8 -l 4 {pat_out}/{file_idx}.pat.gz"
    exe_list.append(command)

    #run_systemstr = 'sbatch --wrap=' + systemstr + command
    #run_systemstr += '" '+sbatch_option+f' -e {cwd}/err_ONT_bam2pat__{file_idx}.log -o {cwd}/out_ONT_bam2pat__{file_idx}.log'
    #fout.write(run_systemstr+'\n')
queue = queue.Queue()
for i in range(8):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

