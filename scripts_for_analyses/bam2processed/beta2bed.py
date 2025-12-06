import sys,os,glob
import queue,threading

class ThreadRBH(threading.Thread):
    def __init__(self,queue):
        threading.Thread.__init__(self)
        self.queue = queue
    def run(self):
        while True:
            systemstr = self.queue.get()
            print (systemstr)
            os.system(systemstr)
            self.queue.task_done()
#~/Programs/wgbs_tools/wgbstools beta2bed --genome hg38 ../celltype_annotation_freeze_06282024/major_celltype_beta/PT.beta
exe_list = []
#os.system("mkdir -p sample_major_beta2bed/")
os.system("mkdir -p ./major_celltype_beta2bed/")
#major_celltype_beta/POD.beta
#flist = glob.glob("../pseudo_bulk/cytosine_report_by_sample/*.CpG_report.merged_CpG_evidence.cov.gz")
#flist = glob.glob("sample_beta/*.beta")# major_celltype_beta/*.beta")
#flist = glob.glob("sample_*_beta/SDK*__*.beta")
#flist = glob.glob("./bam2pat_altered/SDKZ*__*.sorted.beta")# sample_*_beta/SDK*__*.beta")
flist = glob.glob("./major_celltype_beta/*.beta")#SDKZ*__*.beta")

for f in flist:
    sample = f.split("/")[-1].split(".")[0]
    if os.path.exists("major_celltype_beta2bed/"+sample+".bed"):
        continue
    cmd = ' echo ""'
    cmd += " && ~/Programs/wgbs_tools/wgbstools beta2bed --genome hg38 "+f+" | sort -k1,1 -k2,2n >| major_celltype_beta2bed/"+sample+".bed"
    exe_list.append(cmd)

'''
flist = glob.glob("./celltype_beta_altered/*.beta")# sample_*_beta/SDK*__*.beta")
out_dir = "celltype_beta2bed_altered"
os.system(f"mkdir -p {out_dir}")

for f in flist:
    sample = f.split("/")[-1].split(".")[0]
    cmd = ' echo ""'
    cmd += " && ~/Programs/wgbs_tools/wgbstools beta2bed --genome hg38 "+f+f" >| {out_dir}/"+sample+".bed"
    exe_list.append(cmd)
'''
queue = queue.Queue()
for i in range(16):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

