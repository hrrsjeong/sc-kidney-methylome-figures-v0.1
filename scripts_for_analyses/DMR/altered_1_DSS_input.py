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
os.system("mkdir -p DSS_input_altered/table")
os.system("mkdir -p DSS_input_altered/rds")
#os.system("mkdir -p beta2bed/")

#../celltype_annotation_freeze_06282024/tmp_sample_beta2bed_altered/SDKZ0053__PT-healthy.bed
#flist = glob.glob("../pseudo_bulk/cytosine_report_by_sample/*.CpG_report.merged_CpG_evidence.cov.gz")
#flist = glob.glob("../celltype_annotation_freeze_06282024/sample_major_beta/SDKZ*__*.beta")# major_celltype_beta/*.beta")
flist = glob.glob("../celltype_annotation_freeze_06282024/sample_beta2bed_altered/SDKZ*__*.bed")# major_beta/SDKZ*__*.beta")# major_celltype_beta/*.beta")

for f in flist:
    sample = f.split("/")[-1].split(".")[0]
    if os.path.exists("DSS_input_altered/rds/"+sample+".rds"):
        continue
    cmd = ' echo ""'
    #cmd += " && python ~/Programs/pipelines/misc/make_DSS_input_table_from_beta.py beta2bed/"+sample+".bed DSS_input/table "+sample
    cmd += " && python ~/Programs/pipelines/misc/make_DSS_input_table_from_beta.py "+f+" DSS_input_altered/table "+sample
    cmd += " && Rscript ~/Programs/pipelines/misc/make_table2rds.R DSS_input_altered/table/"+sample+"_Info.txt DSS_input_altered/table/"+sample+"_Cov.txt DSS_input_altered/table/"+sample+"_M.txt DSS_input_altered/rds/"+sample+".rds"
    exe_list.append(cmd)

queue = queue.Queue()
for i in range(24):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

