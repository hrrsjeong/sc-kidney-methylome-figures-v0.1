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

exe_list = []

os.system("mkdir -p DSS_results_altered/sample_level/")

flist = glob.glob("DSS_input_altered/rds_combined/*.rds")
for f in flist:
    #continue
    f_base = os.path.basename(f).split(".")[0]
    rds = f
    cov = "sample_covariate_altered.txt"
    state_twolevel_DML = "DSS_results_altered/sample_level/"+f_base+".healthy_vs_altered.DML.sample_level.txt"
    state_twolevel_DMR = "DSS_results_altered/sample_level/"+f_base+".healthy_vs_altered.DMR.sample_level.txt"
    if os.path.exists(state_twolevel_DML) and os.path.exists(state_twolevel_DMR):
        continue
    cmd = "Rscript dss_by_sample_altered.R "+rds+" "+cov+" "+state_twolevel_DML+" "+state_twolevel_DMR
    exe_list.append(cmd)
    #break

f = "DSS_input_altered/rds_combined_celltype/combined.rds"
f_base = "celltype_combined" #os.path.basename(f).split(".")[0]
rds = f
cov = "sample_covariate_altered.txt"
state_twolevel_DML = "DSS_results_altered/sample_level/"+f_base+".healthy_vs_altered.DML.sample_level.txt"
state_twolevel_DMR = "DSS_results_altered/sample_level/"+f_base+".healthy_vs_altered.DMR.sample_level.txt"
celltype_twolevel_DML = "DSS_results_altered/sample_level/"+f_base+".PT_vs_TAL.DML.sample_level.txt"
celltype_twolevel_DMR = "DSS_results_altered/sample_level/"+f_base+".PT_vs_TAL.DMR.sample_level.txt"
cmd = "Rscript dss_by_sample_combined_celltype_altered.R "+rds+" "+cov+" "+state_twolevel_DML+" "+state_twolevel_DMR+" "+celltype_twolevel_DML+" "+celltype_twolevel_DMR
exe_list.append(cmd)

#dss_by_sample_combined_celltype_altered.R

queue = queue.Queue()
for i in range(6):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

