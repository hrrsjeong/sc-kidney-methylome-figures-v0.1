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

os.system("mkdir -p DSS_results_altered/group_level/")

flist = glob.glob("DSS_input_altered/rds_combined_group/*.rds")
for f in flist:
    f_base = os.path.basename(f).split(".")[0]
    rds = f
    state_twolevel_DML = "DSS_results_altered/group_level/"+f_base+".healthy_vs_altered.DML.group_level.txt"
    state_twolevel_DMR = "DSS_results_altered/group_level/"+f_base+".healthy_vs_altered.DMR.group_level.txt"
    cov_filtered = "DSS_results_altered/group_level/"+f_base+".healthy_vs_altered.cov_filtered.group_level.txt"
    if os.path.exists(state_twolevel_DML) and os.path.exists(state_twolevel_DMR):
        continue
    cmd = "Rscript dss_by_group_altered.R "+rds+" "+state_twolevel_DML+" "+state_twolevel_DMR+" "+cov_filtered
    exe_list.append(cmd)
    #break

f = "DSS_input_altered/rds_combined_group_celltype/combined.rds" 
f_base = "celltype_combined" #os.path.basename(f).split(".")[0]
rds = f
healthy_twolevel_DML = "DSS_results_altered/group_level/"+f_base+".PT-healthy_vs_TAL-healthy.DML.group_level.txt"
healthy_twolevel_DMR = "DSS_results_altered/group_level/"+f_base+".PT-healthy_vs_TAL-healthy.DMR.group_level.txt"
altered_twolevel_DML = "DSS_results_altered/group_level/"+f_base+".PT-altered_vs_TAL-altered.DML.group_level.txt"
altered_twolevel_DMR = "DSS_results_altered/group_level/"+f_base+".PT-altered_vs_TAL-altered.DMR.group_level.txt"
cmd = "Rscript dss_by_group_combined_celltype_altered.R "+rds+" "+healthy_twolevel_DML+" "+healthy_twolevel_DMR+" "+altered_twolevel_DML+" "+altered_twolevel_DMR
exe_list.append(cmd)

queue = queue.Queue()
for i in range(5):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

