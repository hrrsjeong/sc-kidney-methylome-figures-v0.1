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

os.system("mkdir -p DSS_results/group_level/")
os.system("mkdir -p DSS_results_excl_inflammation/group_level/")
dd_state = {}

with open("sample_covariate.txt","r") as f:
    header = f.readline().strip().split("\t")
    for line in f:
        line_temp = line.strip().split("\t")
        dd_state[line_temp[0]] = line_temp[1]

'''
cat sample_covariate.txt 
sample  state   age     gender
SDKZ0030        Ref     46      Male
SDKZ0032        CKD     64      Female
'''

flist = glob.glob("DSS_input/rds_combined_group/*.rds")
for f in flist:
    f_base = os.path.basename(f).split(".")[0]
    rds = f
    Disease_twolevel_DML = "DSS_results/group_level/"+f_base+".Disease_DML.group_level.txt"
    Disease_twolevel_DMR = "DSS_results/group_level/"+f_base+".Disease_DMR.group_level.txt"
    if os.path.exists(Disease_twolevel_DML) and os.path.exists(Disease_twolevel_DMR):
        continue
    cmd = "Rscript dss_by_group.R "+rds+" "+Disease_twolevel_DML+" "+Disease_twolevel_DMR
    exe_list.append(cmd)
    #break
flist = glob.glob("DSS_input/rds_combined_group_excl_SDKZ0048_SDKZ0052/*.rds")
for f in flist:
    f_base = os.path.basename(f).split(".")[0]
    rds = f
    Disease_twolevel_DML = "DSS_results_excl_inflammation/group_level/"+f_base+".Disease_DML.group_level.txt"
    Disease_twolevel_DMR = "DSS_results_excl_inflammation/group_level/"+f_base+".Disease_DMR.group_level.txt"
    if os.path.exists(Disease_twolevel_DML) and os.path.exists(Disease_twolevel_DMR):
        continue
    cmd = "Rscript dss_by_group.R "+rds+" "+Disease_twolevel_DML+" "+Disease_twolevel_DMR
    exe_list.append(cmd)
    #break

queue = queue.Queue()
for i in range(5):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

