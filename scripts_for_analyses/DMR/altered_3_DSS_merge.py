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

os.system("mkdir -p DSS_input_altered/rds_combined")
os.system("mkdir -p DSS_input_altered/rds_combined_group")
os.system("mkdir -p DSS_input_altered/rds_combined_celltype")
os.system("mkdir -p DSS_input_altered/rds_combined_group_celltype")
os.system("mkdir -p tmp_sh_altered")

dd_state = {}

with open("sample_covariate_altered.txt","r") as f:
    header = f.readline().strip().split("\t")
    for line in f:
        line_temp = line.strip().split("\t")
        celltype = line_temp[2] # cell state instead of disease state
        celltype_origin = line_temp[3]
        state = line_temp[4]
        sample = line_temp[0].replace("-","_")
        #if state == "AKI" or state == "CKD":
        #    state = "disease"
        #else:
        #    state = "control"
        dd_state[sample] = (celltype_origin,state)
print (dd_state)
#sys.exit()
'''
Sample  Sample_id       Celltype        Celltype_origin Cell_state      State   Age     Gender  Disease_level   Pct_cortex      Age_level
SDKZ0037__TAL-altered   SDKZ0037        TAL-altered     TAL     altered control 53      Female  0       0.95    Mid
SDKZ0030__TAL-altered   SDKZ0030        TAL-altered     TAL     altered control 46      Male    0       0.95    Mid
SDKZ0037__TAL-healthy   SDKZ0037        TAL-healthy     TAL     healthy control 53      Female  0       0.95    Mid
SDKZ0033__PT-healthy    SDKZ0033        PT-healthy      PT      healthy disease 55      Male    1       0.9     Mid

'''

flist = glob.glob("DSS_input_altered/rds/*__*.rds")
dd_combined = {}
for f in flist:
    sample = f.split("/")[-1].split(".")[0]
    sample_id,tmp_group = sample.split("__")
    group = tmp_group.split("-")[0]
    sample = sample.replace("-","_")
    print (group)
    dd_combined.setdefault(group,{})
    #dd_combined[group][sample_id] = f
    dd_combined[group][sample] = f

for group in dd_combined:
    with open("tmp_sh_altered/"+group+"_DSS_merge.R","w") as f:
        f.write("library(DSS)\n")
        state_list = [dd_state[x][1] for x in dd_combined[group]]
        for sample in dd_combined[group]:
            f.write(f'{sample} <- readRDS("{dd_combined[group][sample]}")\n')
        f.write(f'combined_input <- List('+",".join(dd_combined[group].keys())+")\n")
        f.write(f'combined <- combineList(combined_input)\n')
        f.write(f'combined <- sort(combined)\n')
        f.write(f'saveRDS(combined, file = "DSS_input_altered/rds_combined/{group}.rds")\n')
        f.write(f'combined_group <- collapseBSseq(combined, group = c(')
        f.write(",".join([f'"{x}"' for x in state_list]))
        f.write("))\n")
        f.write(f'saveRDS(combined_group, file = "DSS_input_altered/rds_combined_group/{group}.rds")\n')

    exe_list.append("Rscript tmp_sh_altered/"+group+"_DSS_merge.R")

with open("tmp_sh_altered/combined_celltype_DSS_merge.R","w") as f:
    f.write("library(DSS)\n")
    state_list = [dd_state[x][0]+"-"+dd_state[x][1] for x in dd_state]
    for sample in dd_state:
        f_name = dd_combined[sample.split("__")[1].split("_")[0]][sample]
        f.write(f'{sample} <- readRDS("{f_name}")\n')
    f.write(f'combined_input <- List('+",".join(dd_state.keys())+")\n")
    f.write(f'combined <- combineList(combined_input)\n')
    f.write(f'combined <- sort(combined)\n')
    f.write(f'saveRDS(combined, file = "DSS_input_altered/rds_combined_celltype/combined.rds")\n')
    f.write(f'combined_group <- collapseBSseq(combined, group = c(')
    f.write(",".join([f'"{x}"' for x in state_list]))
    f.write("))\n")
    f.write(f'saveRDS(combined_group, file = "DSS_input_altered/rds_combined_group_celltype/combined.rds")\n')

exe_list.append("Rscript tmp_sh_altered/combined_celltype_DSS_merge.R")

queue = queue.Queue()
for i in range(12):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

