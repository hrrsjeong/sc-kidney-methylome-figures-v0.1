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

os.system("mkdir -p DSS_input/rds_combined")
os.system("mkdir -p DSS_input/rds_combined_group")
os.system("mkdir -p DSS_input/rds_combined_excl_SDKZ0048_SDKZ0052")
os.system("mkdir -p DSS_input/rds_combined_group_excl_SDKZ0048_SDKZ0052")
os.system("mkdir -p tmp_sh")

dd_state = {}

with open("sample_covariate.txt","r") as f:
    header = f.readline().strip().split("\t")
    for line in f:
        line_temp = line.strip().split("\t")
        state = line_temp[3]
        #if state == "AKI" or state == "CKD":
        #    state = "disease"
        #else:
        #    state = "control"
        dd_state[line_temp[0]] = state
print (dd_state)
#sys.exit()
'''
cat sample_covariate.txt 
sample  state   age     gender
SDKZ0030        Ref     46      Male
SDKZ0032        CKD     64      Female

head sample_covariate.txt 
Sample  Sample_id       Celltype        State   Age     Gender  Disease_level
EC_SDKZ0033     SDKZ0033        EC      CKD     55      Male    1
IMMU_SDKZ0033   SDKZ0033        IMMU    CKD     55      Male    1
NEU_SDKZ0033    SDKZ0033        NEU     CKD     55      Male    1

'''

flist = glob.glob("DSS_input/rds/*__*.rds")
dd_combined = {}
for f in flist:
    sample = f.split("/")[-1].split(".")[0]
    sample_id,group = sample.split("__")
    dd_combined.setdefault(group,{})
    dd_combined[group][sample] = f
for group in dd_combined:
    with open("tmp_sh/"+group+"_DSS_merge.R","w") as f:
        f.write("library(DSS)\n")
        state_list = [dd_state[x] for x in dd_combined[group]]
        for sample in dd_combined[group]:
            f.write(f'{sample} <- readRDS("{dd_combined[group][sample]}")\n')
        f.write(f'combined_input <- List('+",".join(dd_combined[group].keys())+")\n")
        f.write(f'combined <- combineList(combined_input)\n')
        f.write(f'combined <- sort(combined)\n')
        f.write(f'saveRDS(combined, file = "DSS_input/rds_combined/{group}.rds")\n')
        f.write(f'combined_group <- collapseBSseq(combined, group = c(')
        f.write(",".join([f'"{x}"' for x in state_list]))
        f.write("))\n")
        f.write(f'saveRDS(combined_group, file = "DSS_input/rds_combined_group/{group}.rds")\n')

    exe_list.append("Rscript tmp_sh/"+group+"_DSS_merge.R")

for group in dd_combined:
    with open("tmp_sh/"+group+"_DSS_merge_excl_inflammation.R","w") as f:
        f.write("library(DSS)\n")
        state_list = [dd_state[x] for x in dd_combined[group] if x != "SDKZ0048" and x != "SDKZ0052"]
        sample_list = [x for x in dd_combined[group] if x != "SDKZ0048" and x != "SDKZ0052"]
        for sample in sample_list:
            f.write(f'{sample} <- readRDS("{dd_combined[group][sample]}")\n')
        f.write(f'combined_input <- List('+",".join(sample_list)+")\n")
        f.write(f'combined <- combineList(combined_input)\n')
        f.write(f'combined <- sort(combined)\n')
        f.write(f'saveRDS(combined, file = "DSS_input/rds_combined_excl_SDKZ0048_SDKZ0052/{group}.rds")\n')
        f.write(f'combined_group <- collapseBSseq(combined, group = c(')
        f.write(",".join([f'"{x}"' for x in state_list]))
        f.write("))\n")
        f.write(f'saveRDS(combined_group, file = "DSS_input/rds_combined_group_excl_SDKZ0048_SDKZ0052/{group}.rds")\n')

    exe_list.append("Rscript tmp_sh/"+group+"_DSS_merge_excl_inflammation.R")

queue = queue.Queue()
for i in range(12):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

