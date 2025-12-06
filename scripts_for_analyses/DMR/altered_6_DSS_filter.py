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
os.system("mkdir -p DSS_results_altered/group_level_filtered/") #celltype_combined.PT-altered_vs_TAL-altered.DMR.group_level.txt")
os.system("mkdir -p DSS_results_altered/sample_level_filtered/")#TAL.healthy_vs_altered.DMR.sample_level.txt")# d 
os.system("mkdir -p tmp_sh")

#cat DSS_results_excl_inflammation/sample_level/PT.Disease_DML.sample_level.txt | awk -vOFS='\t' '{if ($5 < 0.1) print $1,$2-1,$2}' | bedtools intersect -a <(tail -n +2 DSS_results_excl_inflammation/sample_level/PT.Disease_DMR.sample_level.txt) -b - -wa | uniq >| DSS_results_excl_inflammation/sample_level/PT.Disease_DMR.sample_level.filtered_by_DML.txt


flist = glob.glob("DSS_results_altered/*_level/*.DMR.*.txt")
for f in flist:
    f_base = os.path.basename(f).split(".txt")[0]
    level = f.split("/")[-2]
    tmp_dir = '/'.join(f.split("/")[0:-1])+'/tmp/'
    os.system("mkdir -p "+tmp_dir)
    DML = f.replace(".DMR.",".DML.") 
    DMR = f
    DML_corrected = tmp_dir+DML.split("/")[-1].replace(".txt",".corrected.txt")
    DMR_corrected = tmp_dir+f_base+".corrected.txt"
    filtered_DMR = f.replace(".txt",".filtered_by_DML.txt").replace(level+"/",level+"_filtered/")
    filtered_DML = f.replace(".DMR.",".DML.").replace(level+"/",level+"_filtered/")
    print (DML)
    with open(DML,'r') as tmp_finp:
        title = tmp_finp.readline().strip().split("\t")
        if "fdr" in title:
            fdr_idx = title.index("fdr") + 1
        elif "fdrs" in title:
            fdr_idx = title.index("fdrs") + 1
        tmp_finp.seek(0)
    print (fdr_idx)
    with open("tmp_sh/"+f_base+".altered.sh","w") as ouf:
        ouf.write("python dmr_filter.py "+DML+" "+DMR+" "+DML_corrected+" "+DMR_corrected+'\n')
        ouf.write("cat "+DML_corrected+" | awk -vOFS='\\t' '{if ($"+str(fdr_idx)+" < 0.05) print $0}' >| "+filtered_DML+'\n')
        ouf.write("cat "+DML_corrected+" | awk -vOFS='\\t' '{if ($"+str(fdr_idx)+" < 0.2) print $1,$2-1,$2+1}' | bedtools intersect -a "+DMR_corrected+" -b - -wa | uniq >| "+filtered_DMR)
    cmd = "bash tmp_sh/"+f_base+".altered.sh"
    exe_list.append(cmd)
    #break

queue = queue.Queue()
for i in range(12):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

