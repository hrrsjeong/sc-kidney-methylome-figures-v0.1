import sys,os,glob
import gzip
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

'''

/weka/hjeong/Projects/restricted/scimet/processing/celltype_annotation_freeze_06282024/bismark_context_file_by_sample_fofn/DCT__SDKZ0032.fofn
/weka/hjeong/Projects/restricted/scimet/processing/metacell/pseudobulk/../../methyl_calling/Human_sciMET_patient/filtered_bc_bismark2extract/hKidney-patient_AGCTGTCA/CpG_context_hKidney-patient_GTATCCTTAT+TATTACTCAT+AGCTGTCA.deduplicated.txt.gz
/weka/hjeong/Projects/restricted/scimet/processing/metacell/pseudobulk/../../methyl_calling/Human_sciMET_patient/filtered_bc_bismark2extract/hKidney-patient_TACTCATA/CpG_context_hKidney-patient_CGAATCTCCT+TGGCCGGCCT+TACTCATA.deduplicated.txt.gz

'''


sample_idx = "hKidney"
os.system("mkdir -p bismark_pseudobulk_bam_file_fofn/")
os.system("mkdir -p bismark_pseudobulk_bam/sorted")
os.system("mkdir -p bismark_pseudobulk_bam/tmp")
os.system("mkdir -p tmp")
dd = {}

#flist = glob.glob("/home/hjeong/Projects/restricted/scimet/processing/pseudo_bulk/bismark_context_file_by_sample_fofn/*.fofn")#TAL-2_SDKZ0122.fofn
#flist = glob.glob("/home/hjeong/Projects/restricted/scimet_mKidney/processing/pseudo_bulk/bismark_context_file_fofn/mKidney_sciMET_c*_mCG_20kb.fofn")
flist = glob.glob("/home/hjeong/Projects/restricted/scimet/processing/celltype_annotation_freeze_06282024/bismark_context_file_by_sample_fofn/*__*.fofn")
dd_sample = {}
#base_kpmp_mouse_scimet_path = "/home/hjeong/hypahub/scimet_mKidney/sciMET_processed/"
base_kpmp_human_scimet_path = "/home/hjeong/hypahub/scimet_hKidney/sciMET_processed/"
for sfile in flist:
    cnt = 0
    leiden_cluster = os.path.basename(sfile).split('.')[0].split('__')[0]
    sample_id = os.path.basename(sfile).split('.')[0].split('__')[1]
    with open(sfile,'r') as fin:
        for line in fin:
            cellbam_file = base_kpmp_human_scimet_path+'/'.join(line.strip().split('/')[-4:]).replace('/filtered_bc_bismark2extract/','/filtered_bc_rmdup_bam/').replace('.txt.gz','.bam').replace('CpG_context_','')
            ##cell_bc = cellbam_file.split('/')[-1].split('.')[0]
            ##sample_id = dd_bc[cell_bc]
            dd_sample.setdefault(leiden_cluster,{})
            dd_sample[leiden_cluster].setdefault(sample_id,0)
            dd_sample[leiden_cluster][sample_id] += 1
            idx = (dd_sample[leiden_cluster][sample_id] // 100 ) + 1
            dd.setdefault(leiden_cluster,{})
            dd[leiden_cluster].setdefault(sample_id,{})
            dd[leiden_cluster][sample_id].setdefault(idx,[])
            dd[leiden_cluster][sample_id][idx].append(cellbam_file)

for leiden in dd:
    for sample_id in dd[leiden]:
        for idx in dd[leiden][sample_id]:
            if os.path.exists(f"bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam.bai"):
                continue
            with open(f"bismark_pseudobulk_bam_file_fofn/{sample_id}__{leiden}__{idx}.fofn",'w') as fout:
                fout.write('\n'.join(dd[leiden][sample_id][idx])+'\n')
            with open(f"tmp/{sample_id}__{leiden}__{idx}.bam_link.sh",'w') as fout:
                fout.write(f"samtools merge -@ 4 bismark_pseudobulk_bam/tmp/{sample_id}__{leiden}__{idx}.bam $(cat bismark_pseudobulk_bam_file_fofn/{sample_id}__{leiden}__{idx}.fofn)\n")
            systemstr = 'echo ""'
            systemstr += f" && bash tmp/{sample_id}__{leiden}__{idx}.bam_link.sh"
            systemstr += f" && samtools sort -@4 -m 4G -o ./bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam ./bismark_pseudobulk_bam/tmp/{sample_id}__{leiden}__{idx}.bam"
            systemstr += f" && samtools index bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam"
            systemstr += f" && samtools view -c bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam >| bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.bam.readcount.txt"
            systemstr += f" && rm -rf bismark_pseudobulk_bam/tmp/{sample_id}__{leiden}__{idx}.bam"

            ##systemstr += f" && java -jar ~/data/net/home/Program/mHapSuite/mHapSuite-2.0-alpha/target/mHapSuite-2.0-jar-with-dependencies.jar convert --cpgPath ~/data/net/home/Program/mHapSuite/DB/hg38_CpG.gz --inputFile bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam --nonDirectional --outPutFile ./mHap/tmp/{sample_id}__{leiden}__{idx}.mhap.gz"

            ##systemstr += f" && samtools sort -n -m 8G -o ./bismark_pseudobulk_bam/sorted/{sample_id}__{leiden}__{idx}.sorted.bam ./bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.bam"
            ##systemstr += f" && bismark_methylation_extractor --comprehensive --merge_non_CpG --no_header --gzip --no_overlap --buffer_size 12G -o ./bismark_pseudobulk_context/tmp/ ./bismark_pseudobulk_bam/sorted/{sample_id}__{leiden}__{idx}.sorted.bam"
            ##systemstr += f" && zcat ./bismark_pseudobulk_context/tmp/CpG_context_{sample_id}__{leiden}__{idx}.sorted.txt.gz | wc -l >| ./bismark_pseudobulk_context/tmp/CpG_context_{sample_id}__{leiden}__{idx}.sorted.txt.gz.readcount"
            ###systemstr += f" && python split_context.py ./bismark_pseudobulk_context/tmp/CpG_context_{sample_id}__{leiden}__{idx}.sorted.txt.gz.readcount"
      
            exe_list.append(systemstr)
#bismark_pseudobulk_context/tmp/CpG_context_TAL__6.sorted.txt.gz

queue = queue.Queue()
for i in range(16):
    t = ThreadRBH(queue)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue.put(exes)
queue.join()

'''
exe_list = []
os.system("mkdir -p ./mHap/celltype/tmp2")
os.system("mkdir -p ./mHap/tmp2")
os.system("mkdir -p ./MHB_discovery/celltype")
os.system("mkdir -p tmp/merged")
os.system("mkdir -p bismark_pseudobulk_bam/merged")
#flist = glob.glob("./mHap/tmp/*__*__*.mhap.gz")
flist = glob.glob("bismark_pseudobulk_bam/{sample_id}__{leiden}__{idx}.sorted.bam")
dd_sample_leiden = {}
for sfile in flist:
    sample_id = os.path.basename(sfile).split('.')[0].split('__')[0]
    leiden_cluster = os.path.basename(sfile).split('.')[0].split('__')[1]
    dd_sample_leiden.setdefault(sample_id,{})
    dd_sample_leiden[sample_id].setdefault(leiden_cluster,[])
    dd_sample_leiden[sample_id][leiden_cluster].append(sfile)

for sample_id in dd_sample_leiden.keys():
    sample_bam_list = []
    for leiden in dd_sample_leiden[sample_id].keys():
        sample_bam_list += dd_sample_leiden[sample_id][leiden]

    with open(f"bismark_pseudobulk_bam_file_fofn/merged/{sample_id}.fofn",'w') as fout:
        fout.write('\n'.join(sample_bam_list)+'\n')
    with open(f"tmp/merged/{sample_id}.bam_link.sh",'w') as fout:
        fout.write(f"samtools merge -@ 4 bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam $(cat bismark_pseudobulk_bam_file_fofn/merged/{sample_id}.fofn)\n")
    
    systemstr = 'echo ""'
    systemstr += f" && bash tmp/merged/{sample_id}.bam_link.sh"
    systemstr += f" && samtools index bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam"
    systemstr += f" && java -jar ~/data/net/home/Program/mHapSuite/mHapSuite-2.0-alpha/target/mHapSuite-2.0-jar-with-dependencies.jar convert --cpgPath ~/data/net/home/Program/mHapSuite/DB/hg38_CpG.gz --inputFile bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam --outPutFile ./mHap/{sample_id}.mhap.gz"

    systemstr += f" && zcat ./mHap/tmp2/{sample_id}.mhap.gz | sort -k1,1 -k2,2n | bgzip >| ./mHap/{sample_id}.mhap.gz && tabix -b 2 -e 3 ./mHap/{sample_id}.mhap.gz"
    systemstr += f" && tabix -b 2 -e 3 ./mHap/{sample_id}.mhap.gz"
    systemstr += f" && java -jar ~/data/net/home/Program/mHapSuite/mHapSuite-2.0-alpha/target/mHapSuite-2.0-jar-with-dependencies.jar MHBDiscovery -mhapPath ./mHap/{sample_id}.mhap.gz -cpgPath ~/data/net/home/Program/mHapSuite/DB/hg38_CpG.gz -window 5 -r2 0.5 -pvalue 0.05 -outputDir MHB_discovery -tag {sample_id}"
    exe_list.append(systemstr)

queue2 = queue2.Queue()
for i in range(4):
    t = ThreadRBH(queue2)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue2.put(exes)
queue2.join()
'''


