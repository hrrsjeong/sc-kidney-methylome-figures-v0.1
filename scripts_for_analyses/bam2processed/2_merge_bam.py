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

#os.system("mkdir -p ./mHap/tmp")
#os.system("mkdir -p ./mHap/celltype")
os.system("mkdir -p bismark_pseudobulk_bam_file_fofn/merged")
os.system("mkdir -p bismark_pseudobulk_bam_file_fofn/celltype")
os.system("mkdir -p tmp/merged")
os.system("mkdir -p tmp/celltype")
os.system("mkdir -p bismark_pseudobulk_bam/merged")
os.system("mkdir -p bismark_pseudobulk_bam/celltype")
#os.system("mkdir -p MHB_mHapSuite")

flist = glob.glob(f"bismark_pseudobulk_bam/*__*__*.sorted.bam")
dd_sample_leiden = {}
for sfile in flist:
    sample_id = os.path.basename(sfile).split('.')[0].split('__')[0]
    leiden_cluster = os.path.basename(sfile).split('.')[0].split('__')[1]
    dd_sample_leiden.setdefault(sample_id,{})
    dd_sample_leiden[sample_id].setdefault(leiden_cluster,[])
    dd_sample_leiden[sample_id][leiden_cluster].append(sfile)
print (dd_sample_leiden)

all_bam_list = []
for sample_id in dd_sample_leiden.keys():
    sample_bam_list = []
    for leiden in dd_sample_leiden[sample_id].keys():
        sample_bam_list += dd_sample_leiden[sample_id][leiden]
        with open(f"bismark_pseudobulk_bam_file_fofn/celltype/{sample_id}__{leiden}.fofn",'w') as fout:
            fout.write('\n'.join(dd_sample_leiden[sample_id][leiden])+'\n')
        with open(f"tmp/celltype/{sample_id}__{leiden}.bam_link.sh",'w') as fout:
            fout.write(f"samtools merge -@ 4 bismark_pseudobulk_bam/celltype/{sample_id}__{leiden}.sorted.bam $(cat bismark_pseudobulk_bam_file_fofn/celltype/{sample_id}__{leiden}.fofn)\n")
        systemstr = 'echo ""'
        systemstr += f" && bash tmp/celltype/{sample_id}__{leiden}.bam_link.sh"
        systemstr += f" && samtools index bismark_pseudobulk_bam/celltype/{sample_id}__{leiden}.sorted.bam"
        #systemstr += f" && export LD_LIBRARY_PATH=/home/hjeong/Programs/mHapTools/htslib-1.10.2/lib"
        #systemstr += f" && /home/hjeong/Programs/mHapTools/mhaptools convert -i ./bismark_pseudobulk_bam/celltype/{sample_id}__{leiden}.sorted.bam -c /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -o ./mHap/tmp/{sample_id}__{leiden}.mhap.gz"
        #systemstr += f" && zcat ./mHap/tmp/{sample_id}__{leiden}.mhap.gz | awk '{{if (length($4) > 1) print $0}}' | sort -k1,1 -k2,2n | bgzip > ./mHap/celltype/{sample_id}__{leiden}.mhap.gz"
        #systemstr += f" && tabix -b 2 -e 3 ./mHap/celltype/{sample_id}__{leiden}.mhap.gz"
        #systemstr += f" && java -jar /home/hjeong/Programs/mHapSuite/target/mHapSuite-2.0-jar-with-dependencies.jar MHBDiscovery -mhapPath ./mHap/celltype/{sample_id}__{leiden}.mhap.gz -cpgPath /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -window 4 -r2 0.5 -pvalue 0.05 -outputDir MHB_mHapSuite -tag {sample_id}__{leiden}_mHapSuite"

        exe_list.append(systemstr)

    '''

    all_bam_list += sample_bam_list

    with open(f"bismark_pseudobulk_bam_file_fofn/merged/{sample_id}.fofn",'w') as fout:
        fout.write('\n'.join(sample_bam_list)+'\n')
    with open(f"tmp/merged/{sample_id}.bam_link.sh",'w') as fout:
        fout.write(f"samtools merge -@ 4 bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam $(cat bismark_pseudobulk_bam_file_fofn/merged/{sample_id}.fofn)\n")
    
    systemstr = 'echo ""'
    systemstr += f" && bash tmp/merged/{sample_id}.bam_link.sh"
    systemstr += f" && samtools index bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam"
    systemstr += f" && export LD_LIBRARY_PATH=/home/hjeong/Programs/mHapTools/htslib-1.10.2/lib"
    systemstr += f" && /home/hjeong/Programs/mHapTools/mhaptools convert -i ./bismark_pseudobulk_bam/merged/{sample_id}.sorted.bam -c /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -o ./mHap/tmp/{sample_id}.mhap.gz"
    systemstr += f" && zcat ./mHap/tmp/{sample_id}.mhap.gz | awk '{{if (length($4) > 1) print $0}}' | sort -k1,1 -k2,2n | bgzip > ./mHap/{sample_id}.mhap.gz"
    systemstr += f" && tabix -b 2 -e 3 ./mHap/{sample_id}.mhap.gz"
    systemstr += f" && java -jar /home/hjeong/Programs/mHapSuite/target/mHapSuite-2.0-jar-with-dependencies.jar MHBDiscovery -mhapPath ./mHap/{sample_id}.mhap.gz -cpgPath /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -window 4 -r2 0.5 -pvalue 0.05 -outputDir MHB_mHapSuite -tag {sample_id}_mHapSuite"

    exe_list.append(systemstr)
    '''

queue2 = queue.Queue()
for i in range(8):
    t = ThreadRBH(queue2)
    t.setDaemon(True)
    t.start()
for exes in exe_list:
    queue2.put(exes)
queue2.join()


'''
#combine cell type
os.system("mkdir -p ./mHap/celltype_combined")
exe_list2 = []
flist = glob.glob(f"bismark_pseudobulk_bam/celltype/*__*.sorted.bam")
dd_celltype_bam = {}
for sfile in flist:
    sample_id = os.path.basename(sfile).split('.')[0].split('__')[0]
    leiden_cluster = os.path.basename(sfile).split('.')[0].split('__')[1]
    dd_celltype_bam.setdefault(leiden_cluster,{})
    dd_celltype_bam[leiden_cluster][sfile] = 1
print (dd_celltype_bam)
for leiden_cluster in dd_celltype_bam.keys():
    with open(f"bismark_pseudobulk_bam_file_fofn/celltype/{leiden_cluster}.fofn",'w') as fout:
        fout.write('\n'.join(dd_celltype_bam[leiden_cluster].keys())+'\n')
    with open(f"tmp/celltype/{leiden_cluster}.bam_link.sh",'w') as fout:
        fout.write(f"samtools merge -@ 4 bismark_pseudobulk_bam/celltype/{leiden_cluster}.sorted.bam $(cat bismark_pseudobulk_bam_file_fofn/celltype/{leiden_cluster}.fofn)\n")
    systemstr = 'echo ""'
    systemstr += f" && bash tmp/celltype/{leiden_cluster}.bam_link.sh"
    systemstr += f" && samtools index bismark_pseudobulk_bam/celltype/{leiden_cluster}.sorted.bam"
    systemstr += f" && export LD_LIBRARY_PATH=/home/hjeong/Programs/mHapTools/htslib-1.10.2/lib"
    systemstr += f" && /home/hjeong/Programs/mHapTools/mhaptools convert -i ./bismark_pseudobulk_bam/celltype/{leiden_cluster}.sorted.bam -c /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -o ./mHap/tmp/{leiden_cluster}.mhap.gz"
    systemstr += f" && zcat ./mHap/tmp/{leiden_cluster}.mhap.gz | awk '{{if (length($4) > 1) print $0}}' | sort -k1,1 -k2,2n | bgzip > ./mHap/celltype_combined/{leiden_cluster}.mhap.gz"
    systemstr += f" && tabix -b 2 -e 3 ./mHap/celltype_combined/{leiden_cluster}.mhap.gz"
    systemstr += f" && java -jar /home/hjeong/Programs/mHapSuite/target/mHapSuite-2.0-jar-with-dependencies.jar MHBDiscovery -mhapPath ./mHap/celltype_combined/{leiden_cluster}.mhap.gz -cpgPath /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -window 4 -r2 0.5 -pvalue 0.05 -outputDir MHB_mHapSuite -tag {leiden_cluster}_mHapSuite"
    exe_list2.append(systemstr)

print (exe_list2)

queue3 = queue.Queue()
for i in range(8):
    t = ThreadRBH(queue3)
    t.setDaemon(True)
    t.start()
for exes in exe_list2:
    queue3.put(exes)
queue3.join()


##combined all###
with open(f"bismark_pseudobulk_bam_file_fofn/merged/all.fofn",'w') as fout:
    bam_sample_level_list = glob.glob("bismark_pseudobulk_bam/merged/*.sorted.bam") #bismark_pseudobulk_bam/merged/
    fout.write('\n'.join(bam_sample_level_list)+'\n')
with open(f"tmp/merged/combined_all.bam_link.sh",'w') as fout:
    fout.write(f"samtools merge -@ 16 bismark_pseudobulk_bam/combined_all.sorted.bam $(cat bismark_pseudobulk_bam_file_fofn/merged/all.fofn)\n")

systemstr = 'echo ""'
systemstr += f" && bash tmp/merged/combined_all.bam_link.sh"
systemstr += f" && samtools index bismark_pseudobulk_bam/combined_all.sorted.bam"
systemstr += f" && export LD_LIBRARY_PATH=/home/hjeong/Programs/mHapTools/htslib-1.10.2/lib"
systemstr += f" && /home/hjeong/Programs/mHapTools/mhaptools convert -i ./bismark_pseudobulk_bam/combined_all.sorted.bam -c /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -o ./mHap/tmp/combined_all.mhap.gz"
systemstr += f" && zcat ./mHap/tmp/combined_all.mhap.gz | awk '{{if (length($4) > 1) print $0}}' | sort -k1,1 -k2,2n | bgzip > ./mHap/combined_all.mhap.gz"
systemstr += f" && tabix -b 2 -e 3 ./mHap/combined_all.mhap.gz"
systemstr += f" && java -jar /home/hjeong/Programs/mHapSuite/target/mHapSuite-2.0-jar-with-dependencies.jar MHBDiscovery -mhapPath ./mHap/combined_all.mhap.gz -cpgPath /home/hjeong/Programs/mHapSuite/DB/hg38_CpG.gz -window 4 -r2 0.5 -pvalue 0.05 -outputDir MHB_mHapSuite -tag combined_all_mHapSuite"
os.system(systemstr)
#exe_list.append(systemstr)
####
'''


