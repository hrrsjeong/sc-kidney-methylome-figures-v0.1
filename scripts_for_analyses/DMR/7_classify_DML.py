import sys,glob,os


flist = glob.glob('DSS_results*/sample_level_filtered/*.*.sample_level.*txt')

for f in flist:
    dir_idx = '/'.join(f.split('/')[0:-1])+'/classified'
    sfile_idx = f.split('/')[-1].split('.txt')[0]
    if not os.path.exists(dir_idx):
        os.makedirs(dir_idx)
    f_hypo = open(dir_idx+'/'+sfile_idx+'.hypo.bed','w')
    f_hyper = open(dir_idx+'/'+sfile_idx+'.hyper.bed','w')
    if "_DML.sample_level" in f:
        with open(f) as fp:
            for line in fp:
                line_temp = line.strip().split('\t')
                if float(line_temp[2]) > 0:
                    f_hyper.write(line_temp[0]+'\t'+str(int(line_temp[1])-1)+'\t'+line_temp[1]+'\t'+'\t'.join(line_temp[2:])+'\n')
                elif float(line_temp[2]) < 0:
                    f_hypo.write(line_temp[0]+'\t'+str(int(line_temp[1])-1)+'\t'+line_temp[1]+'\t'+'\t'.join(line_temp[2:])+'\n')
    elif "_DMR.sample_level" in f:
        with open(f) as fp:
            for line in fp:
                line_temp = line.strip().split('\t')
                if float(line_temp[-1]) > 0:
                    f_hyper.write(line)
                elif float(line_temp[-1]) < 0:
                    f_hypo.write(line)
    f_hyper.close()
    f_hypo.close()


'''
DSS_results_excl_inflammation/sample_level_filtered/TAL.age_DMR.sample_level.filtered_by_DML.txt
chr22   46367527        46368040        514     25      100.076490973148
chr4    186204222       186204477       256     36      -98.7722564348173
chr19   22633686        22633948        263     19      -94.1310461297973
chr6    291982  292251  270     44      -89.7081367270015


DML

chr1    851548  5.29907579011332        1.16390317214197e-07    0.00180088321554619
chr1    930542  4.8364368603112 1.32187235095477e-06    0.0102381426880486
chr1    932038  -5.62824408411107       1.82053352830246e-08    0.000440652031316814
chr1    1099091 4.37146596282703        1.23415091856118e-05    0.0444290451524966')
'''
