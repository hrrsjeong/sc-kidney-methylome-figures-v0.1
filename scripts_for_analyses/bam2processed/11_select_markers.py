import sys,glob,os 
#LOO_major_celltype TAL__SDKZ0030__LOO.beta
input_base = "LOO_major_celltype"
#input_base = "LOO_celltype"
input_idx = f"{input_base}_beta"
marker_input_idx = "markers_" + input_idx
flist = glob.glob(f'{marker_input_idx}/*')#/*/*Markers.POD.bed')
dd = {}
for f in flist:
    LOO_sample = f.split('/')[1]
    dd.setdefault(LOO_sample,{})
    marker_files = [x for x in glob.glob(f+'/*/Markers.*.bed') if '.nan.' not in x]
    for marker_file in marker_files:
        cpg_cov_type = marker_file.split('/')[-2]
        with open(marker_file) as fin:
            title = fin.readline()
            for line in fin:
                if line.startswith('#'):continue
                line_temp = line.strip().split('\t')
                idx = line_temp[6]
                celltype = line_temp[5]
                direction = line_temp[-1]
                tg_mean = float(line_temp[9])
                bg_mean = float(line_temp[10])
                delta_means = float(line_temp[11])
                pval = float(line_temp[-2])
                if int(line_temp[7].split("CpG")[0]) < 5:continue
                flag = 1
                if pval > 0.000001:continue
                if direction == 'U':
                    if tg_mean > 0.25:
                        flag = 0
                    if delta_means < 0.4:
                        flag = 0
                    if bg_mean < 0.7:
                        flag = 0
                if direction == 'M':
                    flag = 0
                    continue
                    if delta_means < 0.4:
                        flag = 0
                    if tg_mean < 0.7:
                        flag = 0
                    if bg_mean > 0.2:
                        flag = 0
                if flag == 0:continue
                dd[LOO_sample].setdefault(celltype,{})
                dd[LOO_sample][celltype][idx] = (pval,line_temp)
os.system(f'mkdir -p sorted_{marker_input_idx}')
os.system('mkdir -p tmp')
os.system(f"mkdir -p Atlas_{marker_input_idx}/tmp")
os.system(f"mkdir -p deconv_{input_base}")
num_top = 250
for LOO_sample in dd:
    print(LOO_sample)
    os.system(f'mkdir -p sorted_{marker_input_idx}/'+LOO_sample)
    for celltype in dd[LOO_sample]:
        #sort by pval
        sorted_idx = sorted(dd[LOO_sample][celltype],key=lambda x:dd[LOO_sample][celltype][x][0])
        with open(f'sorted_{marker_input_idx}/'+LOO_sample+'/'+celltype+'.bed','w') as fout:
            fout.write(title)
            for idx in sorted_idx:
                fout.write('\t'.join(dd[LOO_sample][celltype][idx][1])+'\n')
    with open(f"tmp/{LOO_sample}_top{num_top}.sh",'w') as fout:
        systemstr = f'cat <(head -n 1 sorted_{marker_input_idx}/'+LOO_sample+'/'+celltype+'.bed) '
        systemstr += ' '.join([f'<(head -n {num_top} sorted_{marker_input_idx}/{LOO_sample}/{celltype}.bed | grep -v "#")' for celltype in dd[LOO_sample]])
        systemstr += f' >| sorted_{marker_input_idx}/'+LOO_sample+f'/top{num_top}.bed'
        fout.write(systemstr)
    os.system(f"bash tmp/{LOO_sample}_top{num_top}.sh")

    systemstr = f"~/Programs/UXM_deconv/uxm build --markers sorted_{marker_input_idx}/{LOO_sample}/top{num_top}.bed -o Atlas_{marker_input_idx}/tmp/{LOO_sample}.top{num_top}.tsv --pats {input_base}_pat/*__{LOO_sample}__LOO.pat.gz --tmp_dir sorted_{marker_input_idx}/{LOO_sample}/ -f --threads 32 -l 4"
    os.system(systemstr)
    with open(f"Atlas_{marker_input_idx}/tmp/{LOO_sample}.top{num_top}.tsv",'r') as fin:
        atlas_title = fin.readline().strip().split('\t')
        new_atlas_title = atlas_title[:8] + [x.split('__')[0] for x in atlas_title[8:]]
        with open(f"Atlas_{marker_input_idx}/{LOO_sample}.top{num_top}.tsv",'w') as fout:
            fout.write('\t'.join(new_atlas_title)+'\n')
            for line in fin:
                line_temp = line.strip().split('\t')
                if "\tNA\t" in line:continue
                fout.write(line)
    systemstr = f"~/Programs/UXM_deconv/uxm deconv --atlas Atlas_{marker_input_idx}/{LOO_sample}.top{num_top}.tsv --output deconv_{input_base}/{LOO_sample}.deconv.csv --tmp_dir deconv/{LOO_sample}/ -v -f --threads 16 -l 4 sample_pat/{LOO_sample}.pat.gz"
    os.system(systemstr)
    





'''
markers_LOO_major_celltype_beta/SDKZ0030/cpg15_cov45/Markers.POD.bed
#chr    start   end     startCpG        endCpG  target  region  lenCpG  bp      tg_mean bg_mean delta_means     delta_quants    delta_maxmin    ttest   direction
chr1    1548529 1548903 36202   36215   POD     chr1:1548529-1548903    13CpGs  374bp   0.209   0.866   0.657   0.429   0.384   3.3e-10 U
chr1    9319293 9319731 195089  195107  POD     chr1:9319293-9319731    18CpGs  438bp   0.293   0.902   0.609   0.52    0.505   3.61e-14        U

'''
