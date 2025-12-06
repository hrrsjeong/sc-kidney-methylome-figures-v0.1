import sys,glob,os

flist = glob.glob("bismark_context_file_by_sample_fofn/*__*.fofn")#POD__SDKZ0057.fofn")
dd_true_cnt_celltype = {}
dd_true_cnt_major_celltype = {}
for sfile in flist:
    celltype = sfile.split("/")[-1].split("__")[0]
    sample = sfile.split("/")[-1].split("__")[1].split(".")[0]
    if "-" in celltype:
        major_celltype = celltype.split('-')[0]
    elif celltype == "dPC":
        major_celltype = "PC"
    else:
        major_celltype = celltype
    cnt = 0
    with open(sfile,'r') as fp:
        for line in fp:
            cnt +=1
    dd_true_cnt_celltype.setdefault(sample,{})
    dd_true_cnt_celltype[sample][celltype] = cnt
    dd_true_cnt_major_celltype.setdefault(sample,{})
    dd_true_cnt_major_celltype[sample].setdefault(major_celltype,0)
    dd_true_cnt_major_celltype[sample][major_celltype] += cnt

dd_sample_cnt = {}
for sample in dd_true_cnt_celltype:
    total_cnt = 0
    for celltype in dd_true_cnt_celltype[sample]:
        total_cnt += dd_true_cnt_celltype[sample][celltype]
    dd_sample_cnt[sample] = total_cnt

for stype in ["deconv_LOO_celltype","deconv_LOO_major_celltype"]:
    print (stype)
    dd_deconv = {}
    flist = glob.glob(f"{stype}/*.deconv.csv")
    for sfile in flist:
        sample = sfile.split("/")[-1].split(".")[0]
        with open(sfile,'r') as fp:
            fp.readline()
            for line in fp:
                line_temp = line.strip().split(",")
                celltype = line_temp[0]
                celltype_ratio = float(line_temp[1])
                dd_deconv.setdefault(sample,{})
                dd_deconv[sample][celltype] = celltype_ratio
    if stype == "deconv_LOO_celltype":
        dd_true_cnt = dd_true_cnt_celltype
    else:
        dd_true_cnt = dd_true_cnt_major_celltype
    with open(f"{stype}.prediction.hypo_only.csv",'w') as fp:
        fp.write("Sample\ttotal_num_cells\tCellType\tnum_cells\tsciMET_celltype_ratio\tpredicted_celltype_ratio\n")
        for sample in dd_deconv:
            print (sample)
            for celltype in dd_deconv[sample]:
                if celltype not in dd_true_cnt[sample]:
                    num_cells = 0
                else:
                    num_cells = dd_true_cnt[sample][celltype]
                total_cnt = dd_sample_cnt[sample]
                predicted_celltype_ratio = dd_deconv[sample][celltype]
                sciMET_celltype_ratio = num_cells/total_cnt
                fp.write(f"{sample}\t{total_cnt}\t{celltype}\t{num_cells}\t{sciMET_celltype_ratio}\t{predicted_celltype_ratio}\n")

'''
CellType,SDKZ0053
B,0.0289421
CNT,0.0143560
DCT,0.0329870
'''




