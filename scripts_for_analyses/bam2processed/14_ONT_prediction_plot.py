import sys
import os,glob
'''
more deconv_LOO_major_celltype.prediction.csv 
Sample  total_num_cells CellType        num_cells       sciMET_celltype_ratio   predicted_celltype_ratio
SDKZ0053        1357    B       41      0.03021370670596905     0.0289421
SDKZ0053        1357    CNT     74      0.054532056005895356    0.014356
SDKZ0053        1357    DCT     49      0.03610906411201179     0.032987
SDKZ0053        1357    EC      27      0.01989683124539425     0.0195859
'''
dd = {}
with open("deconv_LOO_major_celltype.prediction.hypo_only.csv") as f:
    title = f.readline().strip().split("\t")
    for line in f:
        line = line.strip().split("\t")
        sample = line[0]
        celltype = line[2]
        predicted_ratio = float(line[5])
        true_ratio = float(line[4])
        dd.setdefault(sample, {})
        dd[sample][celltype] = line
#"./deconv_LOO_ONT/"
flist = glob.glob("deconv_LOO_ONT/*.deconv.csv")
dd_ont = {}
for f in flist:
    sample = os.path.basename(f).split(".")[0].split("_ONT")[0]
    if "hKid" in sample:
        sample_base = sample.split("_")[1]
    else:
        sample_base = sample
    print (sample_base)
    dd_ont.setdefault(sample_base, {})

    with open(f) as f:
        f.readline()
        for line in f:
            line_temp = line.strip().split(",")
            celltype = line_temp[0]
            ratio = line_temp[1]
            dd_ont[sample_base][celltype] = ratio
fout = open("deconv_LOO_major_celltype.prediction.ONT.hypo_only.csv", "w")
fout.write("\t".join(title) + "\tONT_prediction\n")
for sample in dd_ont:
    if "IMR90" in sample:
        continue
    for celltype in dd_ont[sample]:
        fout.write("\t".join(dd[sample][celltype]) + "\t" + dd_ont[sample][celltype] + "\n")
fout.close()

fout = open("IMR90_deconv_LOO_major_celltype.prediction.ONT.summary.hypo_only.csv", "w")
fout.write("Sample\tCellType\tONTPredictedRatio\n")
for sample in dd_ont:
    if "IMR90" not in sample:
        continue
    for celltype in dd_ont[sample]:
        fout.write(sample + "\t" + celltype + "\t" + dd_ont[sample][celltype] + "\n")
fout.close()




'''
head deconv_LOO_ONT/hKid_SDKZ0053_ONT.deconv.csv 
CellType,hKid_SDKZ0053_ONT
B,0.0298512
CNT,0.0926840
DCT,0.0009256
EC,0.0698575
'''






