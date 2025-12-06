import sys, os
import glob

dd = {}
with open("../../sample_id.txt", "r") as f:
    f.readline()
    for line in f:
        line_temp = line.strip().split('\t')
        State,Age,Gender,Pct_cortex = line_temp[5:]
        Pct_cortex = str(float(Pct_cortex)/100)
        Age = int(Age)
        if Age < 40:
            Age_level = "Young"
        elif Age < 60:
            Age_level = "Mid"
        else:
            Age_level = "Old"
        Age = str(Age)
        if State == "Ref":
            Disease_level = "0"
            state_level = "control"
        elif State == "AKI":
            Disease_level = "1"
            state_level = "disease"
        elif State == "CKD":
            Disease_level = "1"
            state_level = "disease"
        if line_temp[4] in ["SDKZ0048","SDKZ0052"]:
            state_level = "control"
            #continue
        print (line)
        dd[line_temp[4]] = [state_level,Age,Gender,Disease_level,Pct_cortex,Age_level]

fout = open("sample_covariate_altered.txt", "w")
fout.write("Sample\tSample_id\tCelltype\tCelltype_origin\tCell_state\tState\tAge\tGender\tDisease_level\tPct_cortex\tAge_level\n")
flist = glob.glob("./DSS_input_altered/rds/*.rds")
for sfile in flist:
    sname = os.path.basename(sfile).split(".")[0]
    #celltype,sample_id = sname.split("_")
    print (sname)
    sample_id,celltype = sname.split("__")
    if sample_id in dd:
        #print(sname, sample_id, celltype, "\t".join(dd[sample_id]))
        fout.write(sname + "\t" + sample_id+'\t'+ celltype +'\t'+celltype.split('-')[0]+'\t'+celltype.split('-')[-1]+'\t'+"\t".join(dd[sample_id]) + "\n")
fout.close() 


'''
kidney-A        HsKidAt_20230316A       3778    SDKZ0030-KID-2-1        SDKZ0030        Ref     46      Male
kidney-B        HsKidAt_20230316B       3825    SDKZ0057-KID-2-1        SDKZ0057        Ref     78      Female
kidney-C        HsKidAt_20230316C       3782    SDKZ0122-KID-1-1        SDKZ0122        AKI     24      Male
A       HsKidAt_20231114A       3477    SDKZ0032-KID-1-2,3,4    SDKZ0032        CKD     64      Female
B       HsKidAt_20231114B       3479    SDKZ0033-KID-1-2,3      SDKZ0033        CKD     55      Male
C       HsKidAt_20231114C       3487    SDKZ0034-KID-1-2,3,4,5  SDKZ0034        CKD     70      Male
D       HsKidAt_20231114D       3716    SDKZ0036-KID-2-1,2      SDKZ0036        CKD     54      Female
'''
