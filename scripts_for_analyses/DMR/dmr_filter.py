import sys
DML = sys.argv[1]
DMR = sys.argv[2]
DML_corrected = sys.argv[3]
DMR_corrected = sys.argv[4]

fout1 = open(DML_corrected, 'w')
fout2 = open(DMR_corrected, 'w')
with open(DML, 'r') as dml:
    dml.readline()
    for line in dml:
        line_temp = line.strip().split('\t')
        fout1.write(line_temp[0]+'\t'+str(int(line_temp[1]))+'\t'+'\t'.join(line_temp[2:])+'\n')

with open(DMR, 'r') as dmr:
    dmr.readline()
    for line in dmr:
        line_temp = line.strip().split('\t')
        fout2.write(line_temp[0]+'\t'+str(int(float(line_temp[1])))+'\t'+str(int(float(line_temp[2])))+'\t'+'\t'.join(line_temp[3:])+'\n')
fout1.close()
fout2.close()



    
