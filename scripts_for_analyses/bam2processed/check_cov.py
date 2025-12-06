import sys,os
total_cpg = 29345332
finp = sys.argv[1]
fout3 = sys.argv[2]
total_cov = 0
sample_id = os.path.basename(finp).split('.')[0]
with open(finp) as f:
    for line in f:
        line_temp = line.strip().split('\t')
        cov = float(line_temp[-1])
        total_cov += cov
avg_cov = total_cov / total_cpg
with open(fout3, 'w') as f:
    f.write(f'{sample_id}\t{avg_cov}\n')
