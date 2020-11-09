
f_rpkm = open("rpkm.txt","w")
f_count = open("counts.txt","w")

with open("pancreas_refseq_rpkms_counts.txt", "r") as f:
	sample_line = f.readline()
	sample_ls = sample_line.strip().split("\t")
	sample_ls.remove("NM_ID")
	f_rpkm.write("\t".join(sample_ls) + "\n")
	f_count.write("\t".join(sample_ls) + "\n")
	for line in f:
		line_ls = line.strip().split("\t")
		rpkm = line_ls[2:3516]
		count = line_ls[3516:]

		f_rpkm.write(line_ls[0] + "\t"+"\t".join(rpkm)+"\n")
		f_count.write(line_ls[0] + "\t"+"\t".join(count)+"\n")


f_rpkm.close()
f_count.close()