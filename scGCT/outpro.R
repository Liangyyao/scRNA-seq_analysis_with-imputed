## last update: 2020/03/28

args <- commandArgs(T)
job_name = args[1]
tmp = args[2]
outdir = args[3]
timestamp = args[4]
library(data.table)
load(paste(tmp,"/original.RData", sep = ""))
file = paste(tmp,"/scIGANs-", job_name,".csv", sep = "")

file = "D:/bigdata/新加实验5.20/imputed/Timecourse_scGCT.csv"
d <- fread(file, header = T)
gcm <- d[,-1]*reads_max_cell


gcm_out <-  t(gcm[,1:geneCount])
colnames(gcm_out) = cellnames
#gcm_out <- cbind(Gene_ID = genenames,gcm_out)
rownames(gcm_out) = genenames

write.csv(gcm_out,file = 'D:/bigdata/新加实验5.20/imputed/Timecourse_scGCT.csv')

outfile = paste(outdir,"/scIGANs_",job_name,".txt", sep = "")
fwrite(as.data.table(gcm_out), file=outfile, col.names = T, row.names = F, sep = "\t", quote = F)
message(paste("\nCompleted!!!", Sys.time()))
message(paste("Imputed matrix:", outfile))

