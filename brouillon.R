
######################

fastq_file = read.table(paste0(root_dir,"/","2022_GEN_OPDM_pilote09_BC2","/","guppy","/","NOTCH2NLC","/","overlap_NOTCH2NLC_all_100.fastq"),sep=" ",quote="",fill = T)
desc_data=fastq_file[seq(1,dim(fastq_file)[1],4),]
seq_data=fastq_file[seq(2,dim(fastq_file)[1],4),1,drop=F]

keep_rows = which((desc_data$V7>20) & (desc_data$V27>90) & (desc_data$V28>90))
desc_data = cbind(setNames(desc_data[keep_rows,c("V1","V7")],c("read.name","qual")),setNames(seq_data[keep_rows,, drop=F],"seq"))

desc_data$rev_seq = lapply(desc_data[,"seq"], reverse_complement)
desc_data$rev_seq = unlist(lapply(desc_data[,"seq"], reverse_complement))

# find first.codon closest to length of primer
desc_data$start=unlist(lapply(desc_data$rev_seq, function(x,len){pos = gregexpr("CCC",x)[[1]]; pos[which.min(abs(pos-len))]},len=100))
# find las repeat codon close to length of ending primer
desc_data$stop=unlist(gregexpr("GGC",desc_data$rev_seq[1])[[1]][which.min(abs(gregexpr("GGC",desc_data$rev_seq[1])[[1]] - (nchar(desc_data$rev_seq[1])-100)))]+2)
# now extract from start to stop
desc_data$excerpt = unlist(apply(desc_data[,c("rev_seq","start","stop")],1,function(x) {substr(x[1],start=x[2],stop=x[3])}))

# find the CCC that is 
