

plot.methyl.fancy <- function(DIR,CALLER,gene,methyl, root="C:/Users/boelle/nextcloud_SU/paper/",
                              from=NULL, to=NULL, type="s", reads=c("short","long")) {
  #methyl is type of methylation, found in files
  # type.methyl 1 or 2 depends on what will be read in modkit
  if (type==1) {type.methyl=1}
  else if (type==2) {type.methyl=2}
  else if (type=="s") {type.methyl=3}
  else if (type=="b") {type.methyl=4}
  else {stop("type must be 1,2, (s)um or (b)oth")}
  
  options(ucscChromosomeNames=FALSE) # necessary when working with txdb as chromosome are not chr natively
  library(Gviz)
  library(rtracklayer)
  
  gtrack <- GenomeAxisTrack()
  # get coordinates of gene
  library(GenomicFeatures)
  
  GENE = genes(txdb, filter=list(gene_id=gene), columns=c("gene_id","tx_chrom"))
  if (length(GENE) != 1) stop("check GENE step")
  start.gene = start(GENE)
  end.gene = end(GENE)
  
  if (is.null(from)) from=start.gene
  if (is.null(to)) to=end.gene
  
  grtrack <- GeneRegionTrack(txdb,
                             chromosome = as.character(chrom(GENE)), 
                             start = from, end = to, 
                             stacking="dense",name="gene")
  
  setwd(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
  library(bsseq)
  # verifier que les fichiers existent
  for (xx in reads) {
    if (!file.exists(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",xx,"_",gene,"_",methyl,".bedmethyl")))
      stop(paste0("file ",xx,"_",gene,"_",methyl,".bedmethyl not found"))
  }
  # quoi faire sur les methylation
  bs.methyl <- read.modkit(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",reads,"_",gene,"_",methyl,".bedmethyl"), rmZeroCov = TRUE)
  
  if (type.methyl==1) {
    bs.methyl=list(bs.methyl[[1]])
    message("using 1st methylation")
  }else if (type.methyl==2) {
    if (length(bs.methyl)>1) {
      bs.methyl=list(bs.methyl[[2]])
      message("using 2nd methylation")
    } else {
      warning("requested 2nd methylation but only 1 is available - defaulting to 1st")
      bs.methyl=list(bs.methyl[[1]])
    }
  } else if (type.methyl==3) {
    if (length(bs.methyl)>1) {
      #add M of 2 to 1
      bs.methyl[[1]]@assays@data@listData$M = bs.methyl[[1]]@assays@data@listData$M +bs.methyl[[2]]@assays@data@listData$M 
      message("added 1st and 2nd methylation together")
      bs.methyl=list(bs.methyl[[1]])
    } else {
      warning("requested add methylation but only 1 methylation available - defaulting to 1st")
    }
  } else if(type.methyl=="b") {
    if (length(bs.methyl)>1) {
      #keep the 2 methylations
      bs.methyl=list(bs.methyl[[1]], bs.methyl[[2]])
    } else {
      warning("requested two methylations but only 1 methylation available - defaulting to 1st")
    }
  }
  
  mat.methyl <- NULL
  mat.cov <- NULL
  names.methyl=NULL
  methyl.sm<-list()
  for (i in 1:length(bs.methyl)) {
    # filter based on coverage - use methyl1 so that get the same ranges
    idx.cov = which(apply(getCoverage(bs.methyl[[1]]),1,sum)>ceiling((max(apply(getCoverage(bs.methyl[[1]]),1,sum))/10)))
    methyl.sm[[i]] <- BSmooth(bs.methyl[[i]][idx.cov,])
    seqlevels(methyl.sm[[i]]@rowRanges) <- as.character(chrom(GENE))
    mat.methyl = rbind( mat.methyl,
                        t(getMeth(methyl.sm[[i]], type="smooth"))  )
    names.methyl = c(names.methyl, paste0(reads,".",i))
    mat.cov = rbind(mat.cov, t(getCoverage(methyl.sm[[i]], type="Cov")))
  }
  # depending on 
  rownames(mat.methyl)=names.methyl
  rownames(mat.cov)=names.methyl
  
  # restrict plot to coverage of reads if set to gene length
  if ((from==start.gene) & (from < min(start(methyl.sm[[1]]@rowRanges)))) { from=min(start(methyl.sm[[1]]@rowRanges)) }
  if ((to==end.gene) & (max(end(methyl.sm[[1]]@rowRanges))<to)) {to = max(end(methyl.sm[[1]]@rowRanges)) }
  
  ctrack <- DataTrack(
    range = methyl.sm[[1]]@rowRanges,
    data = mat.cov ,
    type = "l",
    name = "coverage",
    chromosome = as.character(chrom(GENE)), 
    genome="T2T-CHM13v2.0"
  )
  
  mtrack <- DataTrack(
    range = methyl.sm[[1]]@rowRanges,
    data =  mat.methyl,
    type = "l",
    name = "methylation",
    chromosome = as.character(chrom(GENE)), 
    genome="T2T-CHM13v2.0"
  )
  
  repeat.coords = repeat.coords[repeat.coords$chrom == paste0("chr",as.character(chrom(GENE))) & repeat.coords$start > start.gene & repeat.coords$end < end.gene, ]
  
  rtrack = AnnotationTrack(
    start=repeat.coords$start, width=repeat.coords$end-repeat.coords$start,
    chromosome= as.character(chrom(GENE)),id=repeat.coords$rep,
    groupAnnotation="id", 
    name="repeats" )
  
  pdf(paste0(root,"/",DIR,"/",CALLER,"/methylation_",gene,"_",methyl,"_",type,".pdf"))
  plotTracks(list(gtrack,ctrack,mtrack,grtrack,rtrack), 
             groups=names.methyl,  
             from = from, to = to,
             fontcolor.feature=1, cex.feature=0.5)
  dev.off()  
  
}



txdb = loadDb("C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0") # -> txdb
load("C:/Users/boelle/nextcloud_SU/paper/repeats.rda") # repeat.coords


fix.collapse <- function (DIR,CALLER, gene,methyl,methyl_mod,root_dir) {
  # make sure that bedmethyl files are first combined in strands
  # this is needed when using unfiltered as these are present in +/-
  for (x in c("short","long")) {
    
    data=read.table(file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                x,"_",gene,"_",methyl,".bedmethyl"))
    data2 = aggregate(data[,c("V5",paste0("V",10:18))], 
                      by=list(V1=data$V1,V2=data$V2,V3=data$V3,V4=data$V4,V7=data$V7,V8=data$V8,V9=data$V9), sum)
    data2$V6="-"
    data2=data2[,paste0("V",1:18)]
    write.table(data2,file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                  x,"_",gene,"_",methyl_mod,".bedmethyl"), row.names=F,col.names=F, sep="\t",quote = FALSE)
  }
  
}

# in bedmethyl
# V10 : coverage
# V12 : mod (-> m)
# V13 : canonical
# V13 : other_mod (-> h)


split.and.collapse <- function (DIR,CALLER, gene,methyl,context,root_dir) {
  # split the bedmethyl files in different files according to h/m/a
  # collapse strands at same positions +/-
  # METHYL IS THE TYPE OF METHYL CALLED found in the name (5mC, 5mCG, 5mC_5hmC,5mCG_5hmCG,)
  # context is unfiltered / CPGfiltered
  # files are changed to m/h/a + context
  for (x in c("short","long","overlap")) {
    data=read.table(file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                x,"_",gene,"_",methyl,"_",context,".bedmethyl"))
    for (mod in unique(data$V4)) {
      # this should be h/m/a
      data2 = data[data$V4 == mod,]
      # aggregate over strands 
      data2 = aggregate(data2[,c("V5",paste0("V",10:18))], 
                        by=list(V1=data2$V1,V2=data2$V2,V3=data2$V3,V4=data2$V4,V7=data2$V7,V8=data2$V8,V9=data2$V9), sum)
      
      data2$V6="-"
      # put canonical to canonical + mod_other
      data2$V13=data2$V13 + data2$V14
      # put mod_other to 0
      data2$V14=0
      # valid_coverage
      data$V10 = data$V12 + data$V13
      # remove coverage=0
      data2=data2[data2$V10>0,]
      #reorder
      data2=data2[,paste0("V",1:18)]
      write.table(data2,file=paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/","METHYLATION","/",
                                    x,"_",gene,"_",mod,"_",context,".bedmethyl"), row.names=F,col.names=F, sep="\t",quote = FALSE)
    }
  }
}


check.short.long <- function(DIR,CALLER, root_dir) {
  # find the "ids" files in directories containing read-id and repeat-length
  # make a table showing how many reads in each and median length
  setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/"))
  #  get gene dirs
  genes = list.dirs(recursive = F)
  genes = genes [genes != "./OVERLAPS"]
  res=c()
  for (gene in genes) {
    setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
    liste.files = dir(pattern="(long|short|overlap)[a-zA-Z0-9_]*.ids")
    for (file in liste.files) {
      if (file.size(file)>0) {
        tab <- read.table(file)
        nb.reads = dim(tab)[1]
        len.repeat=median(tab[,2])
      } else {
        nb.reads=0
        len.repeat=NA
      }
      
      res <- rbind(res,data.frame(gene=gene, nb.reads=nb.reads, len.repeat=len.repeat,file=file,dir=DIR))
    }
    setwd(paste0(root_dir,"/",DIR,"/",CALLER,"/"))
  }
  setwd(root_dir)
  res
}


