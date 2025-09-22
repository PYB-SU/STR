

get.gene <- function(gene) {
  GENE = GENES[GENES$gene_id==gene,]
  if (length(GENE) != 1) {
    if (length(GENE)==0) {
      stop(paste("gene",gene,"not found in TxDb file"))
    } else {
      print(GENE)
      stop(paste("gene",gene,"is not unique in TxDb file"))
    }
  }

  return(GENE)
}

get.repeat.coords <- function(GENE) {
  options(ucscChromosomeNames=FALSE) # necessary when working with txdb as chromosome are not chr natively
  repeat.coords.gene = REPEATS[REPEATS$chrom == as.vector(rtracklayer::chrom(GENE)) & REPEATS$start > BiocGenerics::start(GENE) & REPEATS$end < BiocGenerics::end(GENE), ]
  return(repeat.coords.gene)
}

process.methyl.files <- function(DIR, CALLER, gene, methyl, root, repeat.rows=NULL, reads) {
  #compute methylation percentages from raw files
  # this must be streamlined as files are only for one type of methylation
  options(ucscChromosomeNames=FALSE) # necessary when working with txdb as chromosome are not chr natively

  GENE <- get.gene(gene)

  setwd(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
  library(bsseq)
  # verifier que les fichiers existent
  for (xx in reads) {
    if (!file.exists(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",xx,"_",gene,"_",methyl,".bedmethyl")))
      stop(paste0("file ",xx,"_",gene,"_",methyl,".bedmethyl not found in ",root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))
  }
  # quoi faire sur les methylation
  bs.methyl <- bsseq::read.modkit(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/",reads,"_",gene,"_",methyl,".bedmethyl"), rmZeroCov = TRUE)

  mat.methyl.sm <- NULL
  mat.methyl.rw <- NULL
  mat.cov <- NULL
  names.methyl=NULL
  names.cov=NULL

  methyl.sm<-list()
  for (i in 1:length(reads)) {
    # filter based on coverage - use methyl1 so that get the same ranges
    idx.cov = which(apply(bsseq::getCoverage(bs.methyl[,1]),1,sum)>min(10,ceiling((max(apply(bsseq::getCoverage(bs.methyl[,1]),1,sum))/10))))
    methyl.sm[[i]] <- (bsseq::BSmooth(bs.methyl[,i], ns=min(length(idx.cov)/2,70)))[idx.cov,]
    seqlevels(methyl.sm[[i]]@rowRanges) <- as.character(chrom(GENE))
    mat.methyl.sm = rbind( mat.methyl.sm,
                        t(bsseq::getMeth(methyl.sm[[i]], type="smooth"))  )
    mat.methyl.rw = rbind( mat.methyl.rw,
                        t(bsseq::getMeth(methyl.sm[[i]], type="raw"))  )
    names.methyl = c(names.methyl, reads[i])
    mat.cov = rbind(mat.cov, t(bsseq::getCoverage(methyl.sm[[i]], type="Cov")))
    names.cov = c(names.cov, reads[i])
  }
  # get ranges
  row.ranges = methyl.sm[[1]]@rowRanges
  names(row.ranges)=paste0(methyl.sm[[1]]@rowRanges@seqnames,":",methyl.sm[[1]]@rowRanges@ranges@start)
  #d'abord smooth, apres raw
  rownames(mat.methyl.rw)=names.methyl
  rownames(mat.methyl.sm)=names.methyl
  #  rownames(mat.methyl)[rownames(mat.methyl)=="rw.short.1"]="non-expanded"
  #  rownames(mat.methyl)[rownames(mat.methyl)=="rw.long.1"]="expanded"

  if("overlap" %in% reads) {
    names.cov[names.cov=="overlap"]="all"
    names.methyl[names.methyl=="overlap"]="all"
  }
  if("short" %in% reads) {
    names.cov[names.cov=="short"]="non-expand."
    names.methyl[names.methyl=="short"]="non-expand."
  }
  if("long" %in% reads) {
    names.cov[names.cov=="long"]="expanded"
    names.methyl[names.methyl=="long"]="expanded"
  }
  if("H1" %in% reads) {
    names.cov[names.cov=="H1"]="haplo 1"
    names.methyl[names.methyl=="rw.H1"]="haplo 1"
  }
  if("H2" %in% reads) {
    names.cov[names.cov=="H2"]="haplo 2"
    names.methyl[names.methyl=="H2"]="haplo 2"
  }
  rownames(mat.cov)=names.cov
  rownames(mat.methyl.sm)=names.methyl
  rownames(mat.methyl.rw)=names.methyl

  if (is.null(repeat.rows)) {
    repeat.coords.gene = REPEATS[(REPEATS$chrom == rtracklayer::chrom(GENE)) & (REPEATS$start > BiocGenerics::start(GENE)) & (REPEATS$end < BiocGenerics::end(GENE)), ]
  } else {
    # if we want to choose directly which rows to plot
    repeat.coords.gene=REPEATS[repeat.rows,,drop=F]
  }
  return(list(repeat.coords.gene=repeat.coords.gene, mat.cov=mat.cov, mat.methyl.sm=mat.methyl.sm,mat.methyl.rw=mat.methyl.rw,GENE=GENE,row.ranges=row.ranges))
}


plot.methyl <- function(obj, DIR,CALLER,gene,methyl, root,
                              from=NULL, to=NULL,
                              coverage=FALSE,focus=NULL,
                        save="yes", file_name=NULL) {
  # methyl is type of methylation, found in files
  # smooth =TRUE si on veut lisser le % de methylation, sinon smooth=FALSE
  # plot names for substitute of read categoriess
  options(ucscChromosomeNames=FALSE) # necessary when working with txdb as chromosome are not chr natively
#  library(Gviz)
#  library(rtracklayer)
  # to store pdf

  setwd(paste0(root,"/",DIR,"/",CALLER,"/",gene,"/METHYLATION/"))

  gtrack <- GenomeAxisTrack()
  # get coordinates of gene
#  library(GenomicFeatures)
  list.to.plot <- list(gtrack)
  sizes.to.plot <- c(1)

  start.gene = BiocGenerics::start(obj$GENE)
  end.gene = BiocGenerics::end(obj$GENE)
  chromosome = as.vector(rtracklayer::chrom(obj$GENE))
  genome = GenomeInfoDb::genome(obj$GENE)[chromosome]

  if (is.null(from)) from=start.gene
  if (is.null(to)) to=end.gene

  # restrict plot to coverage of reads if set to gene length
  if ((from==start.gene) & (from < min(BiocGenerics::start(obj$row.ranges)))) { from=min(BiocGenerics::start(obj$row.ranges)) }
  if ((to==end.gene) & (max(end(BiocGenerics::end(obj$row.ranges))) < to)) {to = max(BiocGenerics::end(obj$row.ranges)) }


  if (coverage==TRUE) {
    ctrack <- DataTrack(
      range = obj$row.ranges,
      data = obj$mat.cov ,
      type = rep("l",dim(obj$mat.cov)[1]),
      name = "coverage",
      chromosome = chromosome,
      genome=genome,
      groups=rownames(obj$mat.cov)
    )
    list.to.plot <-append(list.to.plot,list(ctrack))
    sizes.to.plot <- c(sizes.to.plot,2)
  } else {
    ctrack <- NULL
  }

  mtrack.sm <- DataTrack(
      range = obj$row.ranges,
      data =  obj$mat.methyl.sm,
      type = rep("l",dim(obj$mat.methyl.sm)[1]),
      name = "methylation",
      chromosome = chromosome,
      genome=genome,
      groups=rownames(obj$mat.methyl.sm),
      ylim = c(0,1)
    )

  mtrack.rw <- DataTrack(
    range = obj$row.ranges,
    data =  obj$mat.methyl.rw,
    type = rep("p",dim(obj$mat.methyl.rw)[1]),
    cex=0.6,
    name = "methylation",
    chromosome = chromosome,
    genome=genome,
    groups=rownames(obj$mat.methyl.rw),
    ylim = c(0,1),
    fontcolor.legend="#FFFFFF"
  )

  mtrack=OverlayTrack(list(mtrack.rw,mtrack.sm), name="methylation")

  list.to.plot <-append(list.to.plot,list(mtrack))
  sizes.to.plot <- c(sizes.to.plot,4)

  for (x in 1:dim(obj$mat.methyl.rw)[1]) {
    mtrack.tmp <- DataTrack(
      range = obj$row.ranges,
      data =  obj$mat.methyl.rw[x,,drop=F],
      type = "gradient",
      gradient=colorRampPalette(c("blue", "white", "red"))(100),
      name = rownames(obj$mat.methyl.rw)[x],
      chromosome = chromosome,
      genome=genome,
      ylim = c(0,1)
    )
    list.to.plot <-append(list.to.plot,list(mtrack.tmp))
    sizes.to.plot <- c(sizes.to.plot,1)
  }

#  grtrack <- GeneRegionTrack(txdb,
#                             chromosome = chromosome,
#                             start = from, end = to,
#                             name="gene", # stacking="dense"
#                             collapseTranscripts = "meta", shape = "arrow",
#                             transcriptAnnotation = "gene")
#
#  list.to.plot <-append(list.to.plot,list(grtrack))
#  sizes.to.plot <- c(sizes.to.plot,2)


  rtrack = AnnotationTrack(
    start=obj$repeat.coords.gene$start, width=obj$repeat.coords.gene$end-obj$repeat.coords.gene$start,
    chromosome= chromosome,
    id=obj$repeat.coords$rep,
    groupAnnotation="id",
    name="repeats" )

  list.to.plot <-append(list.to.plot,list(rtrack))
  sizes.to.plot <- c(sizes.to.plot,1)

  if (!is.null(focus)) {
    from=obj$repeat.coords$start-focus
    to=obj$repeat.coords$end+focus
    if (save=="yes") {
      pdf(paste0(root,"/",DIR,"/",CALLER,"/methylation_",gene,"_",methyl,"_focus.pdf"))
    }
  } else {
    if (save=="yes") {
      if (is.null(file_name)) {
        file_name=paste0(root,"/",DIR,"/",CALLER,"/methylation_",gene,"_",methyl,".pdf")
      } else {
        file_name=paste0(root,"/",DIR,"/",CALLER,"/",file_name)
      }
      pdf(file_name)
    }
  }

  res=plotTracks(list.to.plot,
                 sizes=sizes.to.plot,
                 from = from, to = to,
                 fontcolor.feature=1, cex.feature=0.5)

  if (save=="yes") {
    dev.off()
  }
  setwd(root)
  message(paste("setwd to ",root))
  return(res)
}


fix.collapse <- function (DIR,CALLER, gene,methyl,methyl_mod,root_dir,reads=c("short","long")) {
  # make sure that bedmethyl files are first combined in strands
  # this is needed when using unfiltered as these are present in +/-
  for (x in reads) {

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


split.and.collapse <- function (DIR,CALLER, gene,methyl,context,root_dir, alleles=c("short","long","overall")) {
  # split the bedmethyl files in different files according to h/m/a
  # collapse strands at same positions +/-
  # METHYL IS THE TYPE OF METHYL CALLED found in the name (5mC, 5mCG, 5mC_5hmC,5mCG_5hmCG,)
  # context is unfiltered / CPGfiltered
  # files are changed to m/h/a + context
  for (x in alleles) {
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


