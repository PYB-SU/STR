
val9 <-
  structure(list(V2 = c(474.02335369999997, 524.29398660000004, 
                        453.57659419999999, 540.04266840000003, 464.85949890000001, 497.413118, 
                        443.26412169999998, 526.80962490000002, 449.05519750000002, 510.26957220000003, 
                        427.23132229999999, 524.78619279999998, 467.28340429999997, 497.13356920000001, 
                        445.08893130000001, 519.85179259999995, 481.14169809999999, 521.89437559999999, 
                        471.74168509999998, 555.84872080000002, 489.7109987, 507.95999790000002, 
                        476.43582320000002, 553.75842460000001, 469.2580562, 519.31304969999997, 
                        464.07931760000002, 551.7577, 488.81197220000001, 503.22960979999999, 
                        476.93733120000002, 530.72615340000004, 464.50567189999998, 503.732732, 
                        445.21861039999999, 527.04972280000004, 461.77124800000001, 489.8775023, 
                        453.69806579999999, 517.76841809999996, 444.3112122, 498.21450220000003, 
                        438.95986900000003, 516.70447479999996, 463.21120610000003, 487.12313410000002, 
                        457.88073739999999, 513.363158, 467.79276420000002, 502.66955469999999, 
                        448.487482, 523.29469129999995, 459.56125850000001, 481.70112399999999, 
                        440.4946516, 505.30680130000002, 440.18604809999999, 487.66964680000001, 
                        429.98995839999998, 503.9096318, 454.3520699, 466.56330209999999, 
                        439.77914559999999, 485.49385539999997)), row.names = c("AAA", 
                                                                                "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", 
                                                                                "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", 
                                                                                "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", 
                                                                                "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", 
                                                                                "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", 
                                                                                "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", 
                                                                                "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
                        ), class = "data.frame")

val10 <-
  structure(list(V2 = c(522.94817079999996, 666.71134470000004, 
                        569.75470529999996, 895.68127990000005, 575.50214649999998, 673.5238008, 
                        608.52667629999996, 895.99962719999996, 620.70456479999996, 724.86084170000004, 
                        679.68797759999995, 955.47837570000002, 602.22853580000003, 642.7570144, 
                        624.41357519999997, 870.45266179999999, 864.98512430000005, 771.60865190000004, 
                        924.19338119999998, 883.46408129999998, 836.91231019999998, 801.57243059999996, 
                        881.53114930000004, 916.7638025, 920.64376649999997, 820.74915580000004, 
                        985.33633750000001, 916.17314299999998, 835.01777579999998, 742.50265869999998, 
                        847.6870844, 844.74583389999998, 618.52544850000004, 685.76105370000005, 
                        627.29857049999998, 978.30413750000002, 639.6601422, 684.05162529999996, 
                        661.05637830000001, 966.94902479999996, 705.55530190000002, 721.92711029999998, 
                        735.12771650000002, 1016.6217370000001, 657.26649129999998, 650.26974589999998, 
                        685.19370479999998, 912.22098419999998, 1040.3038959999999, 965.62539030000005, 
                        1091.7344089999999, 980.25918620000004, 1071.3640829999999, 1020.053996, 
                        1097.900862, 984.39739239999994, 1146.1396030000001, 1059.9624839999999, 
                        1222.978177, 1051.245075, 1009.4064059999999, 971.58933560000003, 
                        1023.5851290000001, 941.91063250000002)), row.names = c("AAA", 
                                                                                "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", 
                                                                                "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "CAA", "CAC", "CAG", 
                                                                                "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", 
                                                                                "CTA", "CTC", "CTG", "CTT", "GAA", "GAC", "GAG", "GAT", "GCA", 
                                                                                "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", 
                                                                                "GTG", "GTT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", 
                                                                                "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT"
                        ), class = "data.frame")



make.group <- function(val) {
  groupe=1
  for (x1 in c("A","C","G","T")) {
    for (x2 in c("A","C","G","T")) {
      for (x3 in c("A","C","G","T")) {
        if (is.na(val[paste0(x1,x2,x3),"groupe"])) {
            val[paste0(x1,x2,x3),"groupe"]=groupe
            val[paste0(x2,x3,x1),"groupe"]=groupe
            val[paste0(x3,x1,x2),"groupe"]=groupe
            groupe=groupe+1
        }
      }
    }
  }
  val
}

val9$groupe=NA
val9 = make.group(val9)
val9$codon = rownames(val9)

val10$groupe=NA
val10 = make.group(val10)
val10$codon = rownames(val10)

val9.ag = aggregate(val9[,"V2", drop=F], by=list(val9$groupe),mean)
val10.ag = aggregate(val10[,"V2", drop=F], by=list(val10$groupe),mean)

val9.lab = aggregate(val9[,"codon", drop=F], by=list(val9$groupe),paste,collapse="/")
val10.lab = aggregate(val10[,"codon", drop=F], by=list(val10$groupe),paste,collapse="/")

val9.small = merge(val9.ag, val9.lab, by="Group.1")

val10.small = merge(val10.ag, val10.lab, by="Group.1")

library(ggplot2)
library(ggdendro)


val9$colour='black'
val9$colour[val9$V2 %in%c("GCC","CGC","CCG")] = "red"

dd.row <- as.dendrogram(hclust(dist(val9)))
ddata_x <- dendro_data(dd.row)

p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

labs <- label(ddata_x)
labs$colour = 'black'
labs$colour[labs$label %in%  c("GCC","CGC","CCG")]  = 'red'
  
p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=0, angle=90, hjust=1.1, colour=labs$colour), size=4)+
  theme_bw()+ylim(-20,200)


val10$colour=
val10$colour[val10$V2 %in%] = "red"


dd.row <- as.dendrogram(hclust(dist(val9)))
ddata_x <- dendro_data(dd.row)

p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

labs <- label(ddata_x)
labs$colour = 'black'
labs$colour[labs$label %in%  c("GCC","CGC","CCG")]  = 'red'

p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=0, angle=90, hjust=1, colour=labs$colour))


val9.small$repetition=0
val9.small$repetition[val9.small$Group.1==15]=1
rownames(val9.small)=val9.small$codon

dd.row <- as.dendrogram(hclust(dist(val9.small[,"V2",drop=F])))
ddata_x <- dendro_data(dd.row)

p2 <- ggplot(segment(ddata_x)) +
  geom_segment(aes(x=x, y=y, xend=xend, yend=yend))

labs <- label(ddata_x)
labs$colour = 'black'
labs$colour[labs$label %in%  c("CCG/CGC/GCC")]  = 'red'

p2 + geom_text(data=label(ddata_x),
               aes(label=label, x=x, y=0, angle=90, hjust=1.1, colour=labs$colour), size=4)+
  theme_bw()+ylim(-50,100)

