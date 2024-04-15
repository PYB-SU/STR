
liste_DIR = dir()


liste_DIR.DM1 <-
  c("2022_GEN_DM1_pilote05_02-15235", 
    "2022_GEN_DM1_pilote05_03", "2022_GEN_DM1_pilote05_04",
    "2022_GEN_DM1_pilote05_06", "2022_GEN_DM1_pilote05_07",  
    "2022_GEN_DM1_pilote05_09", "2022_GEn_DM1_pilote06_01", "2022_GEN_DM1_pilote06_02", 
    "2022_GEN_DM1_pilote07_01-15235-homemade", "2022_GEN_DM1_pilote07_02-15235-ONT", 
    "2022_GEN_DM1_WG_pilote01",  
    "DM1_pilote03_BC1", "DM1_pilote03_BC2", "DM1_pilote04")


for (d in liste_DIR.DM1) {
  process.gene.files(d,root_dir,"guppy",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
  process.gene.files(d,root_dir,"dorado_hac",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
  process.gene.files(d,root_dir,"bonito_hac",QUALITY = 20, IDENTITY = 0.9, query=c("CTG", "CAG", "CCG", "CGG", "CTC" ),ncommon=1)
}

liste_DIR.OPDM <-
  c("2022_GEN_OPDM_pilote03-PAM21431", 
    "2022_GEN_OPDM_pilote06_Lib1_BC3", "2022_GEN_OPDM_pilote06_Lib1_BC4", 
    "2022_GEN_OPDM_pilote06_Lib1_BC5", "2022_GEN_OPDM_pilote06_Lib2_BC6", 
    "2022_GEN_OPDM_pilote06_Lib2_BC7", "2022_GEN_OPDM_pilote06_Lib2_BC8", 
    "2022_GEN_OPDM_pilote07_lib1", "2022_GEN_OPDM_pilote07_lib2", 
    "2022_GEN_OPDM_pilote08_Lib1_BC5", "2022_GEN_OPDM_pilote08_Lib1_BC6", 
    "2022_GEN_OPDM_pilote08_Lib1_BC7", "2022_GEN_OPDM_pilote08_Lib2_BC8", 
    "2022_GEN_OPDM_pilote08_Lib2_BC9", "2022_GEN_OPDM_pilote09_BC1", 
    "2022_GEN_OPDM_pilote09_BC2", "2022_GEN_OPDM_pilote10", "2023_GEN_OPDM_pilote05_BC1", 
    "2023_GEN_OPDM_pilote05_BC2", "2023_GEN_OPDM_pilote05_BC3", "2023_GEN_OPDM_pilote05_BC4", 
    "OPDM_pilote03", "OPDM_pilote04_BC1", "OPDM_pilote04_BC2")


for (d in liste_DIR.OPDM) {
  process.gene.files(d,root_dir,"guppy",QUALITY = 20, IDENTITY = 0.9, query=c("CGG", "CCG", "CCT"),ncommon=1)
  process.gene.files(d,root_dir,"dorado_hac",QUALITY = 20, IDENTITY = 0.9 , query=c("CGG", "CCG", "CCT"),ncommon=1)
  process.gene.files(d,root_dir,"bonito_hac",QUALITY = 20, IDENTITY = 0.9 , query=c("CGG", "CCG", "CCT"),ncommon=1)
}


# certain I know : 
# OPDM9BC2 gène Notch2NLC: Triplet mettre en 5'-3' codon GGC bleu, AGG en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat. (modifié) 
process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "guppy",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "dorado_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC2",root_dir = root_dir,gene = "NOTCH2NLC", caller = "bonito_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("AGG", "GGC"),ncommon=0,first.codon="CCC",rev.seq = TRUE, repeat.codon="GGC")

# OPDM9BC1 gène GIPC1: Triplet mettre en 3'-5' codon CCG bleu, GGC en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat.
process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "guppy",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "dorado_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

process.one.gene(DIR = "2022_GEN_OPDM_pilote09_BC1",root_dir = root_dir,gene = "GIPC1", caller = "bonito_hac",QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GGC", "CCG"),ncommon=0,first.codon="CTCC",rev.seq = TRUE, repeat.codon="CCG")

# OPDM10 gène LRP12: Triplet mettre en 3'-5' codon CCG bleu, CAG en orange et other. Faire l'alignement à partir du 1 codon qui juxtapose la repeat.
process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "guppy",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "dorado_hac",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")

process.one.gene(DIR = "2022_GEN_OPDM_pilote10",root_dir = root_dir,gene = "LRP12", caller = "bonito_hac",
                 QUALITY = 20, IDENTITY = 0.9, 
                 query=c("GAC", "GCC"),ncommon=0,first.codon="GAC",rev.seq = TRUE, repeat.codon="GCC")

