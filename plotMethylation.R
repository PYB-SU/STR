library(rtracklayer)

library(GenomeInfoDb)

liste_DIR.DM1 <-
  c("2022_GEN_DM1_pilote05_02-15235", 
    "2022_GEN_DM1_pilote05_03", "2022_GEN_DM1_pilote05_04",
    "2022_GEN_DM1_pilote05_06", "2022_GEN_DM1_pilote05_07",  
    "2022_GEN_DM1_pilote05_09", "2022_GEN_DM1_pilote06_01", "2022_GEN_DM1_pilote06_02", 
    "2022_GEN_DM1_pilote07_01-15235-homemade", "2022_GEN_DM1_pilote07_02-15235-ONT", 
    "2022_GEN_DM1_WG_pilote01",  
    "DM1_pilote03_BC1", "DM1_pilote03_BC2", "DM1_pilote04")


#DM1
res.short.long.DM1 = c()
for (dir in liste_DIR.DM1) {
  res = check.short.long(dir,"dorado_hac", root_dir)
  res.short.long.DM1 = rbind(res.short.long.DM1,res)
}

res.short.long.OPDM =c()
for (dir in liste_DIR.OPDM) {
  res = check.short.long(dir,"dorado_hac", root_dir)
  res.short.long.OPDM = rbind(res.short.long.OPDM,res)
}


split.and.collapse("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000, type="1")
plot.methyl.fancy("2022_GEN_DM1_pilote05_02-15235","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000, type="1")

split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to=48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to=48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_01-15235-homemade","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote07_02-15235-ONT","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_WG_pilote01","dorado_hac","DMPK",methyl = "a_unfiltered",
                  root_dir, from=48594000, to= 48600000)

split.and.collapse("DM1_pilote03_BC1","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("DM1_pilote03_BC1","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "m_unfiltered",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("DM1_pilote03_BC1","dorado_hac","DMPK",methyl = "h_unfiltered",
                  root_dir, from=48594000, to= 48600000)

#05-06
split.and.collapse("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)


plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "m_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_06","dorado_hac","DMPK",methyl = "h_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)


#05-07
split.and.collapse("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "m_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "h_CPGfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)

plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "m_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)
plot.methyl.fancy("2022_GEN_DM1_pilote05_07","dorado_hac","DMPK",methyl = "h_unfiltered",
                  reads="overlap",
                  root_dir, from=48594000, to= 48600000)






# OPDM

split.and.collapse("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir,from=105714000, to = 105718000)

plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote03-PAM21431","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir,from=105714000, to = 105718000)

split.and.collapse("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12","5mCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12","5mCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir,from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote07_lib1","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir,from=105714000, to = 105718000)

split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "m_CPGfiltered",
                  root_dir, from=14620000, to=14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "h_CPGfiltered",
                  root_dir,from=14620000, to = 14625000)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "m_unfiltered",
                  root_dir, from=14620000, to=14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "h_unfiltered",
                  root_dir,from=14620000, to = 14625000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "a_unfiltered",
                  root_dir,from=14620000, to = 14625000)

split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "m_CPGfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "h_CPGfiltered",
                  root_dir,from= 148518000, to=148522000)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "m_unfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "h_unfiltered",
                  root_dir,from= 148518000, to=148522000)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "a_unfiltered",
                  root_dir,from= 148518000, to=148522000)


split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)
split.and.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12","6mA","unfiltered",root_dir = root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)

plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "a_unfiltered",
                  root_dir, from=105714000, to = 105718000)



split.and.collapse("OPDM_pilote03","dorado_hac","LRP12","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote03","dorado_hac","LRP12","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "m_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "h_CPGfiltered",
                  root_dir, from=105714000, to = 105718000)

plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "m_unfiltered",
                  root_dir, from=105714000, to = 105718000)
plot.methyl.fancy("OPDM_pilote03","dorado_hac","LRP12",methyl = "h_unfiltered",
                  root_dir, from=105714000, to = 105718000)





split.and.collapse("OPDM_pilote04_BC1","dorado_hac","GIPC1","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote04_BC1","dorado_hac","GIPC1","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "m_CPGfiltered",
                  root_dir, from=14620000, to = 14625000)
plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "h_CPGfiltered",
                  root_dir, from=14620000, to = 14625000)

plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "m_unfiltered",
                  root_dir, from=14620000, to = 14625000)
plot.methyl.fancy("OPDM_pilote04_BC1","dorado_hac","GIPC1",methyl = "h_unfiltered",
                  root_dir, from=14620000, to = 14625000)


split.and.collapse("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","CPGfiltered",root_dir = root_dir)
split.and.collapse("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC","5mCG_5hmCG","unfiltered",root_dir = root_dir)

plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "m_CPGfiltered",
                  root_dir,from= 148518000, to=148522000)
plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "h_CPGfiltered",
                  root_dir, from= 148518000, to=148522000)

plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "m_unfiltered",
                  root_dir, from= 148518000, to=148522000)
plot.methyl.fancy("OPDM_pilote04_BC2","dorado_hac","NOTCH2NLC",methyl = "h_unfiltered",
                  root_dir, from= 148518000, to=148522000)



fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                  root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","LRP12",methyl = "6mA_collapsed",
                  root_dir)


plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC1","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                  root_dir)

## 9/BC2

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                  root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","LRP12",methyl = "6mA_collapsed",
                  root_dir)


plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote09_BC2","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                  root_dir)

# OPDM10
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","GIPC1",methyl = "6mA_collapsed",
                  root_dir)

plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","LRP12",methyl = "6mA_collapsed",
                  root_dir)


plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)

fix.collapse("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "6mA_unfiltered",methyl_mod = "6mA_collapsed",root_dir = root_dir)
plot.methyl.fancy("2022_GEN_OPDM_pilote10","dorado_hac","NOTCH2NLC",methyl = "6mA_collapsed",
                  root_dir)


#OPDM5

plot.methyl.fancy("2023_GEN_OPDM_pilote05_BC4","dorado_hac","GIPC1",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)


plot.methyl.fancy("2023_GEN_OPDM_pilote05_BC4","dorado_hac","NOTCH2NLC",methyl = "5mCG_5hmCG_CPGfiltered",
                  root_dir)
