cat $1 |  grep -oE "(GT)+(GTCT)+(GC+T)+" | \
    awk '{ligne=$0; match(ligne,/(GT)+/); GT=substr(ligne,RSTART,RLENGTH);ligne2=substr(ligne,RLENGTH+1,length(ligne)-RLENGTH); match(ligne2,/(GTCT)+/); CTGT=substr(ligne2,RSTART,RLENGTH); GCCT=substr(ligne2,RLENGTH+1,length(ligne2)-RLENGTH); if(NR == 1) {print "GT CTGT GCCT subGT subCTGT subGCCT"} print length(GT)/2" "length(CTGT)/4" "length(GCCT)/4" "ligne" "GT" "CTGT" "GCCT}'
