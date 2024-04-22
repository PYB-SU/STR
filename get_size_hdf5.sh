echo "" > size0.txt; for fi in $(ls *.fast5); do h5stat $fi |grep -e "Filename:" -e "  Raw data:" | awk '{if (NR%2==1) {fi=$2} else {si=$3; print fi" "si;}}'>> size0.txt; done

