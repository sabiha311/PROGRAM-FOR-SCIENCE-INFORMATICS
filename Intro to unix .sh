#creating users home directory
cd~
#create directory titled informatics_573
mkdir informatics_573
cd informatics_573
#download all secondary assemblies for human chromosome 1
wget -r-l1-nd-A"chr1_*" -R"chr1.fa.gz"
#unzip all downloaded chromosome1 assemblies
gunzip chr1_*
#creating new empty file called "data_summary.txt:"
touch data_summary.txt
#append detailed information about each file
ls -lh >> data_summary.txt
#append first 10 lines of each assembly
head -n 10 chr1_* >> data_summary.txt
#append the name of each assembly and total number of files
wc chr1_* >> data_summary.txt

