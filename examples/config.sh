#General_parameters
###################################################
threads=24
walkaway=True
cluster=True
samples_alt_dir=False
rm_transit=True

#Input_files
###################################################
lib1_R1=test1_R1.fastq.gz
lib1_R2=test1_R2.fastq.gz
lib1_bc=barcodes_lib1.txt
lib2_R1=test2_R1.fastq.gz
lib2_R2=test2_R2.fastq.gz
lib2_bc=barcodes_lib2.txt

#Tool_parameters
###################################################
front_trim=6
mismatch=1
R1_motif=TGCATA,TGCATC,TGCATT
R2_motif=CATG
end_score=20
window=10
min_len=64
adapter_match=12
q_min=20
q_percent=80

#Visualizations
###################################################
QC_demultiplexed=summary,full
QC_motif_validated=summary,full
QC_end_trimmed=summary,full
QC_adapter_removed=summary,full
QC_final=summary,full
