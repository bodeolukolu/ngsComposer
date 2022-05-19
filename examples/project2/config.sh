#General_parameters
###################################################
threads=16
walkaway=True
cluster=False
samples_alt_dir=False
rm_transit=True

#Input_files
###################################################
lib1_R1=A1_R1.fastq.gz
lib1_R2=A1_R2.fastq.gz
lib1_bc=barcodes_lib1.txt
lib2_R1=A2_R1.fastq.gz
lib2_R2=A2_R2.fastq.gz
lib2_bc=barcodes_lib2.txt
lib3_R1=A3_R1.fastq.gz
lib3_R2=A3_R2.fastq.gz
lib3_bc=barcodes_lib3.txt
lib4_R1=A4_R1.fastq.gz
lib4_R2=A4_R2.fastq.gz
lib4_bc=barcodes_lib4.txt

#Tool_parameters
###################################################
front_trim=2
mismatch=1
end_score=20
window=10
min_len=64
adapter_match=12
q_min=20
q_percent=80

#Visualizations
###################################################
QC_final=summary
