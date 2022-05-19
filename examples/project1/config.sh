#General_parameters
###################################################
threads=24
walkaway=True
cluster=False
samples_alt_dir=False
rm_transit=True

#Input_files
###################################################
lib1_R1=pool384_R1.fastq.gz
lib1_R2=pool384_R2.fastq.gz
lib1_bc=barcodes_lib1.txt

#Tool_parameters
###################################################
front_trim=6
mismatch=1
R1_motif=TCC,TCT
R2_motif=TCC,TCT
end_score=20
window=10
min_len=64
adapter_match=12
q_min=20
q_percent=80

#Visualizations
###################################################
QC_final=summary
