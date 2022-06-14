
if [ -z "$threads" ]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi

if [[ -z $walkaway ]]; then
	walkaway=True
fi
if [[ -z $cluster ]]; then
	cluster=False
fi
if [[ -z $samples_alt_dir ]]; then
	samples_alt_dir=True
fi
if [[ -z $rm_transit ]]; then
	rm_transit=True
fi
if [[ -z $front_trim ]]; then
	front_trim=0
fi
if [[ -z $mismatch ]]; then
	mismatch=1
fi
if [[ -z $R1_motif ]]; then
	R1_motif=""
fi
if [[ -z $R2_motif ]]; then
	R2_motif=""
fi
if [[ -z $non_genomic ]]; then
	non_genomic=0
fi
if [[ -z $end_score ]]; then
	end_score=20
fi
if [[ -z $window ]]; then
	window=10
fi
if [[ -z $min_len ]]; then
	min_len=64
fi
if [[ -z $adapter_match ]]; then
	adapter_match=12
fi
if [[ -z $q_min ]]; then
	q_min=20
fi
if [[ -z $q_percent ]]; then
	q_percent=80
fi
if [[ -z $multithread_demultiplex ]]; then
	multithread_demultiplex=False
fi

cd $projdir

test_lib_R2=$(grep '^lib' config.sh | grep '_R2=' | awk '{gsub(/=/,"\t"); print $2}')
if [[ -z "$test_lib_R2" ]]; then
	test_lib_R2=False
fi

nthreads="$(grep threads config.sh)"
nthreads=${nthreads//*=}
totalk=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
loopthreads=2
if [[ "$threads" -gt 1 ]]; then
	N=$((threads/2))
	ram1=$(($totalk/$N))
else
	N=1 && loopthreads=threads
fi
ram1=$((ram1/1000000))
Xmx1=-Xmx${ram1}G
ram2=$(echo "$totalk*0.00000095" | bc)
ram2=${ram2%.*}
Xmx2=-Xmx${ram2}G
if [[ -z "$threads" ]]; then
	threads=$(nproc --all)
	if [[ "$threads" -ge 4 ]]; then
		threads=$((threads-2))
	fi
fi
if  [[ "$threads" -ge 1 ]]; then
	loopthreads=2
	N=$(($threads/2))
else
	N=1 && loopthreads=$threads
fi
if [[ "$threads" -le 4 ]]; then
	gthreads=$threads
	gN=1
else
	gthreads=4
	gN=$(( threads / gthreads ))
fi


if [[ -d "${projdir}/2_demultiplexed/pe" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/pe/*f* 2> /dev/null)" ]]; then
	for i in ${projdir}/2_demultiplexed/pe/*f*; do var=$(echo $i | awk '{gsub(/_R1/,".R1"); gsub(/_R2/,".R2");}1'); mv $i $var; done
fi
if [[ -d "${projdir}/2_demultiplexed/se" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/se/*f* 2> /dev/null)" ]]; then
	for i in ${projdir}/2_demultiplexed/se/*f*; do var=$(echo $i | awk '{gsub(/_R1/,".R1"); gsub(/_R2/,".R2");}1'); mv $i $var; done
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- performing Intitial QC of library/libraries\n${blue}##############################################################################${white}\n"

main_initial_qc() {

	cd ${projdir}
	test_bc=$(grep '^lib' config.sh | grep '_bc' | awk '{gsub(/=/,"\t"); print $2}')
	test_fq=$(grep '^lib' config.sh | grep '_R' | awk '{gsub(/=/,"\t"); print $2}')
	if [[ -z "$test_bc" || -z "$test_fq" ]]; then
		if [[ -d  2_demultiplexed ]]; then
			if [[ -d "${projdir}/2_demultiplexed/pe" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/pe/*f* 2> /dev/null)" ]]; then
				mkdir -p ${projdir}/2_demultiplexed/pe/qc
				cd ${projdir}/2_demultiplexed/pe
				for i in *.f*; do
					python3 $crinoid -r1 $i -t ${threads} -o ./qc & PIDR1=$!
					wait $PIDR1
				done
				wait
			fi
			if [[ -d "${projdir}/2_demultiplexed/se" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/se/*f* 2> /dev/null)" ]]; then
				mkdir -p ${projdir}/2_demultiplexed/se/qc
				cd ${projdir}/2_demultiplexed/se
				for i in *.f*; do
					python3 $crinoid -r1 $i -t ${threads} -o ./qc & PIDR1=$!
					wait $PIDR1
				done
				wait
			fi
		fi
	else

		cd ${projdir}
		mkdir -p 1_initial_qc
		list_lib=$(grep '^lib' config.sh | grep '_R1=' | awk '{gsub(/=/,"\t"); print $2}')
		for li in $list_lib; do
			if [[ "$test_lib_R2" != False ]]; then
				python3 $crinoid -r1 ./samples/${li} -t $((threads/2)) -o ./1_initial_qc & PIDR1=$!
			else
				python3 $crinoid -r1 ./samples/${li} -t $threads -o ./1_initial_qc & PIDR1=$!
			fi
			if [[ "$test_lib_R2" != False ]]; then
				lj=$(echo $li | awk '{gsub(/R1/,"R2"); print}')
				python3 $crinoid -r1 ./samples/${lj} -t $((threads/2)) -o ./1_initial_qc & PIDR2=$!
			fi
			wait $PIDR1
			wait $PIDR2
		done
	fi

	cd ${projdir}/1_initial_qc

	qscore_files=$(ls qscores*.R1*fastq.gz.csv)
	nucleotides_files=$(ls nucleotides*.R1*fastq.gz.csv)
	awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
	awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1' > qscores_initial_qc_R1_summary.csv &&
	awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
	awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1' > nucleotides_initial_qc_R1_summary.csv &&
	Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R1_summary.csv nucleotides_final_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
	wait $PIDR1

	qscore_files=$(ls qscores*.R2*fastq.gz.csv 2> /dev/null/)
	nucleotides_files=$(ls nucleotides*.R2*fastq.gz.csv 2> /dev/null/)
	if [[ ! -z "$qscore_files" ]]; then
		awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
		awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_initial_qc_R2_summary.csv &&
		wait
	fi
	if [[ ! -z "$nucleotides_files" ]]; then
		awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
		awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_initial_qc_R2_summary.csv &&
		wait
	fi
	Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R2_summary.csv nucleotides_final_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages 2> /dev/null & PIDR1=$!
	wait $PIDR1

	for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
	for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
	wait

	mkdir summary
	mv *_summary* ./summary/
	find . -type d -empty -delete

	echo "intial QC complete" > ${projdir}/1_initial_qc_complete
}
cd $projdir
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Do you want to perform intial QC? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping initial QC ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing initial QC ${white}\n"
		time main_initial_qc &>> log.out
	fi
fi
if [ "$walkaway" == True ]; then
	if [ "$initial_qc" == 1 ]; then
		echo -e "${magenta}- performing initial QC ${white}\n"
		time main_initial_qc &>> log.out
	else
		echo -e "${magenta}- skipping initial QC ${white}\n"
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- Demultiplexing of library/libraries\n${blue}##############################################################################${white}\n"

main_demultiplex() {
  cd ${projdir}
  mkdir -p 2_demultiplexed
	mkdir -p ./2_demultiplexed/pe; mkdir -p ./2_demultiplexed/se; mkdir -p ./2_demultiplexed/unknown
  list_lib=$(grep '^lib' config.sh | grep '_R1=' | awk '{gsub(/=/,"\t"); print $2}')

  for li in $list_lib; do
		cd ${projdir}
    bc_matrix=$(grep -h "$li" config.sh | awk '{gsub(/_R1/,"\t"); print $1"_bc"}')
    bc_matrix=$(grep -h "$bc_matrix" config.sh | awk '{gsub(/=/,"\t"); print $2}')
    if [[ "$test_lib_R2" != False ]]; then
			lj=$(echo $li | awk '{gsub(/_R1/,"_R2");gsub(/.R1/,".R2"); }1')
		fi

    # obtain number of forward and reverse barcodes; and then annotate each cell with row and column position
    rows=$(wc -l $bc_matrix | awk '{print $1-1}')
    cols=$(awk -F'\t' '{gsub(/\t\t/,"\t"); gsub(/\n\n/,"\n"); print}' $bc_matrix | awk -F'\t' '{$1=""}1' | awk '{print NF; exit}')
    >| holdbc.txt
    awk '{gsub(/^\t/,"X\t"); gsub(/^ /,"X\t"); print}' $bc_matrix > temp
    column=`head -n 1 temp | wc -w`
    for (( i=1; i <= $column; i++))
    do
      awk '{printf ("%s%s", tab, $'$i'); tab="\t"} END {print ""}' temp
    done >> holdbc.txt
    wait && rm temp
    for i in $(seq 1 $rows); do
      awk -v col=$((i+1)) -v suf=$i 'NR>1{$col=$col"_Row"suf}1' holdbc.txt > temp
      mv temp holdbc.txt; wait
    done
    mv holdbc.txt temp
    column=`head -n 1 temp | wc -w`
    for (( i=1; i <= $column; i++))
    do
      awk '{printf ("%s%s", tab, $'$i'); tab="\t"} END {print ""}' temp
    done >> holdbc.txt
    for i in $(seq 1 $cols); do
      awk -v col=$((i+1)) -v suf=$i 'NR>1{$col=$col"_Column"suf}1' holdbc.txt > temp
      mv temp holdbc.txt; wait
    done


    # trim barcodes to minimum length for demuliplexing
    Min_Flen=$(awk 'NR>1{print $1}' holdbc.txt | awk '$0!=""' | awk '{print $0 "\t" length($0)}' | sort -n -k2,2 | head -n1 | awk '{print $2}')
    Min_Rlen=$(awk 'NR==1' holdbc.txt | awk '{gsub(/\t/,"\n")}1' | awk 'NR>1{print $0 "\t" length($0)}' | sort -n -k2,2 | head -n1 | awk '{print $2}')
    cp holdbc.txt temp && :> ${bc_matrix%.txt}_flush.txt
    column=`head -n 1 temp | wc -w`
    for (( i=1; i <= $column; i++)); do
      awk '{printf ("%s%s", tab, $'$i'); tab="\t"} END {print ""}' temp
    done >> ${bc_matrix%.txt}_flush.txt
    awk -F "\t" -v min=$Min_Rlen 'BEGIN {OFS=FS}; {$1=substr($1, 1, min); print}' ${bc_matrix%.txt}_flush.txt > temp
    rm ${bc_matrix%.txt}_flush.txt
    column=`head -n 1 temp | wc -w`
    for (( i=1; i <= $column; i++)); do
      awk '{printf ("%s%s", tab, $'$i'); tab="\t"} END {print ""}' temp
    done >> ${bc_matrix%.txt}_flush.txt
    awk -F "\t" -v min=$Min_Flen 'BEGIN {OFS=FS}; {$1=substr($1, 1, min); print}' ${bc_matrix%.txt}_flush.txt > temp
    awk '{gsub(/X/,"",$1); gsub(/ /,"\t"); print}' temp > ${bc_matrix%.txt}_flush.txt
		rm temp


    # convert matrix to 3-column dataframe and generate length of trimmed off bases to account for variable length barcode
    awk '{gsub(/X/,"",$1); gsub(/ /,"\t"); print}' holdbc.txt | \
    awk 'NR==1{n=split($0,c);next}{for(i=1;i<=n;i++)s[++t]=$1 FS c[i] FS $(i+1)}END{for(i=1;i<=t;i++){print s[i]}}' | \
    sort -k2,2 | awk '{gsub(/_Row/,"\tRow"); gsub(/_Column/,"\tColumn"); print}' | awk '$4!=""' | sort -n -k3,3 | sort -n -k4,4 | \
    awk -v truncF=$Min_Flen -v truncR=$Min_Rlen 'BEGIN {OFS=FS}; {print substr($1,truncF)"\t"substr($2,truncR)"\t"$3"\t"$4"\t"$5}' | \
    awk '{print $0 "\t" length($1)-1}' | awk '{print $0 "\t" length($2)-1}' | awk '{print $3"_"$4"_"$5"\t"$6"\t"$7}' > ${bc_matrix%.txt}_fringe.txt


		cd ./2_demultiplexed
		if [[ $multithread_demultiplex == False ]]; then
			awk '{gsub(/_Row/,"\t"); gsub(/_Column/,"\t"); print}' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk '{print $1}' | sort | uniq > ${projdir}/cat_RC.txt
			if [[ "$test_lib_R2" != False ]]; then
				python3 $scallop -r1 ${projdir}/samples/"$li" -f $front_trim -o ./ & PIDR1=$!
				python3 $scallop -r1 ${projdir}/samples/"$lj" -f $front_trim -o ./ & PIDR2=$!
				wait $PIDR1
				wait $PIDR2
			else
				python3 $scallop -r1 ${projdir}/samples/"$li" -f $front_trim -o ./ & PIDR1=$!
				wait $PIDR1
			fi
			$gzip trimmed_se* 2> /dev/null/ && wait

			if [[ "$test_lib_R2" != False ]]; then
				python3 $anemone -r1 ./trimmed_se.${li} -r2 ./trimmed_se.${lj} -m $mismatch -c ${projdir}/${bc_matrix%.txt}_flush.txt -o ./ & PIDR1=$!
				wait $PIDR1
			else
				python3 $anemone -r1 ./trimmed_se.${li} -m $mismatch -c ${projdir}/${bc_matrix%.txt}_flush.txt -o ./ & PIDR1=$!
				wait $PIDR1
			fi
			wait

			rm trimmed_se*
			wait
			for sid in $(ls *.R1.fastq | grep -v unknown); do
				fringelen=$( awk -F'\t' -v sampid=${sid%.R1.fastq} '$1 == sampid' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk -F'\t' '{print $2}' )
				if [[ "$fringelen" -gt 0 ]]; then
					python3 $scallop -r1 $sid -f $fringelen && mv ./trimmed_se.${sid} ${sid}
				fi
				$gzip ${sid}
				wait
			done
			if [[ "$test_lib_R2" != False ]]; then
				for sid in $(ls *.R2.fastq | grep -v unknown); do
					fringelen=$( awk -F'\t' -v sampid=${sid%.R2.fastq} '$1 == sampid' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk -F'\t' '{print $3}' )
					if [[ "$fringelen" -gt 0 ]]; then
						python3 $scallop -r1 $sid -f $fringelen && mv ./trimmed_se.${sid} ${sid}
					fi
					$gzip ${sid}
					wait
				done
			fi
			# Now combine fastq files with the same sample_ID
			while IFS="" read -r p || [ -n "$p" ]; do
				find -type f -wholename "./${p}_Row*_Column*R1*" | xargs cat > ${p}.R1.fastq.gz
				if [[ "$test_lib_R2" != False ]]; then
					find -type f -wholename "./${p}_Row*_Column*R2*" | xargs cat > ${p}.R2.fastq.gz
				fi
				rm ${p}_Row*_Column*
			done < ${projdir}/cat_RC.txt
			rm NA.R1.fastq.gz 2> /dev/null &&
			rm NA.R2.fastq.gz 2> /dev/null &&
			rm na.R1.fastq.gz 2> /dev/null &&
			rm na.R2.fastq.gz 2> /dev/null &&
			wait

			mkdir -p unknown
			if [[ ! -z $(ls unknown*.fastq 2> /dev/null) ]]; then
				for i in unknown*.fastq; do $gzip $i 2> /dev/null; done
				for i in unknown*.fastq.gz; do cat ${i} >> ./unknown/${i}; rm ${i} 2> /dev/null; done
			fi
			wait
			if [[ "$test_lib_R2" != False ]]; then
				mv *.fastq.gz ./pe/
			fi
			wait
			if [[ "$test_lib_R2" == False ]]; then
				mv *.fastq.gz ./se/
			fi
			wait

		else
			$zcat ${projdir}/samples/"$li" | awk 'NR%40000000==1{x="R1_chunk"++i".fastq";}{print > x}' - & PIDR1=$!
			wait $PIDR1
			if [[ "$test_lib_R2" != False ]]; then
				$zcat ${projdir}/samples/"$lj" | awk 'NR%40000000==1{x="R2_chunk"++i".fastq";}{print > x}' - & PIDR2=$!
			fi
			wait $PIDR2
			for f in R1_chunk*; do
				subdir=${f%.fastq}
				subdir=${subdir##*_}
				mkdir -- "$subdir"
				mv "R1_${subdir}.fastq" "$subdir"
				if [[ "$test_lib_R2" != False ]]; then mv "R2_${subdir}.fastq" "$subdir"; fi
			done
			wait
			awk '{gsub(/_Row/,"\t"); gsub(/_Column/,"\t"); print}' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk '{print $1}' | sort | uniq > ${projdir}/cat_RC.txt
			for ck in chunk*; do (
				cd $ck
				if [[ "$test_lib_R2" != False ]]; then
					python3 $scallop -r1 ./R1_${ck}.fastq -f $front_trim && rm R1_${ck}.fastq
					python3 $scallop -r1 ./R2_${ck}.fastq -f $front_trim && rm R2_${ck}.fastq
				else
					python3 $scallop -r1 ./R1_${ck}.fastq -f $front_trim && rm R1_${ck}.fastq
				fi

				if [[ "$test_lib_R2" != False ]]; then
					python3 $anemone -r1 ./trimmed_se.R1_${ck}.fastq -r2 ./trimmed_se.R2_${ck}.fastq -m $mismatch -c ${projdir}/${bc_matrix%.txt}_flush.txt -o ./
				else
					python3 $anemone -r1 ./trimmed_se.R1_${ck}.fastq -m $mismatch -c ${projdir}/${bc_matrix%.txt}_flush.txt -o ./
				fi

				rm trimmed_se*
				wait
				mkdir -p unknown
				if [[ ! -z $(ls unknown*.fastq 2> /dev/null) ]]; then
					for i in unknown*.fastq; do $gzip $i 2> /dev/null; done
					for i in unknown*.fastq.gz; do cat ${i} >> ./unknown/${i}; rm ${i} 2> /dev/null; done
				fi
				wait
				for sid in $(ls *.R1.fastq | grep -v unknown); do
					fringelen=$( awk -F'\t' -v sampid=${sid%.R1.fastq} '$1 == sampid' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk -F'\t' '{print $2}' )
					if [[ "$fringelen" -gt 0 ]]; then
						python3 $scallop -r1 $sid -f $fringelen && mv ./trimmed_se.${sid} ${sid}
					fi
					gzip ${sid}
					wait
				done
				if [[ "$test_lib_R2" != False ]]; then
					for sid in $(ls *.R2.fastq | grep -v unknown); do
						fringelen=$( awk -F'\t' -v sampid=${sid%.R2.fastq} '$1 == sampid' ${projdir}/${bc_matrix%.txt}_fringe.txt | awk -F'\t' '{print $3}' )
						if [[ "$fringelen" -gt 0 ]]; then
							python3 $scallop -r1 $sid -f $fringelen && mv ./trimmed_se.${sid} ${sid}
						fi
						gzip ${sid}
						wait
					done
				fi
				# Now combine fastq files with the same sample_ID
				while IFS="" read -r p || [ -n "$p" ]; do
					find -type f -wholename "./${p}_Row*_Column*R1*" | xargs cat > ${p}.R1.fastq.gz
					if [[ "$test_lib_R2" != False ]]; then
						find -type f -wholename "./${p}_Row*_Column*R2*" | xargs cat > ${p}.R2.fastq.gz
					fi
					rm ${p}_Row*_Column*
				done < ${projdir}/cat_RC.txt
				cd ../
				) &
				if [[ $(jobs -r -p | wc -l) -ge $N ]]; then
					wait
				fi
			done
			wait
			for ck in chunk*; do mv $ck ${bc_matrix%.txt}_${ck}; done
		fi

		rm ${projdir}/${bc_matrix%.txt}_fringe.txt
		rm ${projdir}/${bc_matrix%.txt}_flush.txt
		rm ${projdir}/holdbc.txt
		rm ${projdir}/cat_RC.txt

  done
  wait


	if [[ $multithread_demultiplex == False ]]; then
		:
	else
		rm NA.R1.fastq.gz 2> /dev/null &&
		rm NA.R2.fastq.gz 2> /dev/null &&
		rm na.R1.fastq.gz 2> /dev/null &&
		rm na.R2.fastq.gz 2> /dev/null &&
		wait
		find -type f -wholename "./*chunk*/unknown.R1.fastq.gz" | xargs cat > ./unknown/unknown.R1.fastq.gz & PIDR1=$!
		wait $PIDR1
		if [[ "$test_lib_R2" != False ]]; then
			find -type f -wholename "./*chunk*/unknown.R2.fastq.gz" | xargs cat > ./unknown/unknown.R2.fastq.gz  & PIDR2=$!
		fi
		wait $PIDR2
		rm ./*chunk*/unknown.R1.fastq.gz ./*chunk*/unknown.R2.fastq.gz
		wait
		samples_r1=$(find -type f -wholename "./*/*R1*" | awk '{gsub(/\//,"\t"); print}' | awk '{print $3}' | sort | uniq | grep -v 'unknown' | grep -v 'qc')
		for f in $samples_r1; do (
			if [[ "$(ls -A ./*chunk*/*R2.fastq.gz 2> /dev/null)" ]]; then
				find ./*chunk*/${f} | xargs cat > ./pe/${f}
				find ./*chunk*/${f%.R1.fastq.gz}.R2.fastq.gz | xargs cat > ./pe/${f%.R1.fastq.gz}.R2.fastq.gz
			else
				find ./*chunk*/${f} | xargs cat > ./se/${f}
			fi
			wait
			rm ./*chunk*/${f} ./*chunk*/${f%.R1.fastq.gz}.R2.fastq.gz
			wait ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
	fi


	find . -type d -empty -delete
	find ./pe -size 0 -delete 2> /dev/null
	find ./se -size 0 -delete 2> /dev/null
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done; fi

	if [[ "${QC_demultiplexed}" =~ summary || "${QC_demultiplexed}" =~ full ]]; then
		cd unknown
		mkdir -p qc
		python3 $crinoid -r1 ./unknown.R1.fastq.gz -t "${threads}" -o ./qc & PIDR1=$!
		wait $PIDR1
		if [[ "$test_lib_R2" != False ]]; then
			python3 $crinoid -r1 ./unknown.R2.fastq.gz -t "${threads}" -o ./qc & PIDR1=$!
			wait $PIDR1
		fi
		for i in ./qc/nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
		for i in ./qc/qscores.*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
		wait

		if [[ "$test_lib_R2" != False ]]; then
			cd ../pe
		else
			cd ../se
		fi
		mkdir -p qc
		for f in *.R1.fastq.gz; do (
			python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
		wait
		if [[ "$test_lib_R2" != False ]]; then
			for f in *.R2.fastq.gz; do (
				python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
			done
			wait
		fi

		cd ./qc
		if [[ ! "${QC_demultiplexed}" =~ full ]]; then
			rm *.png
		fi
		if [[ "${QC_demultiplexed}" =~ summary ]]; then
			qscore_files=$(ls qscores*.R1.fastq.gz.csv)
			nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
			awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
			awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_demultiplexed_R1_summary.csv &&
			awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
			awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_demultiplexed_R1_summary.csv &&
			Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_demultiplexed_R1_summary.csv nucleotides_demultiplexed_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages  & PIDR1=$!
			wait $PIDR1

			if [[ "$test_lib_R2" != False ]]; then
				qscore_files=$(ls qscores*.R2.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_demultiplexed_R2_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_demultiplexed_R2_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_demultiplexed_R2_summary.csv nucleotides_demultiplexed_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1
			fi
		fi
		for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
		for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
		wait

		mkdir full summary
		mv *fastq* ./full/
		mv *_summary* ./summary/
		find . -type d -empty -delete
	fi
	echo "demultiplexed complete" > ${projdir}/2_demultiplexed_complete

}
cd $projdir
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Do you want to perform demultiplexing? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping demultiplexing ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing demultiplexing ${white}\n"
		if [[ ! -d ${projdir}/2_demultiplexed ]]; then
			time main_demultiplex &>> log.out
		else
			echo -e "${magenta}- ${projdir}/2_demultiplexed already exist ${white}"
			echo -e "${magenta}- skipping demultiplexing in 10 seconds ${white}\n"
			sleep 10
			if [[ "$(ls ${projdir}/2_demultiplexed/*.f*)" =~ R2 ]]; then
				mkdir -p ${projdir}/2_demultiplexed/pe; mv *.f* ${projdir}/2_demultiplexed/pe/
				wait
			else
				mkdir -p ${projdir}/2_demultiplexed/se; mv *.f* ${projdir}/2_demultiplexed/se/
				wait
			fi
		fi
	fi
fi
if [ "$walkaway" == True ]; then
	if [ "$demultiplex" == 1 ]; then
		echo -e "${magenta}- performing demultiplexing ${white}\n"
		if [[ ! -d ${projdir}/2_demultiplexed ]]; then
			time main_demultiplex &>> log.out
		else
			echo -e "${magenta}- ${projdir}/2_demultiplexed already exist ${white}"
			echo -e "${magenta}- skipping demultiplexing in 10 seconds ${white}\n"
			sleep 10
			if [[ "$(ls ${projdir}/2_demultiplexed/*.f*)" =~ R2 ]]; then
				mkdir -p ${projdir}/2_demultiplexed/pe; mv *.f* ${projdir}/2_demultiplexed/pe/
				wait
			else
				mkdir -p ${projdir}/2_demultiplexed/se; mv *.f* ${projdir}/2_demultiplexed/se/
				wait
			fi
		fi
	else
		echo -e "${magenta}- skipping demultiplexing ${white}\n"
	fi
fi
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Demultiplexing was performed with a mismatch="$mismatch" ${white}"
	echo -e "${magenta}- Check output and QC plot(s) ${white}"
	echo -e "${magenta}- Do you want to change the mismatch value and re-run demultiplexing? ${white}"
	read -p "- n(NO) or <enter_value>? " -n 1 -r
  if [[ $REPLY =~ ^[Nn]$ ]]; then
    printf "\n"
    echo -e "${magenta}- ngsComposer will proceed to the next analytical step ${white}"
  else
    if [[ $REPLY -eq "$mismatch" ]]; then
      printf "\n"
      echo -e "${magenta}- The same mismatch value specified ${white}"
      echo -e "${magenta}- ngsComposer will proceed to the next analytical step ${white}"
      sleep 5
    else
			echo -e "${magenta}- mismatch=$REPLY ${white}"
      mismatch=$REPLY; printf "\n"; time main_demultiplex &>> log.out
    fi
  fi
fi




######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- performing filtering based on known motif validation \n${blue}##############################################################################${white}\n"

main_motif_validation() {

  cd ${projdir}
  mkdir -p 3_motif_validated
  cd 3_motif_validated
	mkdir -p pe se
  motifR1=$( echo $R1_motif | awk '{gsub(/,/," "); print}')
  motifR2=$( echo $R2_motif | awk '{gsub(/,/," "); print}')

	if [[ -d "${projdir}/2_demultiplexed/pe" ]] && [[ ! -z "$motifR1" ]] && [[ ! -z "$motifR2" ]]; then
		for mot in ${projdir}/2_demultiplexed/pe/*.R1.fastq.gz; do (
		    python3 $rotifer -r1 ${mot} -r2 ${mot%.R1.fastq.gz}.R2.fastq.gz -m1 $motifR1 -m2 $motifR2 -o ./pe &&
				motgz=${mot#*/pe/}; motgz=${motgz%.gz} &&
				cd pe &&
				$gzip ./*${motgz} && $gzip ./*${motgz%.R1.fastq}.R2.fastq &&
				cd ../ &&
				wait ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
		done
		wait
	fi

	if [[ -d "${projdir}/2_demultiplexed/pe" ]] && [[ ! -z "$motifR1" ]] && [[ -z "$motifR2" ]]; then
		for mot in ${projdir}/2_demultiplexed/pe/*.R1.fastq.gz; do (
				python3 $rotifer -r1 ${mot} -r2 ${mot%.R1.fastq.gz}.R2.fastq.gz -m1 $motifR1 -o ./pe &&
				motgz=${mot#*/pe/}; motgz=${motgz%.gz} &&
				cd pe &&
				$gzip ./*${motgz} && $gzip ./*${motgz%.R1.fastq}.R2.fastq &&
				cd ../ &&
				wait ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
		done
		wait
	fi

	if [[ -d "${projdir}/2_demultiplexed/pe" ]] && [[ -z "$motifR1" ]] && [[ ! -z "$motifR2" ]]; then
		for mot in ${projdir}/2_demultiplexed/pe/*.R1.fastq.gz; do (
				python3 $rotifer -r1 ${mot} -r2 ${mot%.R1.fastq.gz}.R2.fastq.gz -m2 $motifR2 -o ./pe &&
				motgz=${mot#*/pe/}; motgz=${motgz%.gz} &&
				cd pe &&
				$gzip ./*${motgz} && $gzip ./*${motgz%.R1.fastq}.R2.fastq &&
				cd ../ &&
				wait ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
		done
		wait
	fi

	if [[ -d "${projdir}/2_demultiplexed/se" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/se 2> /dev/null)" ]]; then
		for mot in ${projdir}/2_demultiplexed/se/*.R1.fastq.gz; do (
				python3 $rotifer -r1 ${mot} -m1 $motifR1 -o ./se &&
				motgz=${mot#*/se/}; motgz=${motgz%.gz} &&
				cd se &&
				$gzip ./*${motgz} &&
				cd ../ &&
				wait ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
		done
		wait
	fi
	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i; fi; done &&
	if [[ -d "${projdir}/2_demultiplexed/se" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/se 2> /dev/null)" ]]; then
		for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done &&
		wait
	else
		:
	fi


	mv se pre_se
	mkdir -p se
	mot_se=$( ls ./pe/se* | cat - <(ls ./pre_se/se* 2> /dev/null) | awk '{gsub(/ /,"\n"); gsub(/.\/pe\//,""); gsub(/.\/pre_se\//,""); gsub(/se./,"");}1' | sort | uniq)
	if [[ ! "$(ls -A ./pre_se/ 2> /dev/null)" ]]; then
		rmdir pre_se
	fi
	for i in $mot_se; do
		if [[ -f ./pe/se.${i} && -f ./pre_se/se.${i} ]]; then cat ./pe/se.${i} ./pre_se/se.${i} > ./se/$i; rm ./pe/se.${i}; fi &&
		if [[ -f ./pe/se.${i} && ! -f ./pre_se/se.${i} ]]; then mv ./pe/se.${i} ./se/$i; fi &&
		if [[ ! -f ./pe/se.${i} && -f ./pre_se/se.${i} ]]; then mv ./pre_se/se.${i} ./se/$i; fi &&
		wait
	done
	wait
	rm -rf pre_se 2> /dev/null
	cd ./pe
	if [[ ! -z "$(ls pe.* 2> /dev/null)" ]]; then
		for i in pe.*; do mv $i ${i#pe.}; done
	fi
	cd ../
	find . -type d -empty -delete
	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i; fi; done &&
	if [[ -d "${projdir}/2_demultiplexed/se" ]] && [[ "$(ls -A ${projdir}/2_demultiplexed/se 2> /dev/null)" ]]; then
		for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done &&
		wait
	fi
	wait

	if [[ "${QC_motif_validated}" =~ summary || "${QC_motif_validated}" =~ full ]]; then
		if [[ -d ${projdir}/3_motif_validated/pe ]]; then
			cd ${projdir}/3_motif_validated/pe
			mkdir -p qc
			for f in *.R1.fastq.gz; do (
				python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
				python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
			done
			wait

			cd ./qc
			if [[ ! "${QC_motif_validated}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_motif_validated}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_motif_valid_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_motif_valid_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_motif_valid_R1_summary.csv nucleotides_motif_valid_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_motif_valid_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_motif_valid_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_motif_valid_R2_summary.csv nucleotides_motif_valid_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi

		if [[ -d ${projdir}/3_motif_validated/se && "$test_lib_R2" != False ]]; then
			cd ${projdir}/3_motif_validated/se
			mkdir -p qc
			if [[ "$test_lib_R2" != False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
					python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi
			if [[ -d se && "$test_lib_R2" == False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi

			cd ./qc
			if [[ ! "${QC_motif_validated}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_motif_validated}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_motif_valid_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_motif_valid_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_motif_valid_R1_summary.csv nucleotides_motif_valid_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_motif_valid_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_motif_valid_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_motif_valid_R2_summary.csv nucleotides_motif_valid_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi
	fi
	echo "motif validation complete" > ${projdir}/3_motif_validation_complete

}
cd $projdir
if [[ -z "$motifR1" ]] && [[ -z "$motifR2" ]]; then
	if [ "$walkaway" == False ]; then
		echo -e "${magenta}- Do you want to perform known motif validation? ${white}\n"
		read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
		if [[ ! $REPLY =~ ^[Yy]$ ]]; then
			printf '\n'
			echo -e "${magenta}- skipping known motif validation ${white}\n"
		else
			printf '\n'
			echo -e "${magenta}- performing known motif validation ${white}\n"
			time main_motif_validation &>> log.out
		fi
	fi
	if [ "$walkaway" == True ]; then
		if [ "$motif_validation" == 1 ]; then
			echo -e "${magenta}- performing known motif validation ${white}\n"
			time main_motif_validation &>> log.out
		else
			echo -e "${magenta}- skipping known motif validation ${white}\n"
		fi
	fi

	if [ "$walkaway" == False ]; then
		echo -e "${magenta}- motif validation was performed with R1_motif("$R1_motif") and R2_motif("$R2_motif") ${white}"
		echo -e "${magenta}- Check output and QC plot(s) ${white}"
		echo -e "${magenta}- Do you want to accept R1_motif validation? ${white}"
		read -p "- y(YES) or <enter_new_motif>? " -n 1 -r
		if [[ $REPLY =~ ^[Yy]$ ]]; then
			echo -e "${magenta}- R1_motif=$R1_motif_new ${white}"
			R1_motif_new=""; printf "\n"
		else
			R1_motif_new=$REPLY; printf "\n"
		fi
		echo -e "${magenta}- Do you want to accept R2_motif validation? ${white}"
		read -p "- y(YES) or <enter_new_motif>? " -n 1 -r
		if [[ $REPLY =~ ^[Yy]$ ]]; then
			echo -e "${magenta}- R2_motif=$R2_motif_new ${white}"
			R2_motif_new=""; printf "\n"
		else
			R2_motif_new=$REPLY; printf "\n"
		fi
		if [[ -z "$R1_motif_new" && -z "$R2_motif_new" ]]; then
			echo -e "${magenta}- ngsComposer will replace output from motif_validation with demultiplexing output, i.e. skipping motif_validation ${white}"
			rm -rf ${projdir}/3_motif_validated/pe 2> /dev/null
			rm -rf ${projdir}/3_motif_validated/se 2> /dev/null
			mv ${projdir}/2_demultiplexed/pe ${projdir}/3_motif_validated/
		else
			if [[ "$R1_motif" == "$R1_motif_new" && "$R2_motif" == "$R2_motif_new" ]]; then
				echo -e "${magenta}- ngsComposer will proceed to the next analytical step ${white}"
			else
				R1_motif=$R1_motif_new
				R2_motif=$R2_motif_new
				time main_motif_validation &>> log.out
			fi
		fi
	fi

	if [[ "$rm_transit" == True ]] && [[ -f "${projdir}/3_motif_validation_complete" ]]; then
		rm ${projdir}/2_demultiplexed/pe/*fastq* 2> /dev/null
	fi
fi


######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- performing end-trimming of reads \n${blue}##############################################################################${white}\n"

main_end_trim() {

  cd ${projdir}
  mkdir -p 4_end_trimmed
	cd 4_end_trimmed
	mkdir -p pe se

	if [[ ! -d ${projdir}/3_motif_validated ]]; then mkdir -p ${projdir}/3_motif_validated/pe; mkdir -p ${projdir}/3_motif_validated/se; fi
	if [[ -z "$(ls -A ${projdir}/3_motif_validated/pe 2> /dev/null)" ]]; then mv ${projdir}/2_demultiplexed/pe/*fastq* ${projdir}/3_motif_validated/pe/ 2> /dev/null && 	find ${projdir}/3_motif_validated/pe/ -type d -empty -delete; fi
	if [[ -z "$(ls -A ${projdir}/3_motif_validated/se 2> /dev/null)" ]]; then mv ${projdir}/2_demultiplexed/se/*fastq* ${projdir}/3_motif_validated/se/ 2> /dev/null && 	find ${projdir}/3_motif_validated/se/ -type d -empty -delete; fi
	if [[ "$rm_transit" == True ]]; then
		rm ${projdir}/2_demultiplexed/pe/*fastq* 2> /dev/null
		rm ${projdir}/2_demultiplexed/se/*fastq* 2> /dev/null
	fi
	find ${projdir}/3_motif_validated -type d -empty -delete

	if [[ -d "${projdir}/3_motif_validated/pe" ]]; then
		for etm in ${projdir}/3_motif_validated/pe/*.R1.fastq.gz; do (
			python3 $scallop -r1 ${etm} -r2 ${etm%.R1.fastq.gz}.R2.fastq.gz -e $end_score -w $window -l $min_len -o ./pe/ &&
			etmgz=${etm#*/pe/}; etmgz=${etmgz%.gz} &&
			cd pe &&
			$gzip ./*${etmgz} && $gzip ./*${etmgz%.R1.fastq}.R2.fastq &&
			cd ../ &&
			wait ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
		wait
	fi
	if [[ -d "${projdir}/3_motif_validated/se" ]]; then
		for etm in ${projdir}/3_motif_validated/se/*.R1.fastq.gz; do (
			python3 $scallop -r1 ${etm} -e $end_score -w $window -l $min_len -o ./se/ &&
			etmgz=${etm#*/se/}; etmgz=${etmgz%.gz} &&
			cd se &&
			$gzip ./*${etmgz} &&
			cd ../ &&
			wait

			if [[ -f ${etm%.R1.fastq.gz}.R2.fastq.gz ]]; then
				python3 $scallop -r1 ${etm%.R1.fastq.gz}.R2.fastq.gz -e $end_score -w $window -l $min_len -o ./se/ &&
				cd se &&
				$gzip ./*${etmgz%.R1.fastq}.R2.fastq &&
				cd ../ &&
				wait
			fi  ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
	fi

	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	mv se pre_se
	mkdir -p se
	etm_se=$( ls ./pe/trimmed_se* 2> /dev/null | cat - <(ls ./pre_se/trimmed_se* 2> /dev/null) | awk '{gsub(/ /,"\n"); gsub(/.\/pe\//,""); gsub(/.\/pre_se\//,""); gsub(/trimmed_se./,"");}1' | sort | uniq)
	for i in $etm_se; do
		if [[ -f ./pe/trimmed_se.${i} && -f ./pre_se/trimmed_se.${i} ]]; then cat ./pe/trimmed_se.${i} ./pre_se/trimmed_se.${i} > ./se/$i; fi &&
		if [[ -f ./pe/trimmed_se.${i} && ! -f ./pre_se/trimmed_se.${i} ]]; then mv ./pe/trimmed_se.${i} ./se/$i; fi &&
		if [[ ! -f ./pe/trimmed_se.${i} && -f ./pre_se/trimmed_se.${i} ]]; then mv ./pre_se/trimmed_se.${i} ./se/$i; fi &&
		wait
	done
	wait
	rm -rf pre_se 2> /dev/null
	rm ./pe/trimmed_se* 2> /dev/null
	cd ./pe
	if [[ ! -z "$(ls trimmed_pe.* 2> /dev/null)" ]]; then
		for i in trimmed_pe.*; do mv $i ${i#trimmed_pe.}; done
	fi
	cd ../
	find . -type d -empty -delete
	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	if [[ "${QC_end_trimmed}" =~ summary || "${QC_end_trimmed}" =~ full ]]; then
		if [[ -d ${projdir}/4_end_trimmed/pe ]]; then
			cd ${projdir}/4_end_trimmed/pe
			mkdir -p qc
			for f in *.R1.fastq.gz; do (
				python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
				python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
			done
			wait

			cd ./qc
			if [[ ! "${QC_end_trimmed}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_end_trimmed}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_end_trimmed_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_end_trimmed_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_end_trimmed_R1_summary.csv nucleotides_end_trimmed_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_end_trimmed_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_end_trimmed_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_end_trimmed_R2_summary.csv nucleotides_end_trimmed_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi

		if [[ -d ${projdir}/4_end_trimmed/se && "$test_lib_R2" != False ]]; then
			cd ${projdir}/4_end_trimmed/se
			mkdir -p qc
			if [[ "$test_lib_R2" != False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
					python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi
			if [[ -d se && "$test_lib_R2" == False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi

			cd ./qc
			if [[ ! "${QC_end_trimmed}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_end_trimmed}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_end_trimmed_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_end_trimmed_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_end_trimmed_R1_summary.csv nucleotides_end_trimmed_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_end_trimmed_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_end_trimmed_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_end_trimmed_R2_summary.csv nucleotides_end_trimmed_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi
	fi
	echo "end trimming complete" > ${projdir}/4_end_trimming_complete

}
cd $projdir
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Do you want to perform end-trimming of reads? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping end-trimming of reads ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing end-trimming of reads ${white}\n"
		time main_end_trim &>> log.out
	fi
fi
if [ "$walkaway" == True ]; then
	if [ "$end_trim" == 1 ]; then
		echo -e "${magenta}- performing end-trimming of reads ${white}\n"
		time main_end_trim &>> log.out
	else
		echo -e "${magenta}- skipping end-trimming of reads ${white}\n"
	fi
fi
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- end trimming was performed with end_score="$end_score", window="$window", and min_len="$min_len" ${white}"
	echo -e "${magenta}- Check output and QC plot(s) ${white}"
	echo -e "${magenta}- Do you want to change end_score value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		echo -e "${magenta}- end_score_new=$REPLY ${white}"
		end_score_new=$REPLY
	else
		end_score_new=$end_score
	fi
	echo -e "${magenta}\n- Do you want to change window size? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		echo -e "${magenta}- window_new=$REPLY ${white}"
		window_new=$REPLY
	else
		window_new=$window
	fi
	echo -e "${magenta}\n- Do you want to change min_len value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		echo -e "${magenta}- min_len_new=$REPLY ${white}"
		min_len_new=$REPLY
	else
		min_len_new=$min_len
	fi
	if [[ "$end_score" == "$end_score_new" && "$window" == "$window_new" && "$min_len" == "$min_len_new" ]]; then
		echo -e "${magenta}- ngsComposer will proceed to the next analytical step ${white}"
	else
		end_score=$end_score_new
		window=$window_new
		min_len=$min_len_new
		time main_end_trim &>> log.out
	fi
fi
if [[ "$rm_transit" == True ]] && [[ -f "${projdir}/4_end_trimming_complete" ]]; then
	rm ${projdir}/3_motif_validated/pe/*fastq* 2> /dev/null
	rm ${projdir}/3_motif_validated/se/*fastq* 2> /dev/null
fi




######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- performing adapter removal \n${blue}##############################################################################${white}\n"

main_adapter_remove() {

  cd ${projdir}
	if [[ -f adapters_R1.txt ]]; then mv adapters_R1.txt adapters.R1.txt; fi
	if [[ -f adapters_R2.txt ]]; then mv adapters_R2.txt adapters.R2.txt; fi
  mkdir -p 5_adapter_removed
	cd 5_adapter_removed
	mkdir -p pe se

	if [[ ! -d ${projdir}/4_end_trimmed/pe ]]; then mkdir -p ${projdir}/4_end_trimmed/pe; mkdir -p ${projdir}/4_end_trimmed/se; fi
	if [[ -z "$(ls -A ${projdir}/4_end_trimmed/pe 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/3_motif_validated/pe 2> /dev/null)" ]] ; then
			mv ${projdir}/2_demultiplexed/pe/*fastq* ${projdir}/4_end_trimmed/pe/ 2> /dev/null && 	find ${projdir}/4_end_trimmed/pe -type d -empty -delete
		else
			mv ${projdir}/3_motif_validated/pe/*fastq* ${projdir}/4_end_trimmed/pe/ 2> /dev/null && 	find ${projdir}/4_end_trimmed/pe -type d -empty -delete
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/4_end_trimmed/se 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/3_motif_validated/se 2> /dev/null)" ]] ; then
			mv ${projdir}/2_demultiplexed/se/*fastq* ${projdir}/4_end_trimmed/se/ 2> /dev/null && 	find ${projdir}/4_end_trimmed/se -type d -empty -delete
		else
			mv ${projdir}/3_motif_validated/se/*fastq* ${projdir}/4_end_trimmed/se/ 2> /dev/null && 	find ${projdir}/4_end_trimmed/se -type d -empty -delete
		fi
	fi

	if [[ "$rm_transit" == True ]]; then
		rm ${projdir}/2_demultiplexed/pe/*fastq* 2> /dev/null
		rm ${projdir}/2_demultiplexed/se/*fastq* 2> /dev/null
		rm ${projdir}/3_motif_validated/pe/*fastq* 2> /dev/null
		rm ${projdir}/3_motif_validated/se/*fastq* 2> /dev/null
	fi
	find ${projdir}/4_end_trimmed -type d -empty -delete


	if [[ -d "${projdir}/4_end_trimmed/pe" ]]; then
		for adp in ${projdir}/4_end_trimmed/pe/*.R1.fastq.gz; do (
			python3 $porifera -r1 ${adp} -r2 ${adp%.R1.fastq.gz}.R2.fastq.gz -a1 ${projdir}/adapters.R2.txt -a2 ${projdir}/adapters.R1.txt -m $adapter_match -l $min_len -o ./pe/ &&
			adpgz=${adp#*/pe/}; adpgz=${adpgz%.gz} &&
	    cd pe &&
	    $gzip ./*${adpgz} && $gzip ./*${adpgz%.R1.fastq}.R2.fastq &&
	    cd ../ &&
			wait ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
		wait
	fi
	if [[ -d "${projdir}/4_end_trimmed/se" ]]; then
		for adp in ${projdir}/4_end_trimmed/se/*.R1.fastq.gz; do (
			python3 $porifera -r1 ${adp} -a1 ${projdir}/adapters.R2.txt -m $adapter_match -l $min_len -o ./se/ &&
			adpgz=${adp#*/se/}; adpgz=${adpgz%.gz} &&
	    cd se &&
	    $gzip ./*${adpgz} &&
	    cd ../ &&
	    wait
	    if [[ -f ${adp%.R1.fastq.gz}.R2.fastq.gz ]]; then
				python3 $porifera -r1 ${adp%.R1.fastq.gz}.R2.fastq.gz -a1 ${projdir}/adapters.R1.txt -m $adapter_match -l $min_len -o ./se/ &&
	      cd se &&
	      $gzip ./*${adpgz%.R1.fastq}.R2.fastq &&
	      cd ../ &&
	      wait
	    fi  ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
	fi

	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	if [[ -d "${projdir}/4_end_trimmed/se" ]]; then
		mv se pre_se
		mkdir -p se
		adp_se=$( ls ./pe/se.adapted* 2> /dev/null | cat - <(ls ./pre_se/adapted*) | awk '{gsub(/ /,"\n"); gsub(/.\/pe\//,""); gsub(/.\/pre_se\//,""); gsub(/se.adapted./,"");}1' | awk '{gsub(/adapted./,"");}1' | sort | uniq)
		for i in $adp_se; do
			if [[ -f ./pe/se.adapted.${i} && -f ./pre_se/adapted.${i} ]]; then cat ./pe/se.adapted.${i} ./pre_se/adapted.${i} > ./se/$i; fi &&
			if [[ -f ./pe/se.adapted.${i} && ! -f ./pre_se/adapted.${i} ]]; then mv ./pe/se.adapted.${i} ./se/$i; fi &&
			if [[ ! -f ./pe/se.adapted.${i} && -f ./pre_se/adapted.${i} ]]; then mv ./pre_se/adapted.${i} ./se/$i; fi &&
			wait
		done
	fi

	rm -rf pre_se 2> /dev/null
	rm ./pe/se.adapted*  2> /dev/null
	cd ./pe
	if [[ ! -z "$(ls pe.adapted.* 2> /dev/null)" ]]; then
		for i in pe.adapted.*; do mv $i ${i#pe.adapted.}; done
	fi
	cd ../se
	if [[ -d "${projdir}/4_end_trimmed/se" ]]; then for i in se.adapted.*; do mv $i ${i#se.adapted.} 2> /dev/null; done; fi
	cd ../
	find . -type d -empty -delete
	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	if [[ "${QC_adapter_removed}" =~ summary || "${QC_adapter_removed}" =~ full ]]; then
		if [[ -d ${projdir}/5_adapter_removed/pe ]]; then
			cd ${projdir}/5_adapter_removed/pe
			mkdir -p qc
			for f in *.R1.fastq.gz; do (
				python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
				python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
			done
			wait

			cd ./qc
			if [[ ! "${QC_adapter_removed}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_adapter_removed}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_adapter_removed_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_adapter_removed_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_adapter_removed_R1_summary.csv nucleotides_adapter_removed_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_adapter_removed_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_adapter_removed_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_adapter_removed_R2_summary.csv nucleotides_adapter_removed_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi

		if [[ -d ${projdir}/5_adapter_removed/se && "$test_lib_R2" != False ]]; then
			cd ${projdir}/5_adapter_removed/se
			mkdir -p qc
			if [[ "$test_lib_R2" != False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
					python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi
			if [[ -d se && "$test_lib_R2" == False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi

			cd ./qc
			if [[ ! "${QC_adapter_removed}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_adapter_removed}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_adapter_removed_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_adapter_removed_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_adapter_removed_R1_summary.csv nucleotides_adapter_removed_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_adapter_removed_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_adapter_removed_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_adapter_removed_R2_summary.csv nucleotides_adapter_removed_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi
	fi
	echo "adapter removal complete" > ${projdir}/5_adapter_removal_complete

}
cd $projdir
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Do you want to perform adapter removal? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping adapter removal ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing adapter removal ${white}\n"
		time main_adapter_remove &>> log.out
	fi
fi
if [ "$walkaway" == True ]; then
	if [ "$adapter_remove" == 1 ]; then
		echo -e "${magenta}- performing adapter removal ${white}\n"
		time main_adapter_remove &>> log.out
	else
		echo -e "${magenta}- skipping adapter removal ${white}\n"
	fi
fi
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- adapter removal was performed with adapter_match="$adapter_match" min_len="$min_len" ${white}"
	echo -e "${magenta}- Check output and QC plot(s) ${white}"
	echo -e "${magenta}- Do you want to change adapter_match value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		echo -e "${magenta}- adapter_match_new=$REPLY ${white}"
		adapter_match_new=$REPLY
	else
		adapter_match_new=$adapter_match
	fi
	echo -e "${magenta}\n- Do you want to change min_len value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		echo -e "${magenta}- min_len_new=$REPLY ${white}"
		min_len_new=$REPLY
	else
		min_len_new=$min_len
	fi
	if [[ "$adapter_match" == "$adapter_match_new" && "$min_len" == "$min_len_new" ]]; then
		echo -e "${magenta}\n- ngsComposer will proceed to the next analytical step ${white}"
	else
		adapter_match=$adapter_match_new
		min_len=$min_len_new
		time main_adapter_remove &>> log.out
	fi
fi
if [[ "$rm_transit" == True ]] && [[ -f "${projdir}/5_adapter_removal_complete" ]]; then
	rm ${projdir}/4_end_trimmed/pe/*fastq* 2> /dev/null
	rm ${projdir}/4_end_trimmed/se/*fastq* 2> /dev/null
fi



######################################################################################################################################################
echo -e "${blue}\n############################################################################## ${yellow}\n- performing read quality-filtering \n${blue}##############################################################################${white}\n"

main_quality_filter() {

  cd ${projdir}
  mkdir -p 6_quality_filtered_final
	cd 6_quality_filtered_final
	mkdir -p pe se

	if [[ ! -d ${projdir}/5_adapter_removed/pe ]]; then mkdir -p ${projdir}/5_adapter_removed/pe; mkdir -p ${projdir}/5_adapter_removed/se; fi
	if [[ -z "$(ls -A ${projdir}/5_adapter_removed/pe 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/4_end_trimmed/pe 2> /dev/null)" ]] ; then
			if [[ -z "$(ls -A ${projdir}/3_motif_validated/pe 2> /dev/null)" ]] ; then
				mv ${projdir}/2_demultiplexed/pe/*fastq* ${projdir}/5_adapter_removed/pe/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/pe -type d -empty -delete
			else
				mv ${projdir}/3_motif_validated/pe/*fastq* ${projdir}/5_adapter_removed/pe/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/pe -type d -empty -delete
			fi
		else
			mv ${projdir}/4_end_trimmed/pe/*fastq* ${projdir}/5_adapter_removed/pe/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/pe -type d -empty -delete
		fi
	fi
	if [[ -z "$(ls -A ${projdir}/5_adapter_removed/se 2> /dev/null)" ]]; then
		if [[ -z "$(ls -A ${projdir}/4_end_trimmed/se 2> /dev/null)" ]] ; then
			if [[ -z "$(ls -A ${projdir}/3_motif_validated/se 2> /dev/null)" ]] ; then
				mv ${projdir}/2_demultiplexed/se/*fastq* ${projdir}/5_adapter_removed/se/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/se -type d -empty -delete
			else
				mv ${projdir}/3_motif_validated/se/*fastq* ${projdir}/5_adapter_removed/se/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/se -type d -empty -delete
			fi
		else
			mv ${projdir}/4_end_trimmed/se/*fastq* ${projdir}/5_adapter_removed/se/ 2> /dev/null && 	find ${projdir}/5_adapter_removed/se -type d -empty -delete
		fi
	fi

	if [[ "$rm_transit" == True ]]; then
		rm ${projdir}/2_demultiplexed/pe/*fastq* 2> /dev/null
		rm ${projdir}/2_demultiplexed/se/*fastq* 2> /dev/null
		rm ${projdir}/3_motif_validated/pe/*fastq* 2> /dev/null
		rm ${projdir}/3_motif_validated/se/*fastq* 2> /dev/null
		rm ${projdir}/4_end_trimmed/pe/*fastq* 2> /dev/null
		rm ${projdir}/4_end_trimmed/se/*fastq* 2> /dev/null
	fi
	find ${projdir}/5_adapter_removed -type d -empty -delete

	if [[ -d "${projdir}/5_adapter_removed/pe" ]]; then
		for fin in ${projdir}/5_adapter_removed/pe/*.R1.fastq.gz; do (
			python3 $krill -r1 ${fin} -r2 ${fin%.R1.fastq.gz}.R2.fastq.gz -q $q_min -p $q_percent -o ./pe/ &&
			fingz=${fin#*/pe/}; fingz=${fingz%.gz} &&
			cd pe &&
			$gzip ./*${fingz} && $gzip ./*${fingz%.R1.fastq}.R2.fastq &&
			cd ../ &&
			wait ) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
		wait
	fi
	if [[ -d "${projdir}/5_adapter_removed/se" ]]; then
		for fin in ${projdir}/5_adapter_removed/se/*R1.fastq.gz; do (
			python3 $krill -r1 ${fin} -q $q_min -p $q_percent -o ./se/ &&
			fingz=${fin#*/se/}; fingz=${fingz%.gz}
			cd se &&
			$gzip ./*${fingz} &&
			cd ../ &&
			wait

			if [[ -f ${fin%.R1.fastq.gz}.R2.fastq.gz ]]; then
				python3 $krill -r1 ${fin%.R1.fastq.gz}.R2.fastq.gz -q $q_min -p $q_percent -o ./se/ &&
				cd se &&
				$gzip ./*${fingz%.R1.fastq}.R2.fastq &&
				cd ../ &&
				wait
			fi
			) &
			if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
				wait
			fi
		done
		wait
	fi

	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	mv se pre_se
	mkdir -p se
	fin_se=$( ls ./pe/se.* 2> /dev/null | cat - <(ls ./pre_se/se.* 2> /dev/null ) | awk '{gsub(/ /,"\n"); gsub(/.\/pe\//,""); gsub(/.\/pre_se\//,""); gsub(/se./,"");}1' | sort | uniq)
	for i in $fin_se; do
		if [[ -f ./pe/se.${i} && -f ./pre_se/se.${i} ]]; then cat ./pe/se.${i} ./pre_se/se.${i} > ./se/$i; rm ./pe/se.${i}; fi &&
		if [[ -f ./pe/se.${i} && ! -f ./pre_se/se.${i} ]]; then mv ./pe/se.${i} ./se/$i; fi &&
		if [[ ! -f ./pe/se.${i} && -f ./pre_se/se.${i} ]]; then mv ./pre_se/se.${i} ./se/$i; fi &&
		wait
	done

	rm -rf pre_se
	cd ./pe
	if [[ ! -z "$(ls pe.* 2> /dev/null)" ]]; then
		for i in pe.*; do mv $i ${i#pe.}; done
	fi
	cd ../
	find . -type d -empty -delete
	find ./pe -size 0 -delete 2> /dev/null &&
	find ./se -size 0 -delete 2> /dev/null &&
	if [[ -d pe ]] && [[ "$(ls -A ./pe 2> /dev/null)" ]]; then for i in ./pe/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi
	if [[ -d se ]] && [[ "$(ls -A ./se 2> /dev/null)" ]]; then for i in ./se/*gz; do if [[ $(zcat $i | head -n 10 | wc -l ) -le 4 ]]; then  rm $i 2> /dev/null; fi; done && wait; fi


	if [[ "${QC_final}" =~ summary || "${QC_final}" =~ full ]]; then
		if [[ -d ${projdir}/6_quality_filtered_final/pe ]]; then
			cd ${projdir}/6_quality_filtered_final/pe
			mkdir -p qc
			for f in *.R1.fastq.gz; do (
				python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
				python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
				if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
					wait
				fi
			done
			wait

			cd ./qc
			if [[ ! "${QC_final}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_final}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1' > qscores_final_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1' > nucleotides_final_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R1_summary.csv nucleotides_final_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_final_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_final_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R2_summary.csv nucleotides_final_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages 2> /dev/null & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi

		if [[ -d ${projdir}/6_quality_filtered_final/se && "$test_lib_R2" != False ]]; then
			cd ${projdir}/6_quality_filtered_final/se
			mkdir -p qc
			if [[ "$test_lib_R2" != False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc &&
					python3 $crinoid -r1 ${f%.R1.fastq.gz}.R2.fastq.gz -t "$gthreads" -o ./qc 2> /dev/null ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi
			if [[ -d se && "$test_lib_R2" == False ]]; then
				for f in *.R1.fastq.gz; do (
					python3 $crinoid -r1 $f -t "$gthreads" -o ./qc ) &
					if [[ $(jobs -r -p | wc -l) -ge gN ]]; then
						wait
					fi
				done
				wait
			fi

			cd ./qc
			if [[ ! "${QC_final}" =~ full ]]; then
				rm *.png
			fi
			if [[ "${QC_final}" =~ summary ]]; then
				qscore_files=$(ls qscores*.R1.fastq.gz.csv)
				nucleotides_files=$(ls nucleotides*.R1.fastq.gz.csv)
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_final_R1_summary.csv &&
				awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
				awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_final_R1_summary.csv &&
				Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R1_summary.csv nucleotides_final_R1_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages & PIDR1=$!
				wait $PIDR1

				if [[ "$test_lib_R2" != False ]]; then
					qscore_files=$(ls qscores*.R2.fastq.gz.csv)
					nucleotides_files=$(ls nucleotides*.R2.fastq.gz.csv)
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $qscore_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > qscores_final_R2_summary.csv &&
					awk -F',' '{for (i=1;i<=NF;i++) total[FNR","i]+=$i;} END{for (j=1;j<=FNR;j++) {for (i=1;i<=NF;i++) printf "%3i ",total[j","i]; print "";}}' $nucleotides_files | \
					awk '{gsub(/  /,",");}1' | awk '{gsub(/ /,",");}1' | awk '{gsub(/,,/,",");}1' | awk '{gsub(/^/,""); gsub(/^/,"");}1'  > nucleotides_final_R2_summary.csv &&
					Rscript "${ngsComposer_dir}"/tools/helpers/qc_summary_plots.R qscores_final_R2_summary.csv nucleotides_final_R2_summary.csv ${ngsComposer_dir}/tools/helpers/R_packages 2> /dev/null & PIDR1=$!
					wait $PIDR1
				fi
			fi
			for i in nucleotides*.csv; do cat <(printf "A,C,G,T,N\n") $i > ${i}.tmp; mv ${i}.tmp $i; done
			for i in qscores*.csv; do for n in $(seq 1 $(awk -F',' '{print NF; exit}' $i)); do printf "q${n},"; done | awk '{gsub(/,$/,"");}1' - | cat - $i > ${i}.tmp; mv ${i}.tmp $i; done
			wait

			mkdir full summary
			mv *fastq* ./full/
			mv *_summary* ./summary/
			find . -type d -empty -delete
		fi
	fi
	echo "quality filtering complete" > ${projdir}/6_quality_filtered_complete

}
cd $projdir
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- Do you want to perform read quality-filtering? ${white}\n"
	read -p "- y(YES) or n(NO) " -t 36000 -n 1 -r
	if [[ ! $REPLY =~ ^[Yy]$ ]]; then
		printf '\n'
		echo -e "${magenta}- skipping read quality-filtering ${white}\n"
	else
		printf '\n'
		echo -e "${magenta}- performing read quality-filtering ${white}\n"
		time main_quality_filter &>> log.out
	fi
fi
if [ "$walkaway" == True ]; then
	if [ "$quality_filter" == 1 ]; then
		echo -e "${magenta}- performing read quality-filtering ${white}\n"
		time main_quality_filter &>> log.out
	else
		echo -e "${magenta}- skipping read quality-filtering ${white}\n"
	fi
fi
if [ "$walkaway" == False ]; then
	echo -e "${magenta}- quality filtering was performed with q_min="$q_min" q_percent="$q_percent" ${white}"
	echo -e "${magenta}- Check output and QC plot(s) ${white}"
	echo -e "${magenta}- Do you want to change q_min value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		q_min_new=$REPLY
	else
		q_min_new=$q_min
	fi
	echo -e "${magenta}- Do you want to change q_percent value? ${white}"
	read -p "- n(NO) or <enter_new_value>? " -n 1 -r
	if [[ ! $REPLY =~ ^[Nn]$ ]]; then
		q_percent_new=$REPLY
	else
		q_percent_new=$q_percent
	fi
	if [[ "$q_min" == "$q_min_new" && "$q_percent" == "$q_percent_new" ]]; then
		echo -e "${magenta}- ngsComposer will proceed to the next analytical step ${white}"
	else
		q_min=$q_min_new
		q_percent=$q_percent
		time main_quality_filter &>> log.out
	fi
fi
if [[ "$rm_transit" == True ]] && [[ -f "${projdir}/6_quality_filtered_complete" ]]; then
	rm ${projdir}/5_adapter_removed/pe/*fastq* 2> /dev/null
	rm ${projdir}/5_adapter_removed/se/*fastq* 2> /dev/null
fi


######################################################################################################################################################

cd ${projdir}
rm 1_initial_qc_complete 2_demultiplexed_complete 3_motif_removal_complete 4_end_trimming_complete 5_adapter_removal_complete 6_quality_filtered_complete
touch Analysis_Complete
echo -e "${magenta}- Run Complete. ${white}\n"
