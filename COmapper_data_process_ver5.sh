#!/bin/bash

#2023.09.07
#for processing Nanopore sequenicng data to COmapper input files
#updates
#ver.1 (23.08.31): merging command line codes, make input and output arguments
#ver.2 (23.09.04): make working directory, integrate argument names, add number of threads argument
#ver.3 (23.09.06): split read filtered sam files to parallel base filtering works, running test with test files
#ver.4 (23.09.07): move user defined argument and resource argument to front,


#making workspace
mkdir -p ${dir_home}/output
cd $dir_home

#unzip
#copy=$workdir/$name.fq
gzip -d -c $dir_home/raw/$input > $dir_home/raw/${input%.gz}
#cp ${org%.gz} $copy
echo "unzip done"

input_prefix=${input%.fq.gz}

#minimap2 alignment
#minimap2 -t $n_parallel -ax map-ont --cs $ref $dir_home/raw/${input_prefix}.fq > $dir_home/output/${input_prefix}.sam
minimap2 -t $n_parallel -ax map-ont --cs $ref $dir_home/raw/${input_prefix}.fq | \
    samtools view -@ $n_parallel -Sb - > $dir_home/output/${input_prefix}.bam
echo "alignment done"
#standard alignment, TAIR10 reference genome


#sort, qulaity filter
samtools sort $dir_home/output/${input_prefix}.bam -@ $n_parallel -o $dir_home/output/${input_prefix}.sort.bam -O bam
samtools index $dir_home/output/${input_prefix}.sort.bam -@ $n_parallel
sambamba view -h -t $n_parallel -f bam \
		-F "not unmapped and not duplicate and mapping_quality >= 30 and sequence_length >= 100" \
		$dir_home/output/${input_prefix}.sort.bam Chr1 Chr2 Chr3 Chr4 Chr5 | \
        samtools view -@ $n_parallel - > ${dir_home}/output/${input_prefix}.sort.sbb.sam
# retain
# 1. not unmapped
# 2. not duplicate
# 3. mapping quality >= 30
# 4. seqeucne length >= 100
echo "quality filter done"

samsbb=${dir_home}/output/${input_prefix}.sort.sbb.sam
mkdir -p $dir_home/output/temp
mkdir -p $dir_home/output/tsv

# Split sam file for parallel processing
split -l 200000 -d $samsbb $dir_home/output/temp/${input_prefix}.sbbsplit. --additional-suffix=.sam

linenum=`wc -l $samsbb | awk '{print $1}'`
let num=$(($linenum/200000))
echo $num

i=1
for var in `seq -f "%02g" 0 $num`
do
    [ $((i%n_parallel)) == 0 ] && wait; i=$((i+1))
    samtmp=$dir_home/output/temp/${input_prefix}.sbbsplit.$var.sam
    filtertmp=$dir_home/output/temp/${input_prefix}.filtered.$var.sam
    $awkloc $samtmp > $filtertmp &
done
WORK_PID=`jobs -l | awk '{print $2}'`
wait $WORK_PID

echo "basefilter done"

i=1
for var in `seq -f "%02g" 0 $num`
do
    [ $((i%n_parallel)) == 0 ] && wait; i=$((i+1))
    filtermp=${dir_home}/output/temp/${input_prefix}.filtered.$var.sam
    outputmp=${dir_home}/output/tsv/${input_prefix}.input.$var.tsv
    echo -e "chr\tpos\tcigar\tseq" > $outputmp
    awk -F "\t" '{if($3!="ChrC" && $3!="ChrM" && $4!=0) print $3"\t"$4"\t"$6"\t"$10}' $filtermp >> $outputmp &
done
WORK_PID2=`jobs -l | awk '{print $2}'`
wait $WORK_PID2

echo "output file made"

#remove intermediate files (optional, for free space)

## copied fq file
#rm $copy
#sam file
#rm $sam
#bam file
#rm $bam
## sorted bam, bai (index) file
#rm $sorted
#rm $sorted.bai
## read filtered bam, bai (index) file
#rm $sbbout
#rm $sbbout.bai
## read filtered sam file
#rm $samsbb
## splited read filtered sam files
#rm -r $name.sbbsplit.*.sam
## splited base filtered sam files
#rm -r $name.filtered.*.sam
