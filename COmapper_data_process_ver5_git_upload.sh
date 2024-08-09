#!/bin/bash

#2023.09.07
#for processing Nanopore sequenicng data to COmapper input files
#updates
#ver.1 (23.08.31): merging command line codes, make input and output arguments
#ver.2 (23.09.04): make working directory, integrate argument names, add number of threads argument
#ver.3 (23.09.06): split read filtered sam files to parallel base filtering works, running test with test files
#ver.4 (23.09.07): move user defined argument and resource argument to front,

#------------------------------------------------------------------------------------------------------------------#

#User defined arguments
#folder location to make output files, ensure free space more than 2Tb (~200Gb fq.gz file standard)
workdir=/datasets/data_2/dohwan/Nanopore/240415_1000seedling/WT/raw
#number of threads
corenum=10
#name of dataset
name=WTseedling
#location of original fq.gz file
org=/datasets/data_2/dohwan/Nanopore/240415_1000seedling/WT/Col_Ler_F2_cleandata.fq.gz

#resource file location
#location of masksam.awk file
awkloc=/datasets/data_1/dohwan/Nanopore/masksam.awk
#location of reference genome mmi file
ref=/home/dohwan/refgenome/TAIR10/TAIR10.mmi

#------------------------------------------------------------------------------------------------------------------#

#making workspace
mkdir -p $workdir/nanoplot
mkdir -p ${workdir%raw}output
cd $workdir

#unzip
copy=$workdir/$name.fq
gzip -d $org
cp ${org%.gz} $copy
echo "unzip done"

#minimap2 alignment
minimap2 -t $corenum -ax map-ont --cs $ref $copy > ${copy%.fq}.sam
echo "alignment done"
#standard alignment, TAIR10 reference genome

#samtobam
sam=$workdir/$name.sam
bam=$workdir/$name.bam
samtools view -@ $corenum -Sb $sam > $bam
echo "bamconvert done"

#sort, qulaity filter
sorted=$workdir/$name.sort.bam
sbbout=$workdir/$name.sort.sbb.bam
samtools sort $bam  -@ $corenum -o $sorted -O bam
samtools index $sorted -@ $corenum
sambamba view -h -t $corenum -f bam \
		-F "not unmapped and not duplicate and mapping_quality >= 30 and sequence_length >= 100" \
		$sorted Chr1 Chr2 Chr3 Chr4 Chr5 > $sbbout
# retain
# 1. not unmapped
# 2. not duplicate
# 3. mapping quality >= 30
# 4. seqeucne length >= 100
echo "quality filter done"

#basefilter
samsbb=$workdir/$name.sort.sbb.sam
samtools view -@ $corenum $sbbout > $samsbb
cd $workdir
split -l 200000 -d $samsbb $name.sbbsplit. --additional-suffix=.sam

linenum=`wc -l $samsbb | awk '{print $1}'`
let num=$(($linenum/200000))
echo $num

i=1
for var in `seq -f "%02g" 0 $num`
do
    [ $((i%corenum)) == 0 ] && wait; i=$((i+1))
    sammp=$workdir/$name.sbbsplit.$var.sam
    filtermp=$workdir/$name.filtered.$var.sam
    $awkloc $sammp > $filtermp &
done
WORK_PID=`jobs -l | awk '{print $2}'`
wait $WORK_PID

echo "basefilter done"

i=1
for var in `seq -f "%02g" 0 $num`
do
    [ $((i%corenum)) == 0 ] && wait; i=$((i+1))
    filtermp=$workdir/$name.filtered.$var.sam
    outputmp=${workdir%raw}output/$name.input.$var.tsv
    echo -e "chr\tpos\tcigar\tseq" > $outputmp
    awk -F "\t" '{if($3!="ChrC" && $3!="ChrM" && $4!=0) print $3"\t"$4"\t"$6"\t"$10}' $filtermp >> $outputmp &
done
WORK_PID2=`jobs -l | awk '{print $2}'`
wait $WORK_PID2

echo "output file made"

#remove intermediate files (for free space)

#copied fq file
#rm $copy
#sam file
#rm $sam
#bam file
#rm $bam
#sorted bam, bai (index) file
#rm $sorted
#rm $sorted.bai
#read filtered bam, bai (index) file
#rm $sbbout
#rm $sbbout.bai
#read filtered sam file
#rm $samsbb
#splited read filtered sam files
#rm -r $name.sbbsplit.*.sam
#splited base filtered sam files
#rm -r $name.filtered.*.sam
