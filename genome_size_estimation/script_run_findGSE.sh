# USAGE:  # ./script_run_findGSE.sh <WORKTMP> <FASTQ> <ACC> <KMER> <CORES>

# Parse parameteres
WORKTMP=$1
FASTQ=$2
ACC=$3
KMER=$4
CORES=$5

# Declare other variables
REF=at_ChrC_ChrM_phiX74.fa
ref=$(echo "${REF}" | rev | cut -d'/' -f1 | cut -c 4- | rev)
parallel=$(( CORES / 4 ))

outputTMP=$WORKTMP/2_findGSE/$ACC

# Load conda environment
source activate jellyfish_v2.3.0
# Other required software
#bwa
#samtools
#bamToFastq

# Start
echo ""
echo "Start script"
date

# Create necessary directories
mkdir -p $outputTMP

# READ MAPPING and RETRIEVING UNMAPPED PAIRS

# This script assumes that RAW short-reads have already been trimmed 
echo "Aligning reads to Organelle genomes and PhiX..."
bwa mem -t $CORES -R "@RG\tID:$ACC\tSM:$ACC" $REF \
	$FASTQ.cutadapt.R1.fastq.gz \
	$FASTQ.cutadapt.R2.fastq.gz \
	| $samtools view -@ $CORES -Sbh - > $outputTMP/$ACC.$ref.uns.bam

echo "Sorting bam file..."
samtools sort -@ $CORES -n -o $outputTMP/$ACC.$ref.sort.bam $outputTMP/$ACC.$ref.uns.bam

echo "Estimating stats..."
samtools stats $outputTMP/$ACC.$ref.uns.bam |grep ^SN | cut -f 2- > $outputTMP/$ACC.$ref.stats.txt
samtools flagstat $outputTMP/$ACC.$ref.sort.bam > $outputTMP/$ACC.$ref.sort.txt

echo "Getting R1 unmapped, R2 mapped..."
samtools view -b -f 4 -F 264 $outputTMP/$ACC.$ref.sort.bam > $outputTMP/$ACC.$ref.unmap_map.bam
samtools flagstat $outputTMP/$ACC.$ref.unmap_map.bam > $outputTMP/$ACC.$ref.unmap_map.txt

echo "Getting R1 mapped, R2 unmapped..."
samtools view -b -f 8 -F 260 $outputTMP/$ACC.$ref.sort.bam > $outputTMP/$ACC.$ref.map_unmap.bam
samtools flagstat $outputTMP/$ACC.$ref.map_unmap.bam > $outputTMP/$ACC.$ref.map_unmap.txt

echo "Getting  R1 & R2 unmapped..."
samtools view -b -f 12 -F 256 $outputTMP/$ACC.$ref.sort.bam > $outputTMP/$ACC.$ref.unmap_unmap.bam
samtools flagstat $outputTMP/$ACC.$ref.unmap_unmap.bam > $outputTMP/$ACC.$ref.unmap_unmap.txt

echo "Merging both combinations of SE-mapped..."
rm $outputTMP/$ACC.$ref.TMP1_unmapped.bam
samtools merge $outputTMP/$ACC.$ref.TMP1_unmapped.bam $outputTMP/$ACC.$ref.unmap_map.bam $outputTMP/$ACC.$ref.map_unmap.bam

echo "Discarding supplementary alignments..."
samtools view -b -F 2048 $outputTMP/$ACC.$ref.TMP1_unmapped.bam > $outputTMP/$ACC.$ref.TMP2_unmapped.bam
samtools sort -n -o $outputTMP/$ACC.$ref.single_unmapped.bam $outputTMP/$ACC.$ref.TMP2_unmapped.bam

echo "Transforming BAM to FASTQ..."
bamToFastq -i $outputTMP/$ACC.$ref.single_unmapped.bam -fq $outputTMP/$ACC.single_unmapped.R1.fastq -fq2 $outputTMP/$ACC.single_unmapped.R2.fastq
bamToFastq -i $outputTMP/$ACC.$ref.unmap_unmap.bam -fq $outputTMP/$ACC.unmap_unmap.R1.fastq -fq2 $outputTMP/$ACC.unmap_unmap.R2.fastq

echo "Compressing FASTQ files..."
gzip $outputTMP/$ACC.*.R*

# Remove unnecesary intermediate files
rm $outputTMP/$ACC.$ref.uns.bam
rm $outputTMP/$ACC.$ref.unmapped.bam
rm $outputTMP/$ACC.$ref.TMP*
rm $outputTMP/$ACC.$ref.*map*.bam

# GENOME SIZE ESTIMATION

echo "Counting k-mers..."
zcat $outputTMP/$ACC.*.fastq.gz | jellyfish count /dev/fd/0 -C -o $outputTMP/$ACC."$KMER"mer -m $KMER -t $CORES -s 5G

echo "Generating histogram..."
rm $outputTMP/$ACC."$KMER"mer.histo
jellyfish histo -h 3000000 -o $outputTMP/$ACC."$KMER"mer.histo $outputTMP/$ACC."$KMER"mer

echo "Running findGSE..."
# Remove previous directory output
rm -rf $output/R_output_"$KMER"mer
Rscript --vanilla $PWD/findGSE.R $outputTMP/$ACC."$KMER"mer.histo $KMER $output/R_output_"$KMER"mer

# remove unnecessary files
rm $outputTMP/$ACC."$KMER"mer

date
echo "End script"


