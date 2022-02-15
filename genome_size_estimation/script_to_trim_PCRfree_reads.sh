# USAGE:  # ./script_to_trim_PCRfree_reads.sh <WORKTMP> <SRA> <FASTQ> <SAMPLE> <CORES> 

# Input variables
WORKTMP=$1
SRA=$2
FASTQ=$3
SAMPLE=$4
CORES=$5

# Other input variables
length=75

inputTMP=$WORKTMP/0_fastq
outputTMP=$WORKTMP/1_fastqc

fastq1="$FASTQ"_R1_001
fastq2="$FASTQ"_R2_001

# Load virtual environment for Cutadapt
source  activate cutadapt_python_v3.6.8

# Create necessary directories
mkdir -p $inputTMP
mkdir -p $outputTMP

# Start pipeline
echo ""
date

echo "Creating soft links..."
ln -s $SRA/$fastq1.fastq.gz $inputTMP/$SAMPLE.R1.fastq.gz
ln -s $SRA/$fastq2.fastq.gz $inputTMP/$SAMPLE.R2.fastq.gz

echo "Trimming adapters..." 
cutadapt -j $CORES -q 20,15 \
	-b TruSeq1=AGATCGGAAGAGC -b Nextera1=CTGTCTCTTATACACATCT -b Nextera1rc=AGATGTGTATAAGAGACAG \
	-B TruSeq2=AGATCGGAAGAGC -B Nextera2=CTGTCTCTTATACACATCT -B Nextera2rc=AGATGTGTATAAGAGACAG \
	--trim-n \
	--minimum-length $length \
	-o $inputTMP/$SAMPLE.cutadapt.R1.fastq.gz --paired-output $inputTMP/$SAMPLE.cutadapt.R2.fastq.gz \
	$inputTMP/$SAMPLE.R1.fastq.gz $inputTMP/$SAMPLE.R2.fastq.gz

# Deactivate virtual environment for Cutadapt
conda deactivate

echo "DONE!"
date
echo ""

