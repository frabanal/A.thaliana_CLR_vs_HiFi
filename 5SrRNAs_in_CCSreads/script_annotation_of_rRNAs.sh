# USAGE:   # ./script_annotation_of_rRNAs.sh <WORKDIR> <ASSEMBLY> <ACC> <CORES>

# Parse parameteres
WORKTMP=$1
REF=$2
ACC=$3
CORES=$4

# Declare other variables
parallel=$(( CORES / 4 ))
LIB=$PWD/rRNAs.fa
outputTMP=$WORKTMP/$ACC

# Other required software
SAMTOOLS=<PATH/TO/samtools>


# Create necessary directories
mkdir -p $outputTMP

cd $outputTMP

# Start
echo ""
echo "START script"
date

################################
######### REPEATS ##############
################################

echo "Preparing Genome..."
# Copy reference genome to output directory (soft links became problematic)
cp $REF $outputTMP/$ACC.fa
chmod +w $outputTMP/$ACC.fa

# Make sure to change all lower-case to UPPER-case
#sed -i -e '/>/b; s/a/A/g; s/c/C/g; s/g/G/g; s/t/T/g' $outputTMP/$ACC.fa

# Index reference for future usage
$SAMTOOLS faidx $outputTMP/$ACC.fa

# Activate virtual environment
source activate /ebio/abt6_projects8/alopecurus_genome/bin/anaconda2/envs/RepeatMasker_v4-0-9

echo "Running RepeatMasker – 1st round..."
RepeatMasker -lib $LIB -pa $parallel -cutoff 200 -nolow -gff -xsmall $outputTMP/$ACC.fa

echo "Re-formatting gff2 to gff3..."

sed "s/\"//g" $outputTMP/$ACC.fa.out.gff > $outputTMP/$ACC.gff2

sed -i "s/Motif:/Repeat_type=/g" $outputTMP/$ACC.gff2

## For rDNA (DNA including ITS and ETS)
grep "=5S_rDNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5S_rDNA", $4, $5, ".", $7, $8, "ID=5S_rDNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5S_rDNA.gff
grep "=45S_rDNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "45S_rDNA", $4, $5, ".", $7, $8, "ID=45S_rDNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.45S_rDNA.gff

## For rRNA (RNA subunits)
grep "=5S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5S_rRNA", $4, $5, ".", $7, $8, "ID=5S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5S_rRNA.gff 
grep "=18S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "18S_rRNA", $4, $5, ".", $7, $8, "ID=18S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.18S_rRNA.gff
grep "=5.8S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5.8S_rRNA", $4, $5, ".", $7, $8, "ID=5.8S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5.8S_rRNA.gff
grep "=25S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "25S_rRNA", $4, $5, ".", $7, $8, "ID=25S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.25S_rRNA.gff

conda deactivate

# Concatenate GFF files from RepeatMasker and SORT them
cat $outputTMP/$ACC.*_rRNA.gff \
	> $outputTMP/$ACC.rRNA.TMP.gff

sort -k1,1 -k4n $outputTMP/$ACC.rRNA.TMP.gff > $outputTMP/$ACC.rRNA.gff

sed -i '1i ##gff-version 3' $outputTMP/$ACC.rRNA.gff
sed -i "2i ##date $(date +%Y-%m-%d)" $outputTMP/$ACC.rRNA.gff


# Concatenate GFF files from RepeatMasker and SORT them
cat $outputTMP/$ACC.*_rDNA.gff \
	> $outputTMP/$ACC.rDNA.TMP.gff

sort -k1,1 -k4n $outputTMP/$ACC.rDNA.TMP.gff > $outputTMP/$ACC.rDNA.gff

sed -i '1i ##gff-version 3' $outputTMP/$ACC.rDNA.gff
sed -i "2i ##date $(date +%Y-%m-%d)" $outputTMP/$ACC.rDNA.gff

date
echo "END script"
echo ""
