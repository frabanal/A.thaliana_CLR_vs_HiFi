# USAGE:  # ./script_annotation_of_repeats.sh <WORKDIR> <LIB> <ASSEMBLY> <ACC> <CORES>

# Parse parameteres
WORKTMP=$1
LIB=$2
ASSEMBLY=$3
ACC=$4
CORES=$5

# Declare other variables
parallel=$(( CORES / 4 ))

ChrM=at_mitochondria.fa
ChrC=at_chloroplast.fa

outputTMP=$WORKTMP/$ACC

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
# Copy reference genome to output directory
cp $ASSEMBLY $outputTMP/$ACC.fa
chmod +w $outputTMP/$ACC.fa

# Make sure to change all lower-case to UPPER-case
sed -i -e '/>/b; s/a/A/g; s/c/C/g; s/g/G/g; s/t/T/g' $outputTMP/$ACC.fa

# Index reference for future usage
# Required Software
SAMTOOLS=<PATH/TO/samtools>
$SAMTOOLS faidx $outputTMP/$ACC.fa

# Activate virtual environment
source activate RepeatMasker_v4-0-9

echo "Running RepeatMasker..."
RepeatMasker -lib $LIB -pa $parallel -cutoff 200 -nolow -gff -xsmall $outputTMP/$ACC.fa

echo "Re-formatting gff2 to gff3..."
sed "s/\"//g" $outputTMP/$ACC.fa.out.gff > $outputTMP/$ACC.gff2
sed -i "s/Motif:/Repeat_type=/g" $outputTMP/$ACC.gff2

# For rDNA (DNA including ITS and ETS)
grep "=5S_rDNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5S_rDNA", $4, $5, ".", $7, $8, "ID=5S_rDNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5S_rDNA.gff
grep "=45S_rDNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "45S_rDNA", $4, $5, ".", $7, $8, "ID=45S_rDNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.45S_rDNA.gff

# For rRNA (RNA subunits)
grep "=5S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5S_rRNA", $4, $5, ".", $7, $8, "ID=5S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5S_rRNA.gff 
grep "=18S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "18S_rRNA", $4, $5, ".", $7, $8, "ID=18S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.18S_rRNA.gff
grep "=5.8S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "5.8S_rRNA", $4, $5, ".", $7, $8, "ID=5.8S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.5.8S_rRNA.gff
grep "=25S_rRNA" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "25S_rRNA", $4, $5, ".", $7, $8, "ID=25S_rRNA_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.25S_rRNA.gff

# Other repeats
grep "=centromere" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "centromere", $4, $5, ".", $7, $8, "ID=centromere_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.centromere.gff
grep "=telomere" $outputTMP/$ACC.gff2 | awk -v OFS='\t' '{print $1, $2, "telomere", $4, $5, ".", $7, $8, "ID=telomere_"NR";"$10";Length="$12-$11";Divergence="$6""  }' > $outputTMP/$ACC.telomere.gff

conda deactivate

################################
######### N-stretches ##########
################################

# This section requires the seqtk software
SEQTK=<PATH/TO/seqtk>


echo "Identifying N-stretches..."
$SEQTK cutN -n 10 -g -p 1000000000000 $outputTMP/$ACC.fa > $outputTMP/$ACC.TMP.txt
awk -F"\t" -v OFS="\t" '{$4=$3-$2}1' $outputTMP/$ACC.TMP.txt > $outputTMP/$ACC.N.txt

echo "Convert previously generated N-stretches file into GFF..."
awk -v OFS='\t' '{print $1, "seqtk", "N_stretch", $2, $3, ".", "+", ".", "ID=N_stretch_"NR";Length="$4""  }' $outputTMP/$ACC.N.txt > $outputTMP/$ACC.N.gff

################################
######### ORGANELLES ###########
################################

# Activate virtual environment
source activate pbbioconda_SMRTlinkv9.0

# Align ORGANELLES against the assembly
echo "Aligning mitochondria against scaffolds/contigs..."
minimap2 -t $CORES -cx asm5 $ChrM $outputTMP/$ACC.fa > $outputTMP/$ACC.ChrM.paf

echo "Aligning mitochondria against scaffolds/contigs..."
minimap2 -t $CORES -cx asm5 $ChrC $outputTMP/$ACC.fa > $outputTMP/$ACC.ChrC.paf

echo "Formating PAF into GFF files..."

awk -v OFS='\t' '{print $1, "minimap2", $6, $3+=1, $4, ".", $5, ".", "ID="$6"_"NR";Repeat_type="$6";Length="$11";Coordinates="$8"-"$9";Identity="$10/$11"" }' $outputTMP/$ACC.ChrM.paf > $outputTMP/$ACC.ChrM.gff
awk -v OFS='\t' '{print $1, "minimap2", $6, $3+=1, $4, ".", $5, ".", "ID="$6"_"NR";Repeat_type="$6";Length="$11";Coordinates="$8"-"$9";Identity="$10/$11"" }' $outputTMP/$ACC.ChrC.paf > $outputTMP/$ACC.ChrC.gff

# Concatenate GFF files from RepeatMasker, N_stretches and minimap2
cat $outputTMP/$ACC.*_rDNA.gff \
	$outputTMP/$ACC.*_rRNA.gff \
	$outputTMP/$ACC.centromere.gff \
	$outputTMP/$ACC.telomere.gff \
	$outputTMP/$ACC.N.gff \
	$outputTMP/$ACC.ChrM.gff \
	$outputTMP/$ACC.ChrC.gff \
	> $outputTMP/$ACC.Repeats.TMP.gff

# Sort GFF files
sort -k1,1 -k4n $outputTMP/$ACC.Repeats.TMP.gff > $outputTMP/$ACC.Repeats.gff

sed -i '1i ##gff-version 3' $outputTMP/$ACC.Repeats.gff
sed -i "2i ##date $(date +%Y-%m-%d)" $outputTMP/$ACC.Repeats.gff

conda deactivate

echo "Removing unnecessary files..."
rm $outputTMP/$ACC.fa.out.gff
rm $outputTMP/$ACC.fa.masked
rm $outputTMP/$ACC.TMP.txt
rm $outputTMP/$ACC.Repeats.TMP.gff


######################################################
################## Fix OVERLAPS ######################
######################################################

# Required Software
BEDTOOLS=<PATH/TO/bedtools>

echo "Merging overlapping annotations – Organelles..."
#ChrM
sort -k1,1 -k4n $outputTMP/$ACC.ChrC.gff > $outputTMP/$ACC.ChrC.sorted.gff
$BEDTOOLS merge -i $outputTMP/$ACC.ChrC.sorted.gff -c 1 -o count > $outputTMP/$ACC.ChrC.out
awk -v OFS='\t' '{print $1, "minimap2", "chloroplast", $2+=1, $3, ".", "+", ".", "ID=chloroplast_"NR";Repeat_type=chloroplast;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.ChrC.out > $outputTMP/$ACC.ChrC.merged.gff

#ChrC
sort -k1,1 -k4n $outputTMP/$ACC.ChrM.gff > $outputTMP/$ACC.ChrM.sorted.gff
$BEDTOOLS merge -i $outputTMP/$ACC.ChrM.sorted.gff -c 1 -o count > $outputTMP/$ACC.ChrM.out
awk -v OFS='\t' '{print $1, "minimap2", "mitochondria", $2+=1, $3, ".", "+", ".", "ID=mitochondria_"NR";Repeat_type=mitochondria;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.ChrM.out > $outputTMP/$ACC.ChrM.merged.gff

echo "Combine intersected and unique chloroplast and mitochondria annotations..."
#Intersect
$BEDTOOLS intersect -a  $outputTMP/$ACC.ChrC.merged.gff -b  $outputTMP/$ACC.ChrM.merged.gff -wa -wb > $outputTMP/$ACC.ChrCM.wab.gff
#Pick largest and remove duplicate entries
awk -v OFS='\t' '{if ($5-$4 >= $14-$13) print $1, $2, $3, $4, $5, $6, $7, $8, $9; else print $10, $11, $12, $13, $14, $15, $16, $17, $18;}' $outputTMP/$ACC.ChrCM.wab.gff > $outputTMP/$ACC.ChrCM.larger.gff
uniq $outputTMP/$ACC.ChrCM.larger.gff >  $outputTMP/$ACC.ChrCM.uniq.gff
#Obtain unique annotation
$BEDTOOLS intersect -a  $outputTMP/$ACC.ChrC.merged.gff -b  $outputTMP/$ACC.ChrM.merged.gff -v > $outputTMP/$ACC.ChrC.v.gff
$BEDTOOLS intersect -a  $outputTMP/$ACC.ChrM.merged.gff -b  $outputTMP/$ACC.ChrC.merged.gff -v > $outputTMP/$ACC.ChrM.v.gff
#Concatenate all organellar hits
cat $outputTMP/$ACC.ChrC.v.gff $outputTMP/$ACC.ChrM.v.gff $outputTMP/$ACC.ChrCM.uniq.gff \
	| sort -k1,1 -k4n \
	> $outputTMP/$ACC.ChrCM.merged.gff

echo "Merging overlapping annotations – rDNA..."
# 45S
$BEDTOOLS merge -i $outputTMP/$ACC.45S_rDNA.gff -c 1 -o count > $outputTMP/$ACC.45S_rDNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "45S_rDNA", $2+=1, $3, ".", "+", ".", "ID=45S_rDNA_"NR";Repeat_type=45S_rDNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.45S_rDNA.out > $outputTMP/$ACC.45S_rDNA.merged.gff
# 5S
$BEDTOOLS merge -i $outputTMP/$ACC.5S_rDNA.gff -c 1 -o count > $outputTMP/$ACC.5S_rDNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "5S_rDNA", $2+=1, $3, ".", "+", ".", "ID=5S_rDNA_"NR";Repeat_type=5S_rDNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.5S_rDNA.out > $outputTMP/$ACC.5S_rDNA.merged.gff

echo "Merging overlapping annotations – rRNAs..."
# 18S
$BEDTOOLS merge -i $outputTMP/$ACC.18S_rRNA.gff -c 1 -o count > $outputTMP/$ACC.18S_rRNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "18S_rRNA", $2+=1, $3, ".", "+", ".", "ID=18S_rRNA_"NR";Repeat_type=18S_rRNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.18S_rRNA.out > $outputTMP/$ACC.18S_rRNA.merged.gff
# 5.8S
$BEDTOOLS merge -i $outputTMP/$ACC.5.8S_rRNA.gff -c 1 -o count > $outputTMP/$ACC.5.8S_rRNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "5.8S_rRNA", $2+=1, $3, ".", "+", ".", "ID=5.8S_rRNA_"NR";Repeat_type=5.8S_rRNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.5.8S_rRNA.out > $outputTMP/$ACC.5.8S_rRNA.merged.gff
# 25S
$BEDTOOLS merge -i $outputTMP/$ACC.25S_rRNA.gff -c 1 -o count > $outputTMP/$ACC.25S_rRNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "25S_rRNA", $2+=1, $3, ".", "+", ".", "ID=25S_rRNA_"NR";Repeat_type=25S_rRNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.25S_rRNA.out > $outputTMP/$ACC.25S_rRNA.merged.gff
# 5S
$BEDTOOLS merge -i $outputTMP/$ACC.5S_rRNA.gff -c 1 -o count > $outputTMP/$ACC.5S_rRNA.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "5S_rRNA", $2+=1, $3, ".", "+", ".", "ID=5S_rRNA_"NR";Repeat_type=5S_rRNA;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.5S_rRNA.out > $outputTMP/$ACC.5S_rRNA.merged.gff

echo "Filtering out organellar ribosomes..."
# Pre-Concatenate ALL rDNA/rRNA annotations prior to overlaps check with mitochondria 
cat $outputTMP/$ACC.*_rDNA.merged.gff \
	$outputTMP/$ACC.*_rRNA.merged.gff \
	> $outputTMP/$ACC.ribosomes.merged.gff

# "Only report those entries in A that have no overlap in B. Restricted by -f and -r."
$BEDTOOLS intersect -a $outputTMP/$ACC.ribosomes.merged.gff -b $outputTMP/$ACC.ChrCM.merged.gff -v > $outputTMP/$ACC.ribosomes.v.gff

echo "Merging overlapping annotations – centromeres..."
# centromere
$BEDTOOLS merge -i $outputTMP/$ACC.centromere.gff -c 1 -o count > $outputTMP/$ACC.centromere.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "centromere", $2+=1, $3, ".", "+", ".", "ID=centromere_"NR";Repeat_type=centromere;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.centromere.out > $outputTMP/$ACC.centromere.merged.gff

echo "Merging overlapping annotations – telomeres..."
# telomere
$BEDTOOLS merge -i $outputTMP/$ACC.telomere.gff -c 1 -o count > $outputTMP/$ACC.telomere.out
awk -v OFS='\t' '{print $1, "RepeatMasker", "telomere", $2+=1, $3, ".", "+", ".", "ID=telomere_"NR";Repeat_type=telomere;Length="$3-$2";Merged="$4"" }' $outputTMP/$ACC.telomere.out > $outputTMP/$ACC.telomere.merged.gff

echo "Concatenating GFF files..."
## Concatenate GFF files after merging steps
cat $outputTMP/$ACC.ribosomes.v.gff \
	$outputTMP/$ACC.centromere.merged.gff \
	$outputTMP/$ACC.telomere.merged.gff \
	$outputTMP/$ACC.N.gff \
	$outputTMP/$ACC.ChrCM.merged.gff \
	> $outputTMP/$ACC.Repeats_merged.TMP.gff

# Sort GFF files
sort -k1,1 -k4n $outputTMP/$ACC.Repeats_merged.TMP.gff > $outputTMP/$ACC.Repeats_merged.gff
sed -i '1i ##gff-version 3' $outputTMP/$ACC.Repeats_merged.gff
sed -i "2i ##date $(date +%Y-%m-%d)" $outputTMP/$ACC.Repeats_merged.gff

echo "Removing unneccesary files..."
rm $outputTMP/$ACC.*.out
rm $outputTMP/$ACC.ChrC.*.gff
rm $outputTMP/$ACC.ChrM.*.gff
rm $outputTMP/$ACC.Repeats_merged.TMP.gff


date
echo "END script"
echo ""
