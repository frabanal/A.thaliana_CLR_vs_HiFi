# USAGE:  # ./script_combine_EDTA_and_REPEATS.sh <WORKTMP> <input_EDTA> <input_REPEATS> <ACC>

WORKTMP=$1
input_EDTA=$2
input_REPEATS=$3
ACC=$4

# Declare other variables
inputDIR=$(echo "${input_EDTA}" | rev | cut --complement -d'/' -f1 | rev)

outputTMP=$WORKTMP/$ACC

# Create necessary directories
mkdir -p $outputTMP

cd $outputTMP

# Start
echo ""
echo "START script"
date

########################################
############ Simplify TEs ##############
########################################

echo "Removing signatures of LTRs..."
grep -v "ID=repeat_" $input_EDTA > $inputDIR/$ACC.EDTA.TMP1.gff3 
grep -v "long_terminal_repeat" $inputDIR/$ACC.EDTA.TMP1.gff3 > $inputDIR/$ACC.EDTA.TMP2.gff3  
grep -v "target_site_duplication" $inputDIR/$ACC.EDTA.TMP2.gff3 > $inputDIR/$ACC.EDTA.TEanno.edit.gff3

echo "Classifying TEs..."
sed -i 's/Copia_LTR_retrotransposon/ClassI_LTR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Gypsy_LTR_retrotransposon/ClassI_LTR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/LTR_retrotransposon/ClassI_LTR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/LINE_element/ClassI_nonLTR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/CACTA_TIR_transposon/ClassII_TIR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/hAT_TIR_transposon/ClassII_TIR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Mutator_TIR_transposon/ClassII_TIR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/PIF_Harbinger_TIR_transposon/ClassII_TIR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3
sed -i 's/Tc1_Mariner_TIR_transposon/ClassII_TIR/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3

sed -i 's/helitron/ClassII_Helitrons/' $inputDIR/$ACC.EDTA.TEanno.edit.gff3

rm $inputDIR/$ACC.EDTA.TMP*


##########################################
############ EDTA + REPEATS ##############
##########################################

# Required Software
BEDTOOLS=<PATH/TO/bedtools>

cd $output

echo "# Pre-processing EDTA.gff to filter out 'TEs' overlapping rRNAs..."

# "Only report those entries in A that have no overlap in B. Restricted by -f and -r."
$BEDTOOLS intersect -v -a $inputDIR/$ACC.EDTA.TEanno.edit.gff3 -b $input_REPEATS > EDTA_rDNAPurged.v.gff

# "Perform a “left outer join”. That is, for each feature in A report each overlap with B. If no overlaps are found, report a NULL feature for B." 
$BEDTOOLS intersect -loj -a $inputDIR/$ACC.EDTA.TEanno.edit.gff3 -b $input_REPEATS > EDTA_rDNAPurged.loj.gff

# Concatenate Repeats and TE annotations (only those that do not overlap) 
cat $input_REPEATS EDTA_rDNAPurged.v.gff > $ACC.Repeats_TEanno.TMP1.gff3

# Remove header and sort
grep -v "^#"  $ACC.Repeats_TEanno.TMP1.gff3 > $ACC.Repeats_TEanno.TMP2.gff3
sort -k1,1 -k4n $ACC.Repeats_TEanno.TMP2.gff3 > $ACC.Repeats_TEanno.gff3

sed -i '1i ##gff-version 3' $ACC.Repeats_TEanno.gff3
sed -i "2i ##date $(date +%Y-%m-%d)" $ACC.Repeats_TEanno.gff3

# BEFORE and AFTER TE count
echo "TEs BEFORE purging:"
grep -c -v "^#" $inputDIR/$ACC.EDTA.TEanno.edit.gff3
echo "TEs AFTER purging:"
grep -c -v "^#" EDTA_rDNAPurged.v.gff

rm $ACC.Repeats_TEanno.TMP*

date
echo "END script"
echo ""
