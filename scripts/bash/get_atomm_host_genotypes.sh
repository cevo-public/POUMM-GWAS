# ==============================================================================
# Coverted PLINK BED file to haploid genotype matrix for ATOMM.
# ==============================================================================

# Run interactivelly with: 
# docker run -it --rm --mount type=bind,src=$PWD/scripts,dst=/scripts --mount type=bind,src=$PWD/output,dst=/output plink
# copy-paste these commands

DATA_PREFIX="spVL.shcs.rs.filtered"
MEM=10000

# Make sure REF allele is the major allele
./plink2 \
--bfile output/$DATA_PREFIX \
--maj-ref 'force' \
--make-bed \
--out output/$DATA_PREFIX.majref

# Get variant dosages from .bed file (binary)
./plink2 \
--bfile output/$DATA_PREFIX.majref \
--export A-transpose \
--memory $MEM \
--out output/$DATA_PREFIX.Av.majref

# Make genotypes haploid
head -1 output/$DATA_PREFIX.Av.majref.traw > output/atomm_sequence_host_header.txt

### TO REMOVE ###
# Down-sample host SNPs for first-pass analysis try
head -1000 output/$DATA_PREFIX.Av.majref.traw > output/$DATA_PREFIX.Av.majref.traw.subsample

# Reference allele is counted, so make missing data = 2 reference alleles
sed 's/\tNA/\t2/g' output/$DATA_PREFIX.Av.majref.traw.subsample| tail -n +2 > output/$DATA_PREFIX.Av.haptmp1.traw
# heterozygotes -> homozygous alternate allele (0)
sed 's/\t1/\t0/g' output/$DATA_PREFIX.Av.haptmp1.traw | tail -n +2 > output/$DATA_PREFIX.Av.haptmp2.traw
# homozygous reference allele -> 1
sed 's/\t2/\t1/g' output/$DATA_PREFIX.Av.haptmp2.traw | tail -n +2 > output/$DATA_PREFIX.Av.hap.traw
rm output/$DATA_PREFIX.Av.haptmp1.traw output/$DATA_PREFIX.Av.haptmp2.traw

# File looks like this now:
# CHR     SNP     (C)M    POS     COUNTED ALT     <host1>     <host2>     <host3> 
# 1       1:768116:A:AGTTTT       0       768116  AGTTTT  A       0       1       1 

# Format for ATOMM
# Make 2nd column (SNP IDs) numeric indices 
# Delete columns 3-6 (tr squishes the resulting field separators to a single space)
# Keep columns 1 and 7+ as-is
# This took ages, at least 6 hours, to run
i=1
while read -r line; do 
    echo $line | index="$i" awk '{$2=ENVIRON["index"]; $3=$4=$5=$6=""; print}' | tr -s ' '
    ((i++))
done < output/$DATA_PREFIX.Av.hap.traw > output/atomm_sequence_host.txt
