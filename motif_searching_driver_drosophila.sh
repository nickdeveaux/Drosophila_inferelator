# Written by Alex Rives, modified by Nick De Veaux

# This script computes a comprehensive motif hit dataset from the reference human motif dataset and an input peak bed
# Inputs: A merged bed file of putative peaks: typically, the merged bedfiles across all samples
# using peakdeck to identify peaks, peakMaxArray.py to find summits, and bedtools slop +/- 25 bp
#!/usr/bin/sh

# To run ./motif_searching_driver_drosophila.sh /mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster/sorted_slopped_genes_1.kb.bed 1000
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/meme/4.10.1/bin:$PATH
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/bedtools/bedtools-2.25.0:$PATH
export PATH=/mnt/xfs1/bioinfoCentos7/software/installs/python/anaconda/bin:$PATH

peaksbed=/mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster/1kb_TSS_nomito_sorted_filtered_short_chroms.bed
# Alternatively, use the 1kb:
# /mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster/introns_final.bed
distance=0

genome="/mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster/dm6.fa"
window_file=$peaksbed
motif_database="/mnt/ceph/users/ndeveaux/reference/drosophila_melanogaster/fly_factor_survey.meme"
# The corresponding metadata file is TF_Information_hg19_em.txt
output_dir="/mnt/ceph/users/ndeveaux/DMel_Inferelator/flyfactorsurvey_motif_finding_output"
bgOrder=1

mkdir -p $output_dir
outDir=$output_dir/$(basename $peaksbed)
mkdir -p $outDir

cores=$(nproc)

function log {
    echo $(date +%F_%T) $$ ${BASHPID} $1
}

# Generate a temp working directory
tmpDirRoot=/tmp
wd=$(mktemp -d ${tmpDirRoot}/${USER}_MOTIF_XXXXXXXXXX)
cd $wd

log "Working in $(pwd) on $(hostname)"

# Create ATAC-sample specific background using the peaksbed input
# TODO: allow two bed inputs: one for all ATAC samples to create background, the other for motif searching

## generate background files if necessary
sortedBkgrdBed=$peaksbed

# 2.  see whether a background fasta file exists
backgroundFasta="${peaksbed/\.bed/.fa}"
if [ ! -s ${backgroundFasta} ]
then
        fastaFromBed -fi \
                ${genome} \
                -bed ${sortedBkgrdBed} \
                -fo ${backgroundFasta}
        log "${backgroundFasta} generated."
fi

# 3. generate markov models of orders specified above
backgroundFile=${peaksbed/\.bed/}_bkgrd_Order${bgOrder}.txt
# test to see whether the gene bed file exists, and, if it doesn't, make it:
if [ ! -s ${backgroundFile} ]
then
        fasta-get-markov -m ${bgOrder} ${backgroundFasta} ${backgroundFile}
        log "${backgroundFile} generated."
fi
cp ${backgroundFile} ${wd}
cp ${backgroundFasta} ${wd}
localBkgFile=${wd}/$(basename ${backgroundFile})
localFastaFile=${wd}/$(basename ${backgroundFasta})


# get a list of motifs from the meme database
grep "MOTIF " $motif_database | awk '{print $2}' > "$wd/motifs.txt"

mkdir "$wd/motifs"
mkdir "$wd/jobs"

fimo=`which fimo`
# search for motifs within the peaks
for motif in $(cat "$wd/motifs.txt"); do

  # Get motif hits:
  echo "$fimo \
                --parse-genomic-coord \
                --text \
                --thresh .001 \
                --bgfile ${localBkgFile} \
                --verbosity 1 \
                --motif ${motif} $motif_database \
                $localFastaFile > $wd/motifs/${motif// /_}.txt \
                2> $wd/motifs/${motif// /_}.err"  >> "$wd/jobs/motifs.sh";
done

cat  "$wd/jobs/motifs.sh" | xargs -I CMD -P $cores bash -c CMD

log "Finished running fimo search, outputting txt files"
# convert motif outputs to bedfiles
for f in $wd/motifs/*.txt; do
  motif=$(basename $f .txt)
  bed=$motif.bed
  awk -F'\t' 'NR>1 {print $2"\t"$3"\t"$4"\t"$1"\t"$7"\t"$6}' $f > $wd/motifs/$bed
done

bedtools=`which bedtools` 

log "Finished converting fimo txt outputs to bed files"

# find the target genes for each motif
mkdir "$wd/targets"
rm "$wd/jobs/targets.sh"
for f in $wd/motifs/*.bed; do
  if [[ -s $f ]]; then
    motif=$(basename $f .bed)
    targetbed=${motif}_targets.bed
    # Assign motifs to genes
    echo "$bedtools window -w $distance -a $f -b $window_file > $wd/targets/$targetbed" >> "$wd/jobs/targets.sh"
  fi
done

cat $wd/jobs/targets.sh | xargs -I CMD -P $cores bash -c CMD 

# python $priors_script

mkdir -p $outDir
rsync -av * $outDir


