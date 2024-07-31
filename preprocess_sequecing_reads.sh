#!/bin/bash
#conda activate RNA
DATA_PATH="Metadata/"
OUT_PATH="Processed"
NATIVE_PATH="native_data"
RNA="5TPY"
# Define forward and reverse read files for pre- and post-selection
pre_forward_reads_file="${RNA}_pre_1.fastq"
pre_reverse_reads_file="${RNA}_pre_2.fastq"
post_forward_reads_file="${RNA}_post_1.fastq"
post_reverse_reads_file="${RNA}_post_2.fastq"
target_pre="${RNA}_pre"
target_post="${RNA}_post"
# Setting your adapter seqs.
seq_adapter3="CACTGGCCGTCGTTTTAC"
seq_adapter5="GTTTTCCCAGTCACGAC"
min_len=20  #Cutoff of sequence length
N_cpu=20    #Number of CPU 
fastx_path="/Users/s3000762/Documents/Software/fastx/fastx_toolkit/src"
fastq_quality_filter="$fastx_path/fastq_quality_filter/fastq_quality_filter"
fastx_collapser="$fastx_path/fastx_collapser/fastx_collapser"

min_quality=20   ## Quality threshold value    （parameter In fastq_quality_filter）
min_prob=90      ## The proportion of high quality bases
if [ ! -d "$OUT_PATH" ]; then
    mkdir -p "$OUT_PATH"
    echo "Directory $OUT_PATH created."
fi
# Function to process sequencing data
process_data () {
    local forward_reads_file=$1
    local reverse_reads_file=$2
    local target=$3
    ###Zip the metadata, if the data is .gz format (e.g, .fastq.gz), else continue
    #gzip -dk ${forward_reads_file}.gz
    #gzip -dk ${reverse_reads_file}.gz
    echo "######## Assembles Illumina paired-end reads using PEAR ########"
    pear -f ${DATA_PATH}/${forward_reads_file} -r ${DATA_PATH}/${reverse_reads_file} -o ${OUT_PATH}/${target} -j $N_cpu -y 1G
    rm ${OUT_PATH}/${target}.unassembled.forward.fastq
    rm ${OUT_PATH}/${target}.unassembled.reverse.fastq
    rm ${OUT_PATH}/${target}.discarded.fastq

    echo "######## Removes 3' adapter sequences ########"
    cutadapt -m $min_len -j $N_cpu -a $seq_adapter3 -o ${OUT_PATH}/${target}_tmp.fq ${OUT_PATH}/${target}.assembled.fastq
    echo "######## removes 5' adapter sequences ########"
    cutadapt -m $min_len -j $N_cpu -g $seq_adapter5 -o ${OUT_PATH}/${target}_cut.fq ${OUT_PATH}/${target}_tmp.fq
    rm ${OUT_PATH}/${target}_tmp.fq

    echo "################################"
    echo "Filter low quality sequences with less than ${min_prob}% nucleotides no more than ${min_quality}"
    echo "################################"
    $fastq_quality_filter -q $min_quality -p $min_prob -i ${OUT_PATH}/${target}_cut.fq -o ${OUT_PATH}/${target}_filtered.fq
    echo "Combine the repeat sequences & perform frequency statistics"
    $fastx_collapser -i ${OUT_PATH}/${target}_filtered.fq -o ${OUT_PATH}/${target}_nt.count

    rm ${OUT_PATH}/${target}.assembled.fastq ${OUT_PATH}/${target}_cut.fq ${OUT_PATH}/${target}_filtered.fq
}

# Process pre-selection data
echo "######## Process pre-selection data ########"
process_data $pre_forward_reads_file $pre_reverse_reads_file $target_pre

# Process post-selection data
echo "######## Process post-selection data ########"
process_data $post_forward_reads_file $post_reverse_reads_file $target_post
echo "#####################################"
echo "Great! Sequencing data processing is complete and the generated .count files (in ${OUT_PATH}) can be used for subsequent fitness calculations."
echo "#####################################"

echo "#####################################"
echo "Fitness calculation -->"  
echo "#####################################"
sed 's/U/T/g' "${NATIVE_PATH}/${RNA}.fasta" > "${NATIVE_PATH}/${RNA}-DNA.fasta"
python getRA.py -n ${NATIVE_PATH}/${RNA}-DNA.fasta -b ${OUT_PATH}/${target_pre}_nt.count -a ${OUT_PATH}/${target_post}_nt.count -o ${OUT_PATH}/${RNA}
echo "#####################################"
echo "Good job! file of "${RNA}.var.ra" in ${OUT_PATH} is all your need to infer the base-pairing."  
echo "#####################################"