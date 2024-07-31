![MobiSeq_logo](https://github.com/user-attachments/assets/e148b59f-b120-4f98-8443-ee50385f742f)
# RNA-MobiSeq
Deep mutational scanning and mobility-based selection for RNA base-pair and tertiary structure inference
## Introduction
-------------
The RNA-MobiSeq pipeline is a method to predict the 3D structure of a target RNA from its sequence by combining deep mutation, mobility-based separation of structurally unstable from structurally stable RNA mutants, co-evolution analysis of high-throughput sequencing data, and energy-based structure prediction. 
## Installation
-------------
1.  Set up a python environment
    * We have tested it on Linux & MAC OS with Anacoda.
    * Install anacoda (https://www.anaconda.com)
    * Create your own environment (e.g., RNA) with python >= 3.7
2.  Requirement
    * Install the following software:
    * a). PEAR (https://github.com/tseemann/PEAR)
        PEAR  assembles  Illumina paired-end reads if the DNA fragment sizes are smaller than  twice  the length of reads.
    * b). Cutadapt (https://github.com/marcelm/cutadapt or https://cutadapt.readthedocs.io/en/stable/installation.html)
        Cutadapt finds and removes adapter sequences, primers, poly-A tails and other types of unwanted sequence from your high-throughput sequencing reads. Install: [pip install cutadapt]
    * c). FASTX-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
        The FASTX-Toolkit is a collection of command line tools for Short-Reads FASTA/FASTQ files preprocessing. Insall: [conda install fastx_toolkit]
    * d). BRiQ (https://github.com/Jian-Zhan/RNA-BRiQ) for 3D structure modeling.
3.  Clone RNA-MobiSeq pipeline (git@github.com:RNA-folding-lab/RNA-MobiSeq.git)

## Usage
-------------
0.  Perform high-throughput experiment (see Ref: Tang et al. RNA-MobiSeq: Deep mutational scanning and mobility-based selection for RNA base-pair and tertiary structure inference. bioRxiv 2024.02.21.581340; doi: https://doi.org/10.1101/2024.02.21.581340)：
	* a). Random Mutagenesis: Generate randomly mutated RNA sequences using error-prone PCR or Oligo pool technology.
	* b). Library Construction: Transfect mutated sequences into DH5α cells to construct DNA libraries.
	* c). In Vitro Transcription: Transcribe the DNA libraries into RNA sequences in vitro.
	* d). Native PAGE Separation: Separate the transcribed RNAs by native PAGE based on their mobility. 
	* e). Band Selection: Excise the bands with mobility matching the wild-type RNA sequences.
	* f). Reverse Transcription: Perform reverse transcription on the selected RNA sequences.
	* g). High-Throughput Sequencing: Sequence both pre- and post-selection RNA sequences using high-throughput sequencing.

1.  Put the high-throughput sequencing data of pre- and post-selection RNA (including forward_reads_file & reverse_reads_file) into "Metadata/" folder (e.g., 5TPY_pre_1.fastq, 5TPY_pre_2.fastq,5TPY_post_1.fastq,5TPY_post_1.fastq)
    * Put the RNA native sequence (e.g., 5TPY.fasta) into "native_data" folder

2.  Preprocess the sequencing reads data
    * Open the preprocess_sequencing_reads.sh file, and change the name of RNA (as well as pre/post_forward/reverse_reads_file and target_pre/post), the seq_adapter3/5, and the fastx_path based on your case.
    * Run preprocess_sequencing_reads.sh
    * '''sh
    * bash preprocess_sequencing_reads.sh
    * '''
    * If the script runs successfully, you will find the xx.var.ra files in the Processed/ directory, which can be used for further base pairing inference.
3.  Inferring base-pairing using CODA2 and predicting 2D & 3D structures using MC & RNA-BRiQ, respectively.
    * Open run.sh file, and change the name of RNA, OUT_PATH, BRiQ_PATH, as well as the several hyperparameters (e.g., C, gamma, sd_cut)
    * '''sh
    * bash run.sh
    * '''


## References
-------------
1.  Zhang J, Kobert K, Flouri T, Stamatakis A. PEAR: a fast and accurate Illumina Paired-End reAd mergeR. Bioinformatics. 2014 Mar 1;30(5):614-20. doi: 10.1093/bioinformatics/btt593. 
2.  Xiong P, Wu R, Zhan J, Zhou Y. Pairing a high-resolution statistical potential with a nucleobase-centric sampling algorithm for improving RNA model refinement. Nat Commun. 2021;12(1):2777. doi: 10.1038/s41467-021-23100-4
3.  Zhang Z, Xiong P, Zhang T, Wang J, Zhan J, Zhou Y. Accurate inference of the full base-pairing structure of RNA by deep mutational scanning and covariation-induced deviation of activity. Nucleic Acids Res. 2020;48(3):1451-1465. doi: 10.1093/nar/gkz1192. 
4.  Tang J, Shi Y, Zhang Z, Luo D, Zhan J, Zhou Y. RNA-MobiSeq: Deep mutational scanning and mobility-based selection for RNA structure inference. bioRxiv, 2024; doi: https://doi.org/10.1101/2024.02.21.581340
