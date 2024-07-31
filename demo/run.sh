#!/bin/sh
RNA="5TPY"
data_file="Processed/${RNA}.var.ra"
NATIVE_PATH="native_data"
seq_file=${NATIVE_PATH}/${RNA}.fasta
native_2D=${NATIVE_PATH}/${RNA}-native.txt
OUT_PATH="./"
BRiQ_PATH="/home/asia/Asia/Software/RNA-BRiQ"
score="${OUT_PATH}score/final_score.txt"
C=0.01
gamma=0.1
sd_cut=3.7
random_seed=0

#source activate base  
echo "#######################"
echo "Inferring base-pairing"
echo "#######################"
python CODA2.py $data_file $seq_file -C $C -g $gamma -s $sd_cut -n $native_2D #-o OUT_PATH

echo "#######################"
echo "Predicting 2D structure"
echo "#######################"
python Predict_2D.py $score $seq_file -n $native_2D -seed $random_seed

echo "#######################"
echo "Predicting 3D structure"
echo "#######################"

sec=$(head -n 1 ${OUT_PATH}2D/2D.dot)
seq=$(sed -n '2p' $seq_file)
length=${#seq}
nwc=$(printf "%.0s." $(seq 1 $length))

mkdir 3D/
cd 3D/
cat > input <<EOF
pdb init.pdb
seq $seq
sec $sec
nwc $nwc
EOF

cat >run_briq.sh <<EOF
#!/bin/sh
export BRiQ_BINPATH=${BRiQ_PATH}/build/bin 
export BRiQ_DATAPATH=${BRiQ_PATH}/BRiQ_data 

INPUT=input     # Input file 
OUTPDB=pred.pdb   # Output: refined structure
RANDOMSEED=$random_seed   # Random seed  

## Generate an initial PDB structure from the given sequence  
\$BRiQ_BINPATH/BRiQ_InitPDB $seq init.pdb  
## Predict RNA structure        
\$BRiQ_BINPATH/BRiQ_Predict \$INPUT \$OUTPDB \$RANDOMSEED
EOF

bash run_briq.sh

cd ../

echo "##################################################"
echo "Great! pred.pdb in 3D/ is the predicted structure"
echo "##################################################"