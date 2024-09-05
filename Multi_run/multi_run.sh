#!/bin/sh
RNA="5TPY"
root_path="~/.../RNA-Mobiseq"       #depending on your case
data_file="${root_path}/Multi_run/Processed/${RNA}.var.ra"  #fitness data dir, can be changed by yourself
NATIVE_PATH="${root_path}/Multi_run/native_data"  #native seq data dir
seq_file=${NATIVE_PATH}/${RNA}.fasta
native_2D=${NATIVE_PATH}/${RNA}-native.txt        #(option, only for prediction evalutation)
OUT_PATH="$(pwd)"
CODA2="${root_path}/CODA2.py"
Sort_Scoring="${root_path}/Multi_run/sort_scoring_by_U.py"
Predict_2D="${root_path}/Predict_2D.py"
BRiQ_PATH="~/.../RNA-BRiQ"       ##RNA-BRiQ path on your computer
score="${OUT_PATH}/score/final_score.txt"
N_MC=1000
N_BRiQ=1000
Max_jobs=10
C=(0.01 0.1 1.0 10.0 100.0)
gamma=(0.001,0.01,0.1,1.0,10.0)
sd_cut=(1.5,2.0,2.5,3.0,3.5,4.0)

#source activate base
###limits the number of concurrent requests
run_with_limit() {
    while [ "$(jobs | wc -l)" -ge "$Max_jobs" ]; do
        sleep 1
    done
    "$@" &
} 

echo "#######################"
echo "Inferring base-pairing"
echo "#######################"
mkdir -p scoring score
cd scoring/

run_scoring() {
  local i=$1
  local j=$2
  local k=$3
  local dir="${i}/${j}/${k}"
  
  mkdir -p "$dir"
  cp "$CODA2" "$dir/"
  #ln -s "$CODA2" "$dir/"
  
  echo "Running: $i-$j-$k"
  (
  	cd "$dir" && 
  	python CODA2.py "$data_file" "$seq_file" -C "$i" -g "$j" -s "$k" > out.txt #-n $native_2D -o OUT_PATH )
  )
  rm "$dir"/*.py
}
for i in "${C[@]}"; do
    for j in "${gamma[@]}"; do
        for k in "${sd_cut[@]}"; do
            run_with_limit run_scoring "$i" "$j" "$k"
        done
    done
done
wait  ##wait all jobs end

cp "$Sort_Scoring" .
for i in "${C[@]}"; do
  for j in ${gamma[@]}; do 
    for k in ${sd_cut[@]}; do 
       path="$i/$j/$k/"
       python sort_scoring_by_U.py "$path" "$seq_file" -n "$native_2D"
      done
   done
done
rm *.py
best_score_path=$(awk 'NR == 1 || $2 < min {min = $2; path = $1} END {print path}' U.txt)
cd ../
cp "scoring/$best_score_path/score/final_score.txt" "score"

echo "#######################"
echo "Predicting 2D structure"
echo "#######################"

for i in $(seq 1 ${N_MC}); do
   run_with_limit python ../Predict_2D.py "$score" "$seq_file" -seed "$i" > out_2D.txt  #-n "$native_2D"
done
wait

###Get top 5 2D structures
python analyze_multi_2D.py ./ "${N_MC}" #-n "$native_2D"
mkdir -p 2D_top5
typical_list=$(grep "Typical:" "cluster_2D.txt" | awk -F'[][]' '{print $2}' | tr -d ' ' | cut -d',' -f1-5)
IFS=',' read -r -a numbers <<< "$typical_list"
for i in "${!numbers[@]}"; do
  num="${numbers[$i]}"
  src_file="2D/${num}/2D.dot"
  dest_file="2D_top5/top$((i+1)).dot"
  if [ -f "$src_file" ]; then
    cp "$src_file" "$dest_file"
  else
    echo "File $src_file not found!"
  fi
done

echo "#######################"
echo "Predicting 3D structure"
echo "#######################"
# Function to run the BRiQ prediction
mkdir -p 3D
cd 3D/
run_briq_prediction() {
    local i=$1
    local seq=$(sed -n '2p' "$seq_file")
    local sec=$(head -n 1 "${OUT_PATH}/2D_top5/top${i}.dot")
    local length=${#seq}
    local nwc=$(printf "%.0s." $(seq 1 $length))
    
    for j in $(seq 1 ${N_BRiQ}); do
        mkdir -p "top${i}/${j}"
        cd "top${i}/${j}"

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
OUTPDB=${j}.pdb   # Output: refined structure
RANDOMSEED=${j}   # Random seed  

## Generate an initial PDB structure from the given sequence  
\$BRiQ_BINPATH/BRiQ_InitPDB $seq init.pdb  
## Predict RNA structure        
\$BRiQ_BINPATH/BRiQ_Predict \$INPUT \$OUTPDB \$RANDOMSEED
rm init.pdb
EOF

        bash run_briq.sh
        cd ../../
    done
}
# Run the function in parallel for each top file (top1.dot to top5.dot)
for i in {1..5}; do
    run_with_limit run_briq_prediction "$i"
done
wait

cd 3D/
##Get Top 5 3D structures
mkdir -p 3D_top5
for TOP in {1..5}; do
   	cd top${TOP}
   	: > Energy_${TOP}.txt
   	for i in $(seq 1 ${N_BRiQ}); do
       	cd $i
       	pdb="${i}.pdb"
       	energy=$(awk 'END{print $2}' "$pdb")
       	echo "${TOP} ${i} ${energy}" >> ../Energy_${TOP}.txt
       	cd ../
   	done   	
   	cat Energy_${TOP}.txt >> ../Energy.txt
   	cd ../
done

sort -k 3 -n Energy.txt -o Energy_sorted.txt
k=1
head -n 5 Energy_sorted.txt | while read -r line; do
    TOP=$(echo $line | awk '{print $1}')
    i=$(echo $line | awk '{print $2}')
    cp top${TOP}/${i}/${i}.pdb 3D_top5/top${k}.pdb
    k=$((k + 1))
done

echo "########################################################"
echo "Great! 3D structures are predicted: Top 5 in 3D/3D_top5"
echo "########################################################"