#!/usr/bin/env bash
#$ -V
#$ -cwd
#$ -N comp_time
#$ -m besa
#$ -q long.q
#$ -t 1-8
#$ -pe thread 64
#$ -M louis.lacoste+migale@agroparistech.fr
#$ -o logs/$JOB_NAME
#$ -e logs/$JOB_NAME

BASE_DIR="/home/$USER/work/code-colbisbm"
LOG_DIR="$BASE_DIR/logs"

mkdir -p "$LOG_DIR"

APPLICATIONS_DIR="$BASE_DIR"

echo "$APPLICATIONS_DIR"

# Convert task id (1-8) to index (0-7)
IDX=$((SGE_TASK_ID - 1))

# Factor levels
MODELS=("iid" "pirho")
MS=(3 10)
NS=(100 300)

# Decode indices
MODEL=${MODELS[$((IDX / 4))]}
M=${MS[$(((IDX % 4) / 2))]}
N=${NS[$((IDX % 2))]}

echo "Model = $MODEL"
echo "M = $M"
echo "n = $N"

Rscript "${APPLICATIONS_DIR}/simulations_computation_time.R" \
    --model "$MODEL" \
    --nb-networks "$M" \
    --nb-nodes "$N" \
    &>> "logs/$JOB_NAME.$SGE_TASK_ID"