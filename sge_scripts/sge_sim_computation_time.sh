#!/usr/bin/env bash
#$ -V
#$ -cwd
#$ -N comp_time
#$ -m besa
#$ -q short.q
#$ -t 1-3
#$ -pe thread 64
#$ -M louis.lacoste+migale@agroparistech.fr
#$ -o logs/$JOB_NAME
#$ -e logs/$JOB_NAME

# Creating log directory if it doesn't exists
BASE_DIR="/home/$USER/work/code-colbisbm"
LOG_DIR=$(echo "$BASE_DIR/logs")

if [ ! -d "$LOG_DIR" ]; then
    mkdir -p $LOG_DIR
fi

# Finding directory
APPLICATIONS_DIR=$(echo "$BASE_DIR")

echo $APPLICATIONS_DIR

ARGID=$(($SGE_TASK_ID % 3))

case $ARGID in
    0)
    echo -n "Model will iid"
    MODE="iid"
    TEST_NULL=FALSE
    ;;
    1)
    echo -n "Model will be pirho with no null block"
    MODE="pirho"
    TEST_NULL=FALSE
    ;;
    2)
    echo -n "Model will be pirho with null blocks!"
    MODE="pirho"
    TEST_NULL=TRUE
    ;;
esac

Rscript "${APPLICATIONS_DIR}/simulations_computation_time.R.R" $MODE $TEST_NULL &>> logs/$JOB_NAME.$SGE_TASK_ID