#!/bin/csh -f
#PBS -q small
#PBS -N inmap_test
#PBS -l nodes=1:ppn=16
#PBS -o test.out
#PBS -e test.err
#PBS -l walltime=200:00:00

setenv NSLOTS `cat $PBS_NODEFILE|wc -l`
ulimit -s unlimited
cat $PBS_NODEFILE > host
echo $NSLOTS > nslots

setenv INMAP_ROOT_DIR /data/home/wuruili/WORK/11.InMAP/inmap161_China/

cd $INMAP_ROOT_DIR

source renv.sh

$INMAP_ROOT_DIR/inmap run steady --config=$INMAP_ROOT_DIR/China/configExample.toml
