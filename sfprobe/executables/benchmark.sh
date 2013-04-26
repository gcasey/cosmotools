#!/bin/sh

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo $SCRIPT_DIR

datadir="$SCRIPT_DIR/../data/"

#cuda scaling over number of points
pointSize="
  31
  47
  63
  79
  95
  111
  127
"
maxcores="32"
logfile=`date '+timings-%Y-%m-%d.csv'`

function doCudaPointScaling {
  for size in $pointSize
  do
    echo "cudaProbe, $size" | tee -a $logfile
    for i in {1..20}
    do
    "$SCRIPT_DIR"/bin/cudaProbe $datadir 6.2831 31 $size 2 | sed -n 's/^Dax Probe Time://p' | tee -a $logfile
    done
  done
}

function doTbbPointScaling {
  for size in $pointSize
  do
    echo "tbbProbe, $size" | tee -a $logfile
    for i in {1..20}
    do
    "$SCRIPT_DIR"/bin/tbbProbe $datadir 6.2831 31 $size 2 $maxcores | sed -n 's/^Dax Probe Time://p' | tee -a $logfile
    done
  done
}

function doTbbSerialScaling {
  for size in $pointSize
  do
    echo "serialProbe, $size" | tee -a $logfile
    for i in {1..20}
    do
    "$SCRIPT_DIR"/bin/probe $datadir 6.2831 31 $size 2 | sed -n 's/^Serial Probe Time://p' | tee -a $logfile
    done
  done
}

function doTbbCoreScaling {
  for cores in {i..$maxcores}
  do
    echo "tbb scaling Probe, $size" | tee -a $logfile
    for i in {1..20}
    do
    "$SCRIPT_DIR"/bin/tbbProbe $datadir 6.2831 31 63 2 $cores | sed -n 's/^Dax Probe Time://p' | tee -a $logfile
    done
  done
}


doCudaPointScaling
doTbbPointScaling
doTbbSerialScaling
doTbbCoreScaling
