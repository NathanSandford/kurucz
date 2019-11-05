#!/bin/bash

#date

khome=/global/scratch/nathan_sandford/kurucz/kurucz_home
krun=/global/scratch/nathan_sandford/kurucz/kurucz_run
bindir=/global/scratch/nathan_sandford/kurucz/kurucz_home/bin

indir="${krun}/grids/$1/atm/"
outdir="${krun}/grids/$1/spec/"
moldir="${krun}/grids/$1/molnden/"

#echo "synthe output Dir: ${outdir}"

arr=(`ls $indir`)

#echo "Moving into temporary working directory...."
mkdir ${krun}/synthe/stmp_$1
cd ${krun}/synthe/stmp_$1

model=${arr[0]}
linedir="Lines_v5_PL"

asc="${model/atm/spec}"
asc="${asc/dat/spec}"

echo "Using ${linedir} for model" $1

#generate synthe-ready input file
echo "Passing model" $1 "to at12tosyn.exe"
$bindir/at12tosyn.exe $indir$model $model

#re-link files each time
cp ${khome}/lines/molecules.dat fort.2
cp ${khome}/lines/continua.dat  fort.17

$bindir/xnfpelsyn.exe < $model  > /dev/null

#link the input files generated from synthe.setup
cp ${khome}/synthe/${linedir}/tfort.12 fort.12
cp ${khome}/synthe/${linedir}/tfort.14 fort.14
cp ${khome}/synthe/${linedir}/tfort.19 fort.19
cp ${khome}/synthe/${linedir}/tfort.20 fort.20
cp ${khome}/synthe/${linedir}/tfort.93 fort.93

cp ${khome}/lines/he1tables.dat fort.18

# change C12/C13 ratio
#$bindir/globaladj.exe<<EOF > /dev/null
#12
#13
#74
#EOF

#mv fort.221 fort.19
#mv fort.222 fort.20
#mv fort.223 fort.12
#mv fort.224 fort.14

#run synthe, the main program
echo "running synthe.exe on model" $1
$bindir/synthe.exe > /dev/null

ln -s $model fort.5
ln -s $khome/infiles/spectrv_std.input fort.25
    
echo "running spectrv.exe on model" $1
$bindir/spectrv.exe > /dev/null

#save the molecular number density profiles
rm -f fort.2 fort.5 fort.25 fort.33

#convert the spectrum into an ascii file
mv fort.7 fort.1

echo "running syntoascanga.exe on model" $1
$bindir/syntoascanga.exe  > /dev/null
mv specfile.dat $outdir$asc

#clean up
rm -f fort.*

#echo "removing working directory"
cd ../
rm -r stmp_$1

#date

