#!/bin/bash

khome=/global/scratch/nathan_sandford/kurucz/kurucz_home
krun=/global/scratch/nathan_sandford/kurucz/kurucz_run

# Asplund et al. (2009) abundances of first 99 elements
elems=(0.92068 0.07837 -10.99 -10.66 -9.34 -3.61 -4.21 -3.35 -7.48 -4.11 -5.80 -4.44 -5.59 -4.53 -6.63 -4.92 -6.54 -5.64 -7.01 -5.70 -8.89 -7.09 -8.11 -6.40 -6.61 -4.54 -7.05 -5.82 -7.85 -7.48 -9.00 -8.39 -9.74 -8.70 -9.50 -8.79 -9.52 -9.17 -9.83 -9.46 -10.58 -10.16 -20.00 -10.29 -11.13 -10.47 -11.10 -10.33 -11.24 -10.00 -11.03 -9.86 -10.49 -9.80 -10.96 -9.86 -10.94 -10.46 -11.32 -10.62 -20.00 -11.08 -11.52 -10.97 -11.74 -10.94 -11.56 -11.12 -11.94 -11.20 -11.94 -11.19 -12.16 -11.19 -11.78 -10.64 -10.66 -10.42 -11.12 -10.87 -11.14 -10.29 -11.39 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -12.02 -20.00 -12.58 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00 -20.00)

#output directory
head=$1
outdir=$krun/grids/$head/

#set the input directory
ihead=$2
indir=$krun/grids/$ihead/atm/

#get the model filename
arr=(`ls $indir`)
model=${arr[0]}

#we want a filename that all the header info stripped off
tmodel=${model#$ihead}
tmodel=${tmodel#_}
teff=${tmodel:2:4}
logg=${tmodel:7:4}
outfile=${head}_t0${teff}g${logg}.atm
infile=${ihead}_t0${teff}g${logg}.atm

### for temperature > 10000K ###
#teff=${tmodel:1:5}
#logg=${tmodel:7:4}
#outfile=${head}_t${teff}g${logg}.atm
#infile=${ihead}_t${teff}g${logg}.atm

#echo " "
#date "+%Y-%m-%d %H:%M:%S"

#------------------------------------------------------------------------------
# convergence indicator
conv_ind=0

# at most run 150 iterations
iter_count=0

# loop until converged
while [ $conv_ind -eq 0 ]; do

    #create the directory
    #echo "creating new directory:" $outdir
    mkdir -p $outdir
    mkdir -p $outdir/atm $outdir/molnden $outdir/outfiles $outdir/spec 
    mkdir -p $outdir/flux $outdir/fail $outdir/crash $outdir/sed

    mkdir $krun/atlas12/tmp_$head
    #echo "Moving into temporary working directory...."
    cd $krun/atlas12/tmp_$head

    rm -f fort.*
    ln -s $khome/lines/molecules.new fort.2
    ln -s $indir$model fort.3
    ln -s $khome/lines/gfpred29dec2014.bin fort.11
    ln -s $khome/lines/lowobsat12.bin fort.111
    ln -s $khome/lines/hilines.bin fort.21
    ln -s $khome/lines/diatomicspacksrt.bin fort.31

    #include all molecules
    echo "including TiO and H2O lines for model" $head
    ln -s $khome/molecules/tio/schwenke.bin fort.41 
    ln -s $khome/molecules/h2o/h2ofastfix.bin fort.51

    #run a model with the original input atm abundances
    echo "Running atlas12.exe on model" $head
    $khome/bin/atlas12.exe<<EOF>at12a.out
MOLECULES ON
READ MOLECULES
READ PUNCH
READ LINES
CONVECTION OVER 1.25 0 36
ITERATIONS 1 PRINT 1 PUNCH 0
BEGIN
END
EOF

    mv fort.12 tmp.bin
    rm -f fort.*
    ln -s $indir$model fort.3
    ln -s $khome/lines/nltelinobsat12.bin fort.19
    ln -s $khome/lines/molecules.new fort.2
    ln -s tmp.bin fort.12

    outf="${outfile/atm/out}"

    $khome/bin/atlas12.exe<<EOF> ${outf}
MOLECULES ON
READ MOLECULES
READ PUNCH
TITLE ATLAS12 l/H=1.25
OPACITY ON LINES
OPACITY ON XLINES
CONVECTION OVER 1.25 0 36
ITERATIONS 30
PRINT 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
PUNCH 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 2
SCALE MODEL 80 -6.875 0.125 ${teff}. ${logg}
VTURB 1.00E+05
BEGIN
END
EOF

    echo "atlas is finished running on model" $head

    iter="${outfile/atm/iter}"
    tcorr="${outfile/atm/tcorr}"

    #test if the atm file was sucessfully created
    if [ -f fort.7 ]; then

        #test convergence
	echo "testing model" $head "for convergence"
	$khome/bin/checkconv.exe
	mv fort.66 $outdir/outfiles/$iter
	mv fort.67 $outdir/outfiles/$tcorr

        #save the atm file
        mv fort.7 $outdir/atm/$outfile

        # if the model has not converged, continue to run
	if [ -f model.failed -a $iter_count -lt 4 ]; then
	    echo "model" $head "has not converged, continue to run..."

            # remove temporary file
	    cd ../
	    rm -rf tmp_$head
	    rm -rf $indir/$infile
	    
	    # write the current atmosphere
	    cp $outdir/atm/$outfile $indir/$infile
	    rm -rf $outdir

	    # at most run 150 iterations
	    iter_count=$[$iter_count+1]

	# if already converged, then stop
	else

            #save the flux output
	    flux="${outfile/atm/flux}"
	    mv fort.8 $outdir/flux/$flux

            #clean up
	    rm -f fort.* tmp.bin

	    # run synthe
	    echo "running synthe on model" $head
	    $krun/synthe/synthe.sh $head
	    spec="${outfile/atm/spec}"
	    echo "model" $head "completed; compressing spectra"
	    gzip $outdir/spec/$spec

	    # break the code
	    conv_ind=1
	fi
	
    # if no atm was made (code break)
    else
	echo "atlas12 apparently crashed on model" $head ", exiting..."

	# break the code
	conv_ind=1
	
	cp fort.3 $outdir/crash/$model

	if [ -f fort.66 ]; then
	    mv fort.66 $outdir/outfiles/$iter
	fi
    fi
done

#clean up
cd ../
rm -rf tmp_$head

#date "+%Y-%m-%d %H:%M:%S"
