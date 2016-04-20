
source setenv.sh

projectName="Example"

output="LQ_pars meanValues"

model="MKM"
calculusType="rapidMKM"

precision=0.015
parallelismType=0

doses="1 2 3 4 5 6"

cellType="HSG"
MKM_alpha0=0.312
MKM_beta0=0.073
MKM_rNucleus=4.611
MKM_rDomain=0.365

ion="H"

trackMode="histogram"

lets="0.25 0.5 0.75 1 1.3 1.6 1.9 2.2 2.5 3 3.5 4 4.5 5 5.5 6 7 8 9 10.5 12 14 17 20 25 30 35 40 45 50 55 60 65 70 75 80 84.3"

./survival	-projectName $projectName \
		-output $output \
		-model $model \
		-calculusType $calculusType \
		-precision $precision \
		-parallelismType $parallelismType \
		-doses $doses \
		-cellType $cellType \
		-MKM_alpha0 $MKM_alpha0 \
		-MKM_beta0 $MKM_beta0 \
		-MKM_rNucleus $MKM_rNucleus \
		-MKM_rDomain $MKM_rDomain \
		-ion $ion \
		-trackMode $trackMode \
		-lets $lets \


