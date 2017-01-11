#! /bin/sh -f


export Run_nb="$1"
export minus_factor="1405000000"
export rand_sec=`date +%s`
echo "rand_sec = $rand_sec"
export rand_sec_trunc="$(expr $rand_sec - $minus_factor)"
echo "rand_sec_trunc=$rand_sec_trunc"
export Seed=$rand_sec_trunc$Run_nb
echo "Seed = $Seed" 
echo "Run number = $Run_nb"
export tau="0.41"					#0.41
export DeltaGammaTau="1.44e-2"      #1.44e-2
#export DeltaM="1.18e-2"
export DeltaM="0.0118" 			#0.0118
export ModLambdaf="1.0"
export PhaseLambdaf="5.0"
export ResModelBias="0.0"
#export ResModelSigma="0.2"  #0.0001
export Omega="0.0"
export DeltaOmega="0.0"
export Nevents="5000"
#export Seed="486238            #891879314043"
export Option="SingleSided"
export OutFile="Plots/Example"

./AppExample $Seed $tau  $DeltaGammaTau  $DeltaM  $ModLambdaf  $PhaseLambdaf  $ResModelBias  $Omega  $DeltaOmega  $Nevents  $Option   $OutFile

