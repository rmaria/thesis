#! /bin/sh -f

export tau="0.41"            #0.41     1
export DeltaGammaTau="1.44e-2"      #1.44e-2
export DeltaM="1.18e-2"
#export DeltaM="6*TMath::Pi()/2"      #1.18e-2   2*TMath::Pi()/2
export ModLambdaf="1.0"
export PhaseLambdaf="5.0"
export ResModelBias="0.0"
export ResModelSigma="0.00001"
export Omega="0.0"
export DeltaOmega="0.0"
export Nevents="12000000"
#export Seed=834628364 
export Seed=$1   
export Option="SingleSided"
export OutFile="Plots/Example"

source /libcern/root/5.34.18/sl6.3-x86_64/setup.sh

chmod +x AppExample

./AppExample  $tau  $DeltaGammaTau  $DeltaM  $ModLambdaf  $PhaseLambdaf  $ResModelBias $ResModelSigma	$Omega  $DeltaOmega  $Nevents  $Seed  $Option   $OutFile  

