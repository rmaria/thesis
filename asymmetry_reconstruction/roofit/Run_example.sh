#! /bin/sh -f

export tau="0.410"                # to test, put 2
export DeltaGammaTau="0"      #1.44e-2
export DeltaM="1.18"            #1.18e-2
#export DeltaM="4*3.144/2"
export ModLambdaf="1.0"
export PhaseLambdaf="10.0"
export ResModelBias="0.0"
export ResModelSigma="0."
export Omega="0.0"
export DeltaOmega="0.0"
export Nevents="500000"
export Seed=77774365
export Option="SingleSided"
export OutFile="Plots/Example"

./AppExample  $tau  $DeltaGammaTau  $DeltaM  $ModLambdaf  $PhaseLambdaf  $ResModelBias  $ResModelSigma  $Omega  $DeltaOmega  $Nevents    $Seed   $Option   $OutFile

