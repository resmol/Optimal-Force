#!/bin/bash
#------------------------------------------------------------------------------#
#      Optimal force calculation                                               #
#------------------------------------------------------------------------------#

# Set up the number of steps (step variable):
MaxForce=$(grep 'FMAX' GENERAL.INPUT | tr " = " "\n" | tail -1)
IncForce=$(grep 'INCFORCE' GENERAL.INPUT | tr " = " "\n" | tail -1)
step=$(echo "${MaxForce}/${IncForce}" | bc)
Total=$step

export HomeDir=$PWD

#formchk  state0.chk
#formchk  state1.chk

#cp state0.chk state0_ini.chk
#cp state1.chk state1_ini.chk
cp state0.fchk state0_ini.fchk
cp state1.fchk state1_ini.fchk
cp state0.log state0_ini.log
cp state1.log state1_ini.log
mkdir Traj_Info
mkdir Traj_Info/General

cp state0.log OUTPUT0.molden
cp state1.log OUTPUT1.molden

#formchk state0.chk
#formchk state1.chk

var1=$(grep 'ANALYTICALSTATE1' GENERAL.INPUT | tr " = " "\n" | tail -1)
var2=$(grep 'ANALYTICALSTATE2' GENERAL.INPUT | tr " = " "\n" | tail -1)

#Loop
while ((step-->0))
do

let suma=$Total-$step
let suma2=$suma-1
i=$(grep -n 'PASO' GENERAL.INPUT | awk -v FS=":" '{print $1}') ; sed -i "${i}s/${suma2}/${suma}/g" GENERAL.INPUT
OptimalForce_CompleteSurface.py
wait

if [ $var1 == 'No' ] || [ $var1 == 'NO' ] || [ $var1 == 'no' ]; then
mv state0.log Traj_Info/General/state0_.log
#Launch_g16 state0.com
Launch_gaussian state0.com
cp state0.com Traj_Info/General/state0_.com
fi
if [ $var2 == 'No' ] || [ $var2 == 'NO' ] || [ $var2 == 'no' ]; then
mv state1.log Traj_Info/General/state1_.log
#Launch_g16 state1.com
Launch_gaussian state1.com
cp state1.com Traj_Info/General/state1_.com
fi

FileFound0=$(find ./ -name state0.log)
FileFound1=$(find ./ -name state1.log)
while [ -z $FileFound0 ] || [ -z $FileFound1 ] 
do
sleep 5s
FileFound0=$(find ./ -name state0.log)
FileFound1=$(find ./ -name state1.log)
done

if [ $var1 == 'No' ] || [ $var1 == 'NO' ] || [ $var1 == 'no' ]; then
mv state0.fchk Traj_Info/General/state0_.fchk
formchk state0.chk
cat state0.log >> OUTPUT0.molden
rm state0.sh* state0.Info state0.chk
fi
if [ $var2 == 'No' ] || [ $var2 == 'NO' ] || [ $var2 == 'no' ]; then
mv state1.fchk Traj_Info/General/state1_.fchk
formchk state1.chk
cat state1.log >> OUTPUT1.molden
rm state1.sh* state1.Info state1.chk
fi
  
cp -r Traj_Info/General Traj_Info/${suma}

done

let suma=$Total-$step
let suma2=$suma-1
i=$(grep -n 'PASO' GENERAL.INPUT | awk -v FS=":" '{print $1}') ; sed -i "${i}s/${suma2}/${suma}/g" GENERAL.INPUT
OptimalForce_CompleteSurface.py
wait

mv state0.log Traj_Info/General/state0_.log
mv state1.log Traj_Info/General/state1_.log
mv state0.fchk Traj_Info/General/state0_.fchk
mv state1.fchk Traj_Info/General/state1_.fchk

cp -r Traj_Info/General Traj_Info/${suma}

#cp OUTPUT0.molden $WorkDir/Prueba0.molden
#cp OUTPUT1.molden $WorkDir/Prueba1.molden

cat > utils.out << BRIHUEGA
grep "Gaussian Energy" OUTPUT | awk {'print $4'}
grep "Predicted Energy" OUTPUT | awk {'print $4'}
grep "Gaussian modForce" OUTPUT | awk {'print $4'}
grep "Predicted modForce" OUTPUT | awk {'print $4'}
for i in $(seq 1 11); do cat "${i}"/jmol_Fopt_disp Jmol >> Jmol ; done
ps -u alexander | grep sh
BRIHUEGA

