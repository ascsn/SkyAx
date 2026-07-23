#!/bin/bash
# Batch file for systematic variation of cutoff and folding width.
# Test case is 120Sn with 51 radial grid points and 17 shell above
# the Fermi energy. Proton pairing is practically switched off.
# The parameters are: 
# i = cutoff 
# hfbwidth = folding widths
# pairStrength = v0neut


nucleus="120Sn"
edf="SVbas"

# Set HFB widths
hfbwidths=("0.0d0" "1.0d0" "1.5d0" "2.0d0" "3.0d0")

# Set BCS pairing strengths corresponding to those hfbwidths
# fixing lower cutoff pairing gap
#pairStrengthsBCS=("700" "1884" "18700" "43510")

# Set HFB pairing strengths corresponding to those hfbwidths
#pairStrengthsHFB=("530" "1610" "13460" "35064")
pairStrengthsHFB=("342" "1610" "10815" "30301")

ipair=(7)    # for BCS
#ipair=(7 8)    # for HFB

# Set hfbwidth
hfbwidth="1.5d0"

for j in "${ipair[@]}"; do
if [ $j -eq 7 ]; then
    cond="bcs"
    pairStrength="890"
    #combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsBCS[@]}"))
else
    cond="hfb"
#    pairStrength="30290."
    combined=$(paste -d "," <(printf "%s\n" "${hfbwidths[@]}") <(printf "%s\n" "${pairStrengthsHFB[@]}"))
fi
#while IFS=, read -r hfbwidth pairStrength; do
echo "cond = ${cond}, hfbwidth=${hfbwidth}, pairStrength = ${pairStrength}" 
for i in $(seq 99 1 100); do
echo "&main outnam='tst',title='test 120Sn',iprint=1000,itmax=1000" > for005
echo "/" >> for005
echo "&force" >> for005
echo " name='none',pairstab=-0D0,HFBwidth=${hfbwidth}," >> for005
echo " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921," >> for005
echo " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530," >> for005
echo " t4=124.633,b4p=34.1116989,power=0.3," >> for005
echo " h2m=20.749811,20.721249," >> for005
echo " chg=1D0,ex=1,riftm=0D0,zpe=1," >> for005
echo " ipair=7,v0prot=674.618,v0neut=${pairStrength},rho0pr=0.2011591," >> for005
echo " pair_cutoff=-${i},tregpair=F" >> for005
echo "/" >> for005
echo "&grids tsym=F,nz=81,nr=41,dr=0.5D0,dz=0.5D0,idrtyp=0," >> for005
echo "/" >> for005
echo "&static nprot=50,nneut=70,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0," >> for005
echo "deocc(1)=3.1D0,deocc(2)=3.1D0,errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0," >> for005
echo "tvaryx0=T,inodiag=9999/" >> for005
echo "&constraint" >> for005 
echo "/" >> for005
echo "&thermal /" >> for005


../skyax.ifort.seq < for005 > /dev/null

fileName="collect_cat_${edf}_${nucleus}_${hfbwidth}_${cond}.dat"
echo "cutoff=-${i},hfbwidth=${hfbwidth},pair.strength=${pairStrength}" >> ${fileName}
var1="grep '# Nmass' ensurf.tst | awk '{print \"Nmass:\", \$3, \"Ncharge:\", \$5}'"
var2="grep -A1 'energy' ensurf.tst | tail -n1 | awk '{print \"energy:\", \$13, \"epairp:\", \$14, \"epairn:\", \$15, \"rmsc:\", \$19}'"
var3="grep -A2 'gap' out_details.tst | awk '/proton:/ {print \"proton gap:\", \$4} /neutron:/ {print \"neutron gap:\", \$4}'"

paste -d ' ' <(eval "$var1") <(eval "$var2") <(eval "$var3") >> ${fileName}

done    
echo "${cond} Done!"
mv ${fileName} collect_cat_120Sn/
#done <<< "${combined}"
done