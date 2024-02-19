for i in 80 82 84 86 88 90 92 94 96 98 100 102

do

#  start from prolate
echo "&main outnam='GdC',title='Gd chain',iprint=1000,itmax=1000/" > for005
echo "&force  name='SV-bas',ipair=7/" >> for005
echo "&grids tsym=T,nr=23,nz=23/" >> for005
echo "&static nprot=64,nneut=$i,e0inv=110D0,x0inv=0.4D0,addrho=0.15D0,addocc=1D0," >> for005
echo "        errenr=1D-9,errvar=0.05D0,betain=0.3D0,tvaryx0=T/" >> for005
echo "&constraint swicon2=0D0/" >> for005
echo "&thermal /" >> for005

../../skyax.ifort.seq < for005 > /dev/null

echo "$i" >> collect.GD_SV-bas
tail -n 1 ensurf.GdC >> collect.GD_SV-bas
mv ensurf.GdC ensurf.Gd_SV-bas_N$i
mv radii.GdC radii.Gd_SV-bas_N$i

#start from oblate
echo "&main outnam='GdC',title='Gd chain',iprint=1000,itmax=1000/" > for005
echo "&force  name='SV-bas',ipair=7/" >> for005
echo "&grids tsym=T,nr=23,nz=23/" >> for005
echo "&static nprot=64,nneut=$i,e0inv=110D0,x0inv=0.4D0,addrho=0.15D0,addocc=1D0," >> for005
echo "        errenr=1D-9,errvar=0.05D0,betain=-0.2D0,tvaryx0=T/" >> for005
echo "&constraint swicon2=0D0/" >> for005
echo "&thermal /" >> for005

../../skyax.ifort.seq < for005 > /dev/null

tail -n 1 ensurf.GdC >> collect.GD_SV-bas
cat ensurf.GdC >> ensurf.Gd_SV-bas_N$i
cat radii.GdC >> radii.Gd_SV-bas_N$i

rm ensurf.GdC
rm radii.GdC


done



