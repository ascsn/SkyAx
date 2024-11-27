echo "&main outnam='s$1',title='170Yb cut40 shells$1',iprint=1000,itmax=20000" > for005.s$1
echo "/" >> for005.s$1
echo "&force" >> for005.s$1
echo " name='none',pairstab=-0D0,HFBwidth=1.0d0," >> for005.s$1
echo " t0=-1879.640018,t1=313.7493427,t2=112.6762700,t3=12527.38921," >> for005.s$1
echo " x0=0.2585452462,x1=-0.3816889952,x2=-2.823640993,x3=0.1232283530," >> for005.s$1
echo " t4=124.633,b4p=34.1116989,power=0.3," >> for005.s$1
echo " h2m=20.749811,20.721249," >> for005.s$1
echo " chg=1D0,ex=1,riftm=0D0,zpe=1," >> for005.s$1
echo " ipair=8,v0prot=1563,v0neut=1448,rho0pr=0.2011591," >> for005.s$1
echo " pair_cutoff=-40,tregpair=F" >> for005.s$1
echo "/" >> for005.s$1
echo "&grids tsym=F,nz=51,nr=25,dr=0.7D0,dz=0.7D0,idrtyp=0," >> for005.s$1
echo "/" >> for005.s$1
echo "&static nprot=70,nneut=100,e0inv=110D0,x0inv=0.3D0,addrho=0.15D0,addocc=0.3D0," >> for005.s$1
echo "deocc(1)=$1,deocc(2)=$1,errenr=1D-14,errvar=0.05,colerr=1D-9,itpini=0," >> for005.s$1
echo "tvaryx0=T,inodiag=9999/" >> for005.s$1
echo "&constraint" >> for005.s$1
echo "/" >> for005.s$1
echo "&thermal /" >> for005.s$1


~/skyax/SkyAx/skyax.gfortran.seq < for005.s$1 > /dev/null
