!		SPARSE NETWORK MODEL 
!		THREE POPULATIONS: se, si, sp
!		CONECTIVITY: K COUPLINGS:
!		FIRST ORDER IMPLICIT INTEGRATION
!		INTERACTION: 1 EXPONENTIAL
!   FULLY OPTIMIZED (CONNECTIVITY, INTEGRATION) OF DEPEIR.F
!
!    E->I' Facilitating   I'->I depressing  E-E depressing
PROGRAM INTFIRE
  implicit none
  integer n,nc,kk,nu,nce,nci
  integer ne,ni,nn,nnn
  real*8 kkee,kkei,kkii,kkie
  real*8 kkee0,kkei0,kkii0,kkie0
  parameter(n=64000)  !poner multiplo de nu
  parameter(nu=1) 
  parameter (kk=300) 
  parameter(nc=2100,nce=2100,nci=800)
  parameter(ne=64000,ni=16000)
  integer i,j,k,iu,it,itmax,iseed1,iseed2
  integer it1,it2,it3,it4,it5,itrans
  integer i1,i2,i3
  real*8, allocatable :: switche(:),switchi(:)
  real*8, allocatable :: se(:),si(:)		   !estados de las poblaciones
  real*8, allocatable :: ue(:),ui(:)		   !estados de las poblaciones
  real*8, allocatable :: semact(:),simact(:)
  real*8, allocatable :: seed(:),seem(:)
  real*8, allocatable :: averhe(:),averhi(:)
  real*8, allocatable :: averhe2(:),averhi2(:)
  real*8, allocatable :: averexcie(:),averexcie2(:)
  real*8, allocatable :: averinhibe(:),averinhibe2(:)
  real*8, allocatable :: averexcii(:),averexcii2(:)
  real*8, allocatable :: averinhibi(:),averinhibi2(:)
  integer, allocatable :: mee(:,:),mei(:,:),mii(:,:)   !matrices de conexion
  integer, allocatable :: mie(:,:)
  integer, allocatable :: ncee(:),ncei(:),ncii(:),ncie(:)
  real*8, allocatable :: he(:),hi(:)    !campos internos
  real*8, allocatable :: ree(:),rie(:),rei(:),rii(:)
  real*8, allocatable :: hee(:),hie(:),hei(:),hii(:)
  real*8, allocatable :: heeampa(:),hieampa(:)
  integer, allocatable :: nspike(:),nspiki(:)
  real*8, allocatable :: fee2(:),fei2(:),fii2(:)
  real*8, allocatable :: fie2(:)
  real*8, allocatable :: feeampa2(:),fieampa2(:)
  real*8 , allocatable ::ke1(:),ki1(:),kp1(:)
  real*8, allocatable :: ke2(:),ki2(:),kp2(:)
  real*8, allocatable :: sebis(:),sibis(:),spbis(:)
  real*8, allocatable :: tspike(:),tspikeold(:)
  real*8, allocatable :: tspiki(:),tspikiold(:)
  real*8, allocatable :: tspikp(:),tspikpold(:)
  real*8, allocatable :: sumspe(:),sumsp2e(:),sumspi(:)
  real*8, allocatable :: sumsp2i(:)
  real*8, allocatable :: sumecv2(:),sumicv2(:),sumpcv2(:)
  real*8 , allocatable :: dte(:),dteold(:),dti(:),dtiold(:),dtp(:),dtpold(:)
  real*8, allocatable :: fire(:),firi(:)
  real*8, allocatable :: cve(:),cvi(:),cvp(:)
  real*8, allocatable :: cv2e(:),cv2i(:),cv2p(:)
  real*8, allocatable :: cose(:),cosi(:),sine(:),sini(:)
  real*8, allocatable :: cose2(:),cosi2(:),sine2(:),sini2(:)
  real*8, allocatable :: thetae(:),thetai(:)
  real*8, allocatable :: c1e(:),c1i(:),s1e(:),s1i(:),m1e(:),m1i(:)
  real*8 , allocatable :: c2e(:),c2i(:),s2e(:),s2i(:),m2e(:),m2i(:)

  real*8 eps,umeps,taufac,taurec,bigu  !valor medio temporal
  real*8 taufacee,taurecee,biguee
  real*8 gee,geeampa,gei,gii,gie,ie,ii    !pesos de conexion
  real*8 ieon,ieoff,iion,iioff
  real*8 g1ee,g1ei,g1ii,g1ie    !pesos de conexion
  real*8 sigmaee,sigmaei,sigmaii,sigmaie    !pesos de conexion
  real*8 gpec
  real*8 tauact,epsact,umact
  real*8 gee0,gampa0,gei0,gii0,gie0,gieampa
  real*8 gpi0
  real*8 gg(-n:n),sumg,sigma 		   !distr.de prob.de conexion
  real*8 ggee(-n:n),sumgee 	   !distr.de prob.de conexion
  real*8 rampanmda,rampanmdaie
  real*8 rescale,rescalexci,rescalinhib
  real*8 xx,yy,zz
  real*8 seed2

  integer ng,ngee				!ancho de gg (en sitios)
  integer ngei,ngie,ngii
  real*8 f,x
  real*8 rnd,ran3,pi,dpi,iesqk,gamma,iisqk,iisqk0,ran0
  real*8 ggubfs
  real*8 iesqk0,iesqk1,iesqk2
  real*8 iesqkon,iesqkoff,iisqkoff,iisqkon
  real*8 ample,ampli,ies,iis,theta
  real*8 thetacue,amplion,ampleon
  real*8 thetasac,amplisac,amplesac
  real*8 thetaoff,amplioff,ampleoff
  real*8 averc,avers,thetav,epsvec
  integer deltit,idir,anneal
  integer ireal,nreal
  character*80 nombfile,nombfile0
  character*5 type_ee,type_ei,type_ie,type_ii
  character*15 type_anneal
  character*94 r(100)

  real*8 tee2,tei2,tii2,tie2,teeampa2
  real*8 s1,s2,s3,dtf,dt2
  real*8 exee2,exei2,exii2,exie2
  real*8 expi2
  real*8 exeeampa2
  real*8 dt,eps9
  real*8 ser,set,sir,sit,spr,spt,dtf2
  real*8	snew,uu
  real*8 taue,taui

  real*8 coeff,coef,coef1,coef2,coefk
  real*8 vpspeenmda,vpspeeampa,vpspei,vpspie,vpspii
  integer je,ji,jp
  real*8 hepop,hipop
  common /distrib/gg,ng


  allocate (switche(ne),switchi(ni))
  allocate (se(0:n),si(0:n))             !estados de las poblaciones
  allocate (ue(nu),ui(nu)  )             !estados de las poblaciones
  allocate (semact(0:n),simact(0:n))
  allocate (seed(0:n),seem(0:n))
  allocate (averhe(n),averhi(n))
  allocate (averhe2(n),averhi2(n))
  allocate (averexcie(n),averexcie2(n))
  allocate (averinhibe(n),averinhibe2(n))
  allocate (averexcii(n),averexcii2(n))
  allocate (averinhibi(n),averinhibi2(n))
  allocate (mee(ne,nce),mei(ni,nce),mii(ni,nci))   !matrices de conexion
  allocate (mie(ne,nci))
  allocate (ncee(n),ncei(n),ncii(n),ncie(n))
  allocate (he(n),hi(n))    !campos internos
  allocate (ree(n),rie(n),rei(n),rii(n))
  allocate (hee(n),hie(n),hei(n),hii(n))
  allocate (heeampa(n),hieampa(n))
  allocate (nspike(n),nspiki(n))
  allocate (fee2(0:n),fei2(0:n),fii2(0:n))
  allocate (fie2(0:n))
  allocate (feeampa2(0:n),fieampa2(0:n))
  allocate (ke1(0:n),ki1(0:n),kp1(0:n))
  allocate (ke2(0:n),ki2(0:n),kp2(0:n))
  allocate (sebis(0:n),sibis(0:n),spbis(0:n))
  allocate (tspike(n),tspikeold(n))
  allocate (tspiki(n),tspikiold(n))
  allocate (tspikp(n),tspikpold(n))
  allocate (sumspe(n),sumsp2e(n),sumspi(n))
  allocate (sumsp2i(n))
  allocate (sumecv2(n),sumicv2(n),sumpcv2(n))
  allocate (dte(n),dteold(n),dti(n),dtiold(n),dtp(n),dtpold(n))
  allocate (fire(n),firi(n))
  allocate (cve(n),cvi(n),cvp(n))
  allocate (cv2e(n),cv2i(n),cv2p(n))
  allocate (cose(ne),cosi(ni),sine(ne),sini(ni))
  allocate (cose2(ne),cosi2(ni),sine2(ne),sini2(ni))
  allocate (thetae(ne),thetai(ni))
  allocate (c1e(nu),c1i(nu),s1e(nu),s1i(nu),m1e(nu),m1i(nu))
  allocate (c2e(nu),c2i(nu),s2e(nu),s2i(nu),m2e(nu),m2i(nu))




  eps9=1e-09
  epsvec=0.95
  uu=0.0001
  coeff=1.
  nreal=1

  !
  open(2,file='depeirtuningresicurcos.ini')
  read(2,*) iseed1
  !	read(2,*) iseed2
  read(2,*) seed2
  read(2,*) itmax
  read(2,*) dtf
  read(2,*) deltit
  read(2,*) gamma
  read(2,*) gee
  read(2,*) gei
  read(2,*) gii
  read(2,*) gie
  read(2,*) g1ee
  read(2,*) g1ei
  read(2,*) g1ii
  read(2,*) g1ie
  read(2,*) sigmaee
  read(2,*) sigmaei
  read(2,*) sigmaii
  read(2,*) sigmaie
  read(2,*) rescale
  read(2,*) rescalexci
  read(2,*) rescalinhib
  read(2,*) rampanmda
  read(2,*) rampanmdaie
  read(2,*) taue
  read(2,*) taui
  read(2,*) teeampa2
  read(2,*) tee2
  read(2,*) tei2
  read(2,*) tii2
  read(2,*) tie2
  read(2,*) ser
  read(2,*) set
  read(2,*) sir
  read(2,*) sit
  !read(2,*) ne
  !read(2,*) ni
  read(2,*) kkee
  read(2,*) kkei
  read(2,*) kkii
  read(2,*) kkie
  read(2,*) ie
  read(2,*) ii
  read(2,*) taufacee
  read(2,*) taurecee
  read(2,*) biguee
  read(2,*) tauact
  read(2,*) it1
  read(2,*) it2
  read(2,*) itrans
  read(2,*) it3
  read(2,*) it4
  read(2,*) it5
  !	read(2,*) iesqk1
  read(2,*) ieon
  read(2,*) iion
  !	read(2,*) iesqk2	
  read(2,*) ieoff
  read(2,*) iioff
  read(2,*) thetacue
  read(2,*) thetasac
  read(2,*) thetaoff
  read(2,*) ampleon
  read(2,*) amplion
  read(2,*) amplesac
  read(2,*) amplisac
  read(2,*) ampleoff
  read(2,*) amplioff
  read(2,*) anneal
  read(2,*) type_anneal
  read(2,*) nnn
  read(2,*) nombfile0
  close(2)

  geeampa=gee*rampanmda/(1.+rampanmda)
  gee=gee/(1.+rampanmda)
  gieampa=gie*rampanmdaie/(1+rampanmdaie)
  gie=gie/(1+rampanmdaie)
  !
  open(4,file='conect-depr-'//nombfile0)
  !
  kkee0=kkee
  kkei0=kkei
  kkie0=kkie
  kkii0=kkii


  coefk=1.
  kkee=kkee*coefk
  kkie=kkie*coefk
  kkei=kkei*coefk
  kkii=kkii*coefk

  coef2=1**2
  coef1=1**2
  coef=1.
  ie=ie*coef
  ii=ii*coef

  pi=dacos(-1.d0)
  dpi=2.d0*pi
  dtf2=dtf/2
  eps=1/taufac
  epsact=1/tauact
  iesqk=ie*dsqrt(1.d0*kkee)
  iisqk=ii*dsqrt(1.d0*kkii)
  iisqk0=iisqk
  iesqk0=iesqk
  !       iesqk1=iesqk1*sqrt(1.d0*kkee)*coef
  iesqkon=ieon*sqrt(1.d0*kkee)*coef
  if(anneal.eq.1) iesqk1=iesqkon
  iesqkoff=ieoff*sqrt(1.d0*kkee)*coef
  iisqkon=iion*sqrt(1.d0*kkii)*coef
  iisqkoff=iioff*sqrt(1.d0*kkii)*coef
  !       iesqk2=iesqk2*sqrt(1.d0*kkee)*coef
  !
  exee2=exp(-dtf/tee2)
  exeeampa2=exp(-dtf/teeampa2)
  exei2=exp(-dtf/tei2)
  exii2=exp(-dtf/tii2)
  exie2=exp(-dtf/tie2)

  umeps=exp(-dtf/taufac)
  umact=exp(-dtf/tauact)

  open(12,file='pr-depr-'//nombfile0)
  write(12,*) 'N,K=',n,kk
  write(12,*) 'iseed1=',iseed1
  !	write(12,*) 'iseed2=',iseed2
  write(12,*) 'seed2=',seed2
  write(12,*) 'itmax=',itmax
  write(12,*) 'ser=',ser
  write(12,*) 'set=',set
  write(12,*) 'sir=',sir
  write(12,*) 'sit=',sit

  write(12,*) 'gamma=',gamma
  write(12,*) 'gee=',gee
  write(12,*) 'geeampa=',geeampa
  write(12,*) 'gee (total)=',gee+geeampa
  write(12,*) 'gieampa=',gieampa
  write(12,*) 'gie (total)=',gie+gieampa
  write(12,*) 'gei=',gei
  write(12,*) 'gii=',gii
  write(12,*) 'gie=',gie
  write(12,*) 'g1ee =',g1ee
  write(12,*) 'g1ei=',g1ei
  write(12,*) 'g1ii=',g1ii
  write(12,*) 'g1ie=',g1ie
  write(12,*) 'sigmaee =',sigmaee
  write(12,*) 'sigmaei=',sigmaei
  write(12,*) 'sigmaii=',sigmaii
  write(12,*) 'sigmaie=',sigmaie
  write(12,*) 'rescale=',rescale
  write(12,*) 'rescalexci=',rescalexci
  write(12,*) 'rescalinhib=',rescalinhib
  write(12,*) 'rampanmda=',rampanmda
  write(12,*) 'rampanmdaie=',rampanmdaie

  write(12,*) 'ne=',ne
  write(12,*) 'ni=',ni
  write(12,*) 'kkee=',kkee0
  write(12,*) 'kkei=',kkei0
  write(12,*) 'kkii=',kkii0
  write(12,*) 'kkie=',kkie0

  write(12,*) 'ie=',ie
  write(12,*) 'it1,it2,itrans=',it1,it2,itrans
  write(12,*) 'it3,it4,it5=',it3,it4,it5
  write(12,*) 'ieon=',ieon
  write(12,*) 'iion=',iion
  write(12,*) 'ieoff=',ieoff
  write(12,*) 'iioff=',iioff
  write(12,*) 'thetacue=',thetacue
  write(12,*) 'thetaoff=',thetaoff
  write(12,*) 'ampleon=',ampleon
  write(12,*) 'amplion=',amplion
  write(12,*) 'amplesac=',amplesac
  write(12,*) 'amplisac=',amplisac
  write(12,*) 'ampleoff=',ampleoff
  write(12,*) 'amplioff=',amplioff
  write(12,*) 'taufacee=',taufacee
  write(12,*) 'taurecee=',taurecee
  write(12,*) 'biguee=',biguee
  write(12,*) 'tauact=',tauact
  write(12,*) 'tee2=',tee2
  write(12,*) 'teeampa2=',teeampa2
  write(12,*) 'tei2=',tei2
  write(12,*) 'tii2=',tii2
  write(12,*) 'tie2=',tie2
  write(12,*) 'taue=',taue
  write(12,*) 'taui=',taui
  write(12,*) 'nombfile=',nombfile0

  write(12,*) 'Connectivity for the simulated size:'
  write(12,*) 'kkee=',kkee
  write(12,*) 'kkei=',kkei
  write(12,*) 'kkii=',kkii
  write(12,*) 'kkie=',kkie
  write(12,*) 'anneal, type_anneal=',anneal,type_anneal
  call flush(12)


  gee=gee/sqrt((kkee))*gamma/tee2
  geeampa=geeampa/sqrt((kkee))*gamma/teeampa2
  gieampa=gieampa/sqrt((kkie))*gamma/teeampa2
  gei=gei/sqrt((kkei))*gamma/tei2
  gii=gii/sqrt((kkii))*gamma/tii2
  gie=gie/sqrt((kkie))*gamma/tie2

  gee=gee*rescale
  geeampa=geeampa*rescale
  gieampa=gieampa*rescale
  gei=gei*rescale
  gii=gii*rescale
  gie=gie*rescale
  ie=ie*rescale
  ii=ii*rescale

  gee=gee*rescalexci
  geeampa=geeampa*rescalexci
  gieampa=gieampa*rescalexci
  gei=gei*rescalinhib
  gii=gii*rescalinhib
  gie=gie*rescalexci

  !       gee=gee*exee2
  !       geeampa=geeampa*exeeampa2
  !       gieampa=gieampa*exeeampa2
  !       gei=gei*exei2
  !       gii=gii*exii2
  !       gie=gie*exie2
  !----------------------------------------------------------
  !	calculation PSP's
  write(12,*) 'Post-synaptic potentials in mV'
  call psp(gee,tee2,taue,set,vpspeenmda)
  write(12,*) 'PSP EE nmada=',vpspeenmda*biguee
  call psp(geeampa,teeampa2,taue,set,vpspeeampa)
  write(12,*) 'PSP EE ampa=',vpspeeampa*biguee
  call psp(gei,tei2,taue,set,vpspei)
  write(12,*) 'PSP EI =',vpspei
  call psp(gii,tii2,taui,sit,vpspii)
  write(12,*) 'PSP II =',vpspii
  call psp(gie,tie2,taui,sit,vpspie)
  write(12,*) 'PSP IE =',vpspie
  do i=1,ne
     thetae(i)=dpi*float(i-ne/2)/float(ne)
     cose(i)=dcos(thetae(i))
     sine(i)=dsin(thetae(i))
     cose2(i)=dcos(2.d0*thetae(i))
     sine2(i)=dsin(2.d0*thetae(i))
  end do
  do i=1,ni
     thetai(i)=dpi*float(i-ni/2)/float(ni)
     cosi(i)=dcos(thetai(i))
     sini(i)=dsin(thetai(i))
     cosi2(i)=dcos(2.d0*thetai(i))
     sini2(i)=dsin(2.d0*thetai(i))
  end do
  !	 --------- Matrices de conexion -------
  call crmatcos(iseed1,ne,ne,sigmaee,kkee,g1ee,nce,'ee',ncee,mee,mei,mie,mii)
  call crmatcos(iseed1,ne,ni,sigmaei,kkei,g1ei,nce,'ei',ncei,mee,mei,mie,mii)
  call crmatcos(iseed1,ni,ne,sigmaie,kkie,g1ie,nci,'ie',ncie,mee,mei,mie,mii)
  call crmatcos(iseed1,ni,ni,sigmaii,kkii,g1ii,nci,'ii',ncii,mee,mei,mie,mii)
  !	call crmat(iseed1,ne,ne,sigmaee,kkee,0.d0,nce,'ee',ncee,mee,mei,mie,mii)
  !       call crmat(iseed1,ne,ni,sigmaei,kkei,0.d0,nce,'ei',ncei,mee,mei,mie,mii)
  !       call crmat(iseed1,ni,ne,sigmaie,kkie,0.d0,nci,'ie',ncie,mee,mei,mie,mii)
  !       call crmat(iseed1,ni,ni,sigmaii,kkii,0.d0,nci,'ii',ncii,mee,mei,mie,mii)


  r(1)=nombfile0(1:nnn)//'-r01'
  ! r(1)=nombfile0(1:nnn)
  r(2)=nombfile0(1:nnn)//'-r02'
  r(3)=nombfile0(1:nnn)//'-r03'
  r(4)=nombfile0(1:nnn)//'-r04'
  r(5)=nombfile0(1:nnn)//'-r05'
  r(6)=nombfile0(1:nnn)//'-r06'
  r(7)=nombfile0(1:nnn)//'-r07'
  r(8)=nombfile0(1:nnn)//'-r08'
  r(9)=nombfile0(1:nnn)//'-r09'
  r(10)=nombfile0(1:nnn)//'-r10'

  do ireal=1,nreal
     open(3,file='ux-depr-'//r(ireal))
     open(26,file='spike-depr-'//r(ireal))
     open(27,file='spiki-depr-'//r(ireal))
     open(31,file='fire-depr-'//r(ireal))
     open(32,file='firi-depr-'//r(ireal))
     open(33,file='modul1-depr-'//r(ireal))
     open(34,file='modul2-depr-'//r(ireal))
     open(35,file='thetapop-depr-'//r(ireal))
     open(36,file='syncurneurE=1-depr-'//r(ireal))
     open(37,file='syncurneurE=6-depr-'//r(ireal))
     open(38,file='popaverage-depr-'//r(ireal))
     open(39,file='synaverage-depr-'//r(ireal))
     open(40,file='synaveragi-depr-'//r(ireal))
     open(41,file='inpute-depr-'//r(ireal))
     open(42,file='inputi-depr-'//r(ireal))
     !	---------- Condicion Inicial-----------

     se(0)=0.  !NOTA: el elemento 0 de todas las poblaciones es 
     si(0)=0.  !      siempre 0 y no se toca. Est�para cuando aparece 
     semact(0)=0.
     simact(0)=0.
     seem(0)=biguee
     seed(0)=1
     averc=0.
     avers=0
     thetav=0.
     do i=1,n
        fee2(i)=0
        fie2(i)=0
        feeampa2(i)=0
        fei2(i)=0
        fii2(i)=0

        he(i)=0 
        hi(i)=0 
        hee(i)=0 
        hei(i)=0 
        averhe(i)=0
        averhi(i)=0
        averhe2(i)=0
        averhi2(i)=0
        averexcie(i)=0
        averexcie2(i)=0
        averinhibe(i)=0
        averinhibe2(i)=0
        averexcii(i)=0
        averexcii2(i)=0
        averinhibi(i)=0
        averinhibi2(i)=0
        hii(i)=0 
        hie(i)=0 

        nspike(i)=0
        nspiki(i)=0
        tspikeold(i)=0
        tspikiold(i)=0
        sumspe(i)=0
        sumsp2e(i)=0
        sumspi(i)=0
        sumsp2i(i)=0

        sumecv2(i)=0
        sumicv2(i)=0

        dte(i)=0
        dteold(i)=0
        dti(i)=0
        dtiold(i)=0

     end do

     write(12,*) 'seed2=',seed2
     do i=1,n
        !	rnd=ran3(iseed2)
 	rnd=ggubfs(seed2)
  !	rnd=ran0(iseed2)
        se(i)=ser+rnd*(set-ser)
        !	rnd=ran3(iseed2)
 	rnd=ggubfs(iseed2)
  !	rnd=ran0(iseed2)
        si(i)=sir+rnd*(sit-sir)
        semact(i)=0
        simact(i)=0
        seed(i)=1.
        seem(i)=biguee
     end do

     !	----------- DINAMICA ----------

     !	open(3,file='ux-'//nombfile)
     !	open(26,file='spike-'//nombfile)
     !	open(27,file='spiki-'//nombfile)

     idir=1
     ample=0
     ampli=0

     do it=1,itmax

        ! 	if(it.gt.5000) ampl=0

	if(anneal.eq.0) then
           iesqk=iesqk0
           iisqk=iisqk0
           ample=0
           ampli=0

           if(it.gt.it1 .and. it.lt.it2) then
              !                       iesqk=iesqk1
              iesqk=iesqkon
              !                       iisqk=iesqk/sqrt(kkee)*sqrt(kkii)
              iisqk=iisqkon
              ample=ampleon
              ampli=amplion
              theta=thetacue
           end if
           if(it.gt.it3.and.it.lt.it4) then
              ample=amplesac
              ampli=amplisac
              theta=thetasac+thetav
           end if
           if(it.gt.it4.and.it.lt.it5) then
              ample=ampleoff
              ampli=amplioff
              theta=thetaoff+thetav
              !                       iesqk=iesqk2
              !                 iesqk=(-iesqkoff+iesqk0)*(it-it4)/float(it5-it4)
              !    .             +iesqkoff
              iesqk=iesqkoff
              !                       iisqk=iesqk2/sqrt(kkee)*sqrt(kkii)
              !                       iisqk=0.65*sqrt(kkii)
              !                 iisqk=(-iisqkoff+iisqk0)*(it-it4)/float(it5-it4)
              !    .  	   +iisqkoff
              iisqk=iisqkoff
              !                       iisqk=iisqk0+(3.*iesqk/sqrt(kkee)*sqrt(kkii)
              !    .    -iisqk0)*exp(-10*(float(it4+it3)/2-it)**2/float(it4-it3)**2)
              !                       iesqk=iesqk0+(0.375*sqrt(kkee)
              !    .  -iesqk0)*exp(-60.*(float(it4+it3)/2-it)**2/float(it4-it3)**2)
              !                       iesqk=2.5*iesqk
           end if

	else
           if(it.ge.itmax/2) idir=-1

           if(type_anneal.eq.'2dir') iesqk=iesqk+idir*(iesqk1-iesqk0)/(itmax/2)

           if(type_anneal.eq.'2direi') then
              iesqk=iesqk+idir*(iesqk1-iesqk0)/(itmax/2)
              iisqk=1.*iesqk/sqrt(kkee)*sqrt(kkii)    
           end if

           if(type_anneal.eq.'1dir') then
              iesqk=iesqk+(iesqk1-iesqk0)/itmax
              iisqk=1.*iesqk/sqrt(kkee)*sqrt(kkii)
           end if
           !	if(type_anneal.eq.'1dirgpe') gpec=gpe/float(itmax)*float(it)


	endif

	do k=1,ne
    !	if(fee2(k).gt.eps9) fee2(k)=fee2(k)*exee2  
    !	if(feeampa2(k).gt.eps9) feeampa2(k)=feeampa2(k)*exeeampa2  
    !	if(fieampa2(k).gt.eps9) fieampa2(k)=fieampa2(k)*exeeampa2  
    !	if(fie2(k).gt.eps9) fie2(k)=fie2(k)*exie2 
           if(semact(k).gt.eps9) semact(k)=semact(k)*umact
	end do
	do k=1,ni
    !	if(fei2(k).gt.eps9) fei2(k)=fei2(k)*exei2 
    !	if(fii2(k).gt.eps9) fii2(k)=fii2(k)*exii2  
           if(simact(k).gt.eps9) simact(k)=simact(k)*umact
	end do
 
	do k=1,ne
           ree(k)=0
           rei(k)=0
           hee(k)=hee(k)*exee2
           heeampa(k)=heeampa(k)*exeeampa2
	end do

	do k=1,ne
           if(switche(k).eq.1) then
              do j=1,ncee(k)
                 nn=mee(k,j)
                 hee(nn)=hee(nn)+fee2(k)	
                 heeampa(nn)=heeampa(nn)+feeampa2(k)	
              end do
           end if
	end do

 !       if(rampanmda.eq.0) go to 190
 !	do j=1,nce
 !	do k=1,ne
 !	heeampa(k)=heeampa(k)+feeampa2(mee(k,j))	
 !	end do
 !	end do

190     continue

	do k=1,ne
           hei(k)=hei(k)*exei2
	end do

	do k=1,ni
           if(switchi(k).eq.1) then
              do j=1,ncei(k)
                 nn=mei(k,j)
                 hei(nn)=hei(nn)+fei2(k)
              end do
           end if
	end do

        if(rampanmda.eq.0) then
           do k=1,ne
              he(k)=hee(k)*gee+hei(k)*gei
           end do
        else       
           do k=1,ne
              he(k)=hee(k)*gee+heeampa(k)*geeampa+hei(k)*gei
           end do
        end if

	if(it.gt.itrans .and. it.lt.it3) then
           do k=1,ne
              ies=iesqk*(1+ample-ample*cos(float(k)*dpi/float(ne)-theta*dpi/360.))
              averexcie(k)=averexcie(k)+hee(k)*gee+heeampa(k)*geeampa+ies
              averexcie2(k)=averexcie2(k)+(hee(k)*gee+heeampa(k)*geeampa+ies)**2
              averinhibe(k)=averinhibe(k)+hei(k)*gei
              averinhibe2(k)=averinhibe2(k)+(hei(k)*gei)**2
              averhe(k)=averhe(k)+he(k)+ies
              averhe2(k)=averhe2(k)+(he(k)+ies)**2
           end do
	end if


	do k=1,ni
           rie(k)=0
           rii(k)=0
           hie(k)=hie(k)*exie2
           hieampa(k)=hieampa(k)*exeeampa2
	end do
	do k=1,ne
           if(switche(k).eq.1) then
              do j=1,ncie(k)
                 nn=mie(k,j)
                 hie(nn)=hie(nn)+fie2(k)
                 hieampa(nn)=hieampa(nn)+fieampa2(k)
              end do
           end if
	end do

 !       do j=1,nce
 !	do k=1,ni
 !                hieampa(k)=hieampa(k)+fieampa2(mie(k,j))
 !	end do
 !	end do

	do k=1,ni
           hii(k)=hii(k)*exii2
	end do
	do k=1,ni
           if(switchi(k).eq.1) then
              do j=1,ncii(k)
                 nn=mii(k,j)
                 hii(nn)=hii(nn)+fii2(k) 
              end do
           end if
	end do


	do k=1,ni
           hi(k)=hie(k)*gie+hii(k)*gii+hieampa(k)*gieampa
	end do
        if(it.gt.itrans .and. it.lt.it3) then
           do k=1,ne
              iis=iisqk*(1+ampli-ampli*cos(float(k)*dpi/float(ni)-theta*dpi/360.))
              averexcii(k)=averexcii(k)+hie(k)*gie+hieampa(k)*gieampa+iis
              averexcii2(k)=averexcii2(k)+(hie(k)*gie+hieampa(k)*gieampa+iis)**2
              averinhibi(k)=averinhibi(k)+hii(k)*gii
              averinhibi2(k)=averinhibi2(k)+(hii(k)*gii)**2
              averhi(k)=averhi(k)+hi(k)+iis
              averhi2(k)=averhi2(k)+(hi(k)+iis)**2
           end do
        end if



        do k=1,ne
           !    ies=iesqk*(1-ampl*cos(float(k)*dpi/float(ne)-theta*dpi/360.))
           ies=iesqk*(1+ample-ample*cos(float(k)*dpi/float(ne)-theta*dpi/360.))
           sebis(k)=se(k)/(1+dtf/taue)+dtf*(he(k)+ies)/(taue+dtf)
        end do
        !	iis=iisqk
        do k=1,ni
           iis=iisqk*(1+ampli-ampli*cos(float(k)*dpi/float(ni)-theta*dpi/360.))
           sibis(k)=si(k)/(1+dtf/taui)+dtf*(hi(k)+iis)/(taui+dtf)
        end do


	do k=1,ne
           switche(k)=0
           if(sebis(k).gt.set) then
              switche(k)=1
              dt2=dtf*(sebis(k)-set)/(sebis(k)-se(k))
              snew=ser+(sebis(k)-set)*(1+dtf/taue*(-ser+se(k))/(sebis(k)-se(k)))
              tspike(k)=it*dtf
              dte(k)=tspike(k)-tspikeold(k)
              seem(k)=seem(k)*(1.-biguee)*dexp(-dte(k)/taufacee)+biguee
              xx=dexp(-dte(k)/taurecee)
              seed(k)=seed(k)*(1.-seem(k))*xx+1.-xx
              fee2(k)=exp(-dt2/tee2)*seem(k)*seed(k)
              !  		fee2(k)=fee2(k)+exp(-dt2/tee2)*seem(k)*seed(k)
              feeampa2(k)=exp(-dt2/teeampa2)*seem(k)*seed(k)
              fieampa2(k)=exp(-dt2/teeampa2)
              !  		feeampa2(k)=feeampa2(k)+exp(-dt2/teeampa2)*seem(k)*seed(k)
              !  		fieampa2(k)=fieampa2(k)+exp(-dt2/teeampa2)
              fie2(k)=exp(-dt2/tie2)
              !  		fie2(k)=fie2(k)+exp(-dt2/tie2)
              !  		fee2(k)=fee2(k)+(1.-dt2/tee2)
              !  		fie2(k)=fie2(k)+(1.-dt2/tie2)
              !	if(mod(k,4).eq.0) write(26,*) (it+(ireal-1)*itmax)*dtf*0.01,
              !	write(26,*) (it+(ireal-1)*itmax)*dtf*0.01,
              write(26,*) (it)*dtf*0.01,k
              call flush(26)
              se(k)=snew
              semact(k)=semact(k)+epsact*exp(-dt2/tauact)
              if(it.gt.itrans .and. it.lt.it3) then
                 nspike(k)=nspike(k)+1
                 sumspe(k)=sumspe(k)+dte(k)
                 sumsp2e(k)=sumsp2e(k)+dte(k)**2
                 sumecv2(k)=sumecv2(k)+ 2*abs(dte(k)-dteold(k))/(dte(k)+dteold(k))
              endif
              dteold(k)=dte(k)
              tspikeold(k)=tspike(k)
           else 
              se(k)=sebis(k)
           end if
	end do

	do k=1,ni
           switchi(k)=0
           if(sibis(k).gt.sit) then
              switchi(k)=1
              dt2=dtf*(sibis(k)-sit)/(sibis(k)-si(k))
              snew=sir+(sibis(k)-sit)*(1+dtf/taui*(-sir+si(k))/(sibis(k)-si(k)))
              fei2(k)=exp(-dt2/tei2)
              fii2(k)=exp(-dt2/tii2)
              !  		fei2(k)=fei2(k)+exp(-dt2/tei2)
              !  		fii2(k)=fii2(k)+exp(-dt2/tii2)
              !  		fei2(k)=fei2(k)+(1.-dt2/tei2)
              !  		fii2(k)=fii2(k)+(1.-dt2/tii2)
              !	 if(mod(k,2).eq.0) write(27,*) (it+(ireal-1)*itmax)*dtf*0.01,
              !	write(27,*) (it+(ireal-1)*itmax)*dtf*0.01,
              write(27,*) (it)*dtf*0.01,k
              call flush(27)
              si(k)=snew
              simact(k)=simact(k)+epsact*exp(-dt2/tauact)
              !		sim(k)=sim(k)+eps*(1-dt2/taufac)

              tspiki(k)=it*dtf
              dti(k)=tspiki(k)-tspikiold(k)
              if(it.gt.itrans .and. it.lt.it3) then
                 nspiki(k)=nspiki(k)+1
                 sumspi(k)=sumspi(k)+dti(k)
                 sumsp2i(k)=sumsp2i(k)+dti(k)**2
                 sumicv2(k)=sumicv2(k)+ 2*abs(dti(k)-dtiold(k))/(dti(k)+dtiold(k))
              endif
              tspikiold(k)=tspiki(k)
              dtiold(k)=dti(k)
           else 
              si(k)=sibis(k)
           end if
	end do

	if(mod(it,deltit).eq.10) then
           do iu=1,nu
              ue(iu)=0.
              ui(iu)=0.
              c2e(iu)=0.
              s2e(iu)=0.
              c2i(iu)=0.
              s2i(iu)=0.
              c1e(iu)=0.
              s1e(iu)=0.
              c1i(iu)=0.
              s1i(iu)=0.
              do i=1,ne/nu
                 ue(iu)=ue(iu)+100*semact((iu-1)*(ne/nu)+i)
                 c1e(iu)=c1e(iu)+100*semact((iu-1)*(ne/nu)+i)*cose(i)
                 s1e(iu)=s1e(iu)+100*semact((iu-1)*(ne/nu)+i)*sine(i)
                 c2e(iu)=c2e(iu)+100*semact((iu-1)*(ne/nu)+i)*cose2(i)
                 s2e(iu)=s2e(iu)+100*semact((iu-1)*(ne/nu)+i)*sine2(i)
              end do
              do i=1,ni/nu
                 ui(iu)=ui(iu)+100*simact((iu-1)*(ni/nu)+i)
                 c1i(iu)=c1i(iu)+100*simact((iu-1)*(ni/nu)+i)*cosi(i)
                 s1i(iu)=s1i(iu)+100*simact((iu-1)*(ni/nu)+i)*sini(i)
                 c2i(iu)=c2i(iu)+100*simact((iu-1)*(ni/nu)+i)*cosi2(i)
                 s2i(iu)=s2i(iu)+100*simact((iu-1)*(ni/nu)+i)*sini2(i)
              end do
              c2i(iu)=c2i(iu)/float(ni)
              s2i(iu)=s2i(iu)/float(ni)
              c2e(iu)=c2e(iu)/float(ne)
              s2e(iu)=s2e(iu)/float(ne)
              c1i(iu)=c1i(iu)/float(ni)
              s1i(iu)=s1i(iu)/float(ni)
              c1e(iu)=c1e(iu)/float(ne)
              s1e(iu)=s1e(iu)/float(ne)
              ue(iu)=ue(iu)/float(ne)
              ui(iu)=ui(iu)/float(ni)
              m1e(iu)=sqrt(c1e(iu)**2+s1e(iu)**2)/ue(iu)
              m1i(iu)=sqrt(c1i(iu)**2+s1i(iu)**2)/ui(iu)
              m2e(iu)=sqrt(c2e(iu)**2+s2e(iu)**2)/ue(iu)
              m2i(iu)=sqrt(c2i(iu)**2+s2i(iu)**2)/ui(iu)
              averc=averc*epsvec+c1e(iu)
              avers=avers*epsvec+s1e(iu)
              if(averc.eq.0.and.avers.eq.0) then
                 thetav=0
              else
                 thetav=acos(averc/sqrt(averc**2+avers**2))/pi*180
              end if
              if(avers.lt.0) thetav=-thetav

              if(anneal.eq.0) then
                 hepop=0
                 hipop=0
                 do i=1,ne
                    hepop=hepop+hee(i)*gee+iesqk
                 end do
                 hepop=hepop/float(ne)
                 do i=1,ni
                    hipop=hipop+hei(i)*gei
                 end do
                 hipop=hipop/float(ni)
                 write(3,*) sngl(it*dtf*0.01d0),sngl(ue(iu)*nu),sngl(ui(iu)*nu),sngl(iesqk),sngl(iisqk)
                 write(36,*) sngl(it*dtf*0.01d0),sngl(hee(1)*gee+iesqk),sngl(hei(1)*gei),sngl(he(1)+ies)
                 write(37,*) sngl(it*dtf*0.01d0),sngl(hepop),sngl(hipop),sngl(hepop+hipop)
                 write(38,*) sngl(it*dtf*0.01d0),sngl(hee(6)*gee+iesqk),sngl(hei(6)*gei),sngl(he(6)+ies)
                 write(33,*) sngl(it*dtf*0.01d0),sngl(m1e(iu)),sngl(m1i(iu)),sngl(iesqk),sngl(iisqk)
                 write(34,*) sngl(it*dtf*0.01d0),sngl(m2e(iu)),sngl(m2i(iu)),sngl(iesqk),sngl(iisqk)
                 write(35,*) sngl(it*dtf*0.01d0),thetav
                 !   write(45,*) sngl(it*dtf*0.01d0),sngl(se(1)),sngl(si(1))
              else
                 if(type_anneal.eq.'1dirgpe') then 
                    write(3,*) sngl(gpec),sngl((ue(iu)*nu)),sngl((ui(iu)*nu))
                    write(33,*) sngl(gpec),sngl(m1e(iu)),sngl(m1i(iu))
                    write(34,*) sngl(gpec),sngl(m2e(iu)),sngl(m2i(iu))
                 else
                    write(3,*) sngl(iesqk/sqrt(1.*kkee)),sngl(iisqk/sqrt(1.*kkii)),sngl((ue(iu)*nu)),sngl((ui(iu)*nu))
                    write(33,*) sngl(iesqk/sqrt(1.*kkee)),sngl(iisqk/sqrt(1.*kkii)),sngl(m1e(iu)),sngl(m1i(iu))
                    write(34,*) sngl(iesqk/sqrt(1.*kkee)),sngl(iisqk/sqrt(1.*kkii)),sngl(m2e(iu)),sngl(m2i(iu))
                 endif
              endif
           end do
           call flush(3)
           call flush(36)
           call flush(37)
           !	call flush(45)

           je=1
           ji=1
           jp=1
           !	write(44,*)  sngl(it*dtf*0.01d0),
           !    .   sngl(abs(ree(je)+iesqk+rei(je))/(ree(je)+iesqk-rei(je))),
           !    .   sngl(abs(rie(ji)+rii(ji)))/(rie(ji)-rii(ji))

           !       call flush(44)

	end if

     end do !it

     close(3)

     do i=1,ne
        fire(i)=nspike(i)/float(it3-itrans)/dtf
	averhe(i)=averhe(i)/float(it3-itrans)
	averhe2(i)=averhe2(i)/float(it3-itrans)
	write(39,*) i,averhe(i),dsqrt(averhe2(i)-averhe(i)**2)
	write(41,*) i,averexcie(i),averinhibe(i)
        if(nspike(i).gt.2) then
           sumspe(i)=sumspe(i)/nspike(i)
           sumsp2e(i)=sumsp2e(i)/nspike(i)
           cve(i)=sqrt(sumsp2e(i)-sumspe(i)**2)/sumspe(i)
           cv2e(i)=sumecv2(i)/nspike(i)
        endif
        if((nspike(i).gt.2)) then
           write(31,*) sngl(i*360.d0/ne),sngl(100*fire(i)), sngl(cve(i)),sngl(cv2e(i))
        else
           write(31,*) sngl(i*360.d0/ne),sngl(100.d0*fire(i)),0,0
        endif
     end do


     do i=1,ni
	averhi(i)=averhi(i)/float(it3-itrans)
	averhi2(i)=averhi2(i)/float(it3-itrans)
	write(40,*) i,averhi(i),dsqrt(averhi2(i)-averhi(i)**2)
	write(42,*) i,averexcii(i),averinhibi(i)
        firi(i)=nspiki(i)/float(it3-itrans)/dtf
        if(nspiki(i).gt.2) then
           sumspi(i)=sumspi(i)/nspiki(i)
           sumsp2i(i)=sumsp2i(i)/nspiki(i)
           cvi(i)=sqrt(sumsp2i(i)-sumspi(i)**2)/sumspi(i)
           cv2i(i)=sumicv2(i)/nspiki(i)
        endif
        if((nspiki(i).gt.2)) then
           write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),sngl(cvi(i)),sngl(cv2i(i))
        else
           write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),0,0
        endif
     end do

     write(12,*) 'ireal=',ireal
     call flush(12)
  end do !ireal


end PROGRAM INTFIRE  !---- FIN DEL PROGRAMA ----

subroutine crmat(iseed,n1,n2,sigma,k12,g0,nm,typcon,ncmjj,mee,mei,mie,mii)
  implicit none
  integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
  real*8 k12
  parameter(n=64000)
  parameter (kk=300)
  parameter(nc=2100,nce=2100,nci=800)
  parameter(ne=64000,ni=16000)
  real*8 gg(-n:n) 	!distr.de prob.de conexion
  integer ng		!ancho de gg (en sitios)
  integer i,j,il,iseed,ncm
  real*8  pp,rnd,ran3,ncmm,sumg,ran0
  real*8 sigma,g0,pi,dpi,p12,alpha
  real*8 theta1(n),theta2(n),deltatheta,sigmarad
  real*8 prob(0:ne/2)
  integer mee(ne,nce),mei(ni,nce),mii(ni,nci)      !matrices de conexion
  integer mie(ne,nci)
  integer ncmjj(n),ij
  character*2 typcon

  integer idum,itemp,jflone,jflmsk
  real ftemp
  equivalence (itemp,ftemp)
  data jflone /z'3f800000'/,jflmsk /z'007fffff'/

  dpi=2.d0*acos(-1.d0)
  pi=acos(-1.d0)
  !       -----------------------------------------------------------------------
  do i=1,n1
     theta1(i)=dpi*float(i)/float(n1)
  end do
  do i=1,n2
     theta2(i)=dpi*float(i)/float(n2)
  end do

  if(typcon.eq.'ee') then
     do i=1,nce
        do j=1,ne
           mee(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ei') then
     do i=1,nce
        do j=1,ni
           mei(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ie') then
     do i=1,nci
        do j=1,ne
           mie(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ii') then
     do i=1,nci
        do j=1,ni
           mii(j,i)=0
        end do
     end do
  end if

  sigmarad=sigma/360.*dpi

  sumg=0.
  alpha=1./2./sigmarad**2
  do i=1,n2/2
     prob(i)=g0+exp(-alpha*theta2(i)**2)
     !	prob(i)=1
     sumg=sumg+prob(i)
  end do
  sumg=1.+g0+2.*sumg
  sumg=k12/sumg


  !	open(4,file='ncmg.dat')

  if(typcon.eq.'ee') then
     do i=1,ne/2
        prob(i)=sumg*prob(i)	
     end do
     prob(0)=sumg*(1+g0)

     do j=1,ne
        ncm=0
        do i=1,ne
           ij=abs(i-j)
           if(ij.gt.n1/2) ij=n1-ij
           p12=prob(ij)

           !       idum=1664525*idum+1013904223
           !       itemp=ior(jflone,iand(jflmsk,idum))
           !       rnd=ftemp-1.



 	   rnd=ran3(iseed) 
     !	   rnd=ran0(iseed) 
           if (rnd.lt.p12) then
              ncm=ncm+1
              mee(j,ncm)=i	
              if (ncm.eq.nm) then
                 write(*,*) 'error en crmat: ncm=nc',i,il
                 stop
              end if
           end if
        end do
	ncmjj(j)=ncm
     end do

  else

     do j=1,n2
        ncm=0
        do i=1,n1
           deltatheta=dabs(theta1(i)-theta2(j))
           if(deltatheta.gt.pi) deltatheta=dpi-deltatheta
           p12=sumg*(g0+dexp(-alpha*deltatheta**2))
           !	p12=sumg
33         continue

           !	idum=1664525*idum+1013904223
           !       itemp=ior(jflone,iand(jflmsk,idum))
           !       rnd=ftemp-1.
           rnd=ran3(iseed)
           !          rnd=ran0(iseed)
           if (rnd.lt.p12) then
              ncm=ncm+1
              if(typcon.eq.'ei') mei(j,ncm)=i
              if(typcon.eq.'ie') mie(j,ncm)=i
              if(typcon.eq.'ii') mii(j,ncm)=i
              if (ncm.eq.nm) then
                 write(*,*) 'error en crmat: ncm=nc',i,il
                 stop
              end if
           end if
        end do
        ncmjj(j)=ncm
     end do

  end if

  do i=1,n2
     if(typcon.eq.'ee') write(4,*) i,mee(i,1)
     if(typcon.eq.'ei') write(4,*) i,mei(i,1)
     if(typcon.eq.'ie') write(4,*) i,mie(i,1)
     if(typcon.eq.'ii') write(4,*) i,mii(i,1)
  end do

  !	open(4,file='histg.dat')
  !	do i=1,nc
  !		write(4,*) i,mm(5000,i)
  !	end do
  !	close(4)

  return
end subroutine crmat


subroutine crmatcos(iseed,n1,n2,sigma,k12,g0,nm,typcon,ncmjj,mee,mei,mie,mii)
  implicit none
  integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
  real*8 k12
  parameter(n=64000)
  parameter (kk=300)
  parameter(nc=2100,nce=2100,nci=800)
  parameter(ne=64000,ni=16000)
  real*8 gg(-n:n) 	!distr.de prob.de conexion
  integer ng		!ancho de gg (en sitios)
  integer i,j,il,iseed,ncm
  real*8  pp,rnd,ran3,ncmm,sumg,ran0
  real*8 sigma,g0,pi,dpi,p12,alpha
  real*8 theta1(n),theta2(n),deltatheta,sigmarad
  real*8 prob(0:ne/2)
  integer mee(ne,nce),mei(ni,nce),mii(ni,nci)      !matrices de conexion
  integer mie(ne,nci)
  integer ncmjj(n),ij
  character*2 typcon

  integer idum,itemp,jflone,jflmsk
  real ftemp
  equivalence (itemp,ftemp)
  data jflone /z'3f800000'/,jflmsk /z'007fffff'/

  dpi=2.d0*acos(-1.d0)
  pi=acos(-1.d0)
  !       -----------------------------------------------------------------------
  do i=1,n1
     theta1(i)=dpi*float(i)/float(n1)
  end do
  do i=1,n2
     theta2(i)=dpi*float(i)/float(n2)
  end do

  if(typcon.eq.'ee') then
     do i=1,nce
        do j=1,ne
           mee(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ei') then
     do i=1,nce
        do j=1,ni
           mei(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ie') then
     do i=1,nci
        do j=1,ne
           mie(j,i)=0
        end do
     end do
  end if

  if(typcon.eq.'ii') then
     do i=1,nci
        do j=1,ni
           mii(j,i)=0
        end do
     end do
  end if

  sumg=0.
  do i=1,n2/2
     prob(i)=1.+g0*dcos(theta2(i))
     sumg=sumg+prob(i)
  end do
  sumg=1.+g0+2.*sumg
  sumg=k12/sumg


  !	open(4,file='ncmg.dat')

  if(typcon.eq.'ee') then
     do i=1,ne/2
        prob(i)=sumg*prob(i)	
     end do
     prob(0)=sumg*(1+g0)

     do j=1,ne
        ncm=0
        do i=1,ne
           ij=abs(i-j)
           if(ij.gt.n1/2) ij=n1-ij
           p12=prob(ij)

 	   rnd=ran3(iseed) 
     !	   rnd=ran0(iseed) 
           if (rnd.lt.p12) then
              ncm=ncm+1
              mee(j,ncm)=i	
              if (ncm.eq.nm) then
                 write(*,*) 'error en crmat: ncm=nc',i,il
                 stop
              end if
           end if
        end do
	ncmjj(j)=ncm
     end do

  else

     do j=1,n2
        ncm=0
        do i=1,n1
           deltatheta=dabs(theta1(i)-theta2(j))
           if(deltatheta.gt.pi) deltatheta=dpi-deltatheta
           p12=sumg*(1+g0*dcos(deltatheta))
33         continue

           !	idum=1664525*idum+1013904223
           !       itemp=ior(jflone,iand(jflmsk,idum))
           !       rnd=ftemp-1.
           rnd=ran3(iseed)
           !          rnd=ran0(iseed)
           if (rnd.lt.p12) then
              ncm=ncm+1
              if(typcon.eq.'ei') mei(j,ncm)=i
              if(typcon.eq.'ie') mie(j,ncm)=i
              if(typcon.eq.'ii') mii(j,ncm)=i
              if (ncm.eq.nm) then
                 write(*,*) 'error en crmat: ncm=nc',i,il
                 stop
              end if
           end if
        end do
        ncmjj(j)=ncm
     end do

  end if

  do i=1,n2
     if(typcon.eq.'ee') write(4,*) i,mee(i,1)
     if(typcon.eq.'ei') write(4,*) i,mei(i,1)
     if(typcon.eq.'ie') write(4,*) i,mie(i,1)
     if(typcon.eq.'ii') write(4,*) i,mii(i,1)
  end do

  !	open(4,file='histg.dat')
  !	do i=1,nc
  !		write(4,*) i,mm(5000,i)
  !	end do
  !	close(4)
  return
end subroutine crmatcos
!------------------------------------------------------------------
FUNCTION ran3(idum)
  INTEGER idum
  INTEGER MBIG,MSEED,MZ
  !     REAL MBIG,MSEED,MZ
  DOUBLE PRECISION ran3,FAC
  PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
  !     PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
  INTEGER i,iff,ii,inext,inextp,k
  INTEGER mj,mk,ma(55)
  !     REAL mj,mk,ma(55)
  SAVE iff,inext,inextp,ma
  DATA iff /0/
  if(idum.lt.0.or.iff.eq.0)then
     iff=1
     mj=MSEED-iabs(idum)
     mj=mod(mj,MBIG)
     ma(55)=mj
     mk=1
     do 11 i=1,54
        ii=mod(21*i,55)
        ma(ii)=mk
        mk=mj-mk
        if(mk.lt.MZ)mk=mk+MBIG
        mj=ma(ii)
11      continue
        do 13 k=1,4
           do 12 i=1,55
              ma(i)=ma(i)-ma(1+mod(i+30,55))
              if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
12            continue
13            continue
              inext=0
              inextp=31
              idum=1
           endif
           inext=inext+1
           if(inext.eq.56)inext=1
           inextp=inextp+1
           if(inextp.eq.56)inextp=1
           mj=ma(inext)-ma(inextp)
           if(mj.lt.MZ)mj=mj+MBIG
           ma(inext)=mj
           ran3=mj*FAC
           return
        END do
        !------------------------------------------------------
	subroutine psp(g,tau,tm,vt,vpsp)
          implicit none
          real*8 g,tau,tm,vt,vpsp,tp
          !
          tp=log(tm/tau)/(1/tau-1/tm)
          vpsp=g*exp(-tp/tau)
          vpsp=20.0*vpsp/vt
          !
          return
	end subroutine psp
 !-----------------------------------------------------
        FUNCTION ran0(idum)
          INTEGER idum,IA,IM,IQ,IR,MASK
          REAL ran0,AM
          PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,MASK=123459876)
          INTEGER k
          idum=ieor(idum,MASK)
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          ran0=AM*idum
          idum=ieor(idum,MASK)
          return
        END FUNCTION ran0
        !-----------------------------------------------------------
   	function ggubfs(seed)

          real*8 ggubfs
          real*8 seed
          real*8 d2p31m,d2p31
          data d2p31m /2147483647.d0/
          data d2p31  /2147483711.d0/
          !
          seed = dmod(16807.d0*seed,d2p31m)
          ggubfs=seed/d2p31
          return
        end function ggubfs
