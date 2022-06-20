	implicit none
	integer n,nc,kk,nu,nce,nci
	integer ne,ni,nn,nnn,sqne,sqni
	real*8 kkee,kkei,kkii,kkie,kklgn
	parameter(n=19600)
	parameter(nu=1) 
	parameter (kk=300) 
        parameter(nc=12000,nce=12000,nci=12000)
        parameter(ne=19600,ni=19600,sqne=140,sqni=140)
	integer i,j,k,iu,it,itmax
	integer it1
	integer i1,i2,i3,switche(ne),switchi(ni)
	real*8, allocatable :: se(:),si(:)		   
	real*8, allocatable :: xe(:,:),xi(:,:),empty(:,:)
	real*8, allocatable :: xoute(:,:),xouti(:,:)
	real*8, allocatable :: ue(:),ui(:)		   
	real*8, allocatable :: semact(:),simact(:)
	real*8, allocatable :: seed(:),seem(:)
	real*8, allocatable :: he(:),hi(:)    
	real*8, allocatable :: inpee(:),inpei(:)    
	real*8, allocatable :: inpeeclamp(:),inpeiclamp(:)    
	real*8, allocatable :: inpie(:),inpii(:)    
	real*8, allocatable :: inpieclamp(:),inpiiclamp(:)    
        real*8, allocatable :: heenmda(:),hienmda(:),hei(:),hii(:)
        real*8, allocatable :: avhee(:),avhie(:),avhei(:),avhii(:),avhetot(:),avhitot(:)
        real*8, allocatable :: avhee2(:),avhie2(:),avhei2(:),avhii2(:),avhetot2(:),avhitot2(:)
        real*8, allocatable :: avgee(:),avgie(:),avgei(:),avgii(:)
	real*8, allocatable :: heeampa(:),hieampa(:)
	real*8, allocatable :: feenmda2(:),fei2(:),fii2(:)
	real*8, allocatable :: fienmda2(:)
	real*8, allocatable :: feeampa2(:),fieampa2(:)
	real*8, allocatable :: ke1(:),ki1(:),kp1(:)
	real*8, allocatable :: ke2(:),ki2(:),kp2(:)
	real*8, allocatable :: sebis(:),sibis(:),spbis(:)
        real*8, allocatable :: tspike(:),tspikeold(:)
        real*8, allocatable :: tspiki(:),tspikiold(:)
        real*8, allocatable :: tspikp(:),tspikpold(:)
        real*8, allocatable :: sumspe(:),sumsp2e(:),sumspi(:)
	real*8, allocatable :: sumsp2i(:)
	real*8, allocatable :: sumecv2(:),sumicv2(:),sumpcv2(:)
	real*8, allocatable :: dte(:),dteold(:),dti(:),dtiold(:)
        real*8, allocatable :: fire(:),firi(:)
        real*8, allocatable :: fires(:),firis(:)
        real*8, allocatable :: firec(:),firic(:)
        real*8, allocatable :: cve(:),cvi(:),cvp(:)
        real*8, allocatable :: cv2e(:),cv2i(:),cv2p(:)
	real*8,allocatable :: cose(:),cosi(:),sine(:),sini(:)
	real*8,allocatable :: cose2(:),cosi2(:),sine2(:),sini2(:)
	real*8,allocatable :: thetae(:),thetai(:)
	real*8,allocatable :: thetare(:),thetari(:)
	real*8,allocatable :: phiosc1e(:),phiosc1i(:)
	real*8,allocatable :: phiosc2e(:),phiosc2i(:)
	real*8,allocatable :: phiosc3e(:),phiosc3i(:)
	real*8,allocatable :: c1e(:),c1i(:),s1e(:),s1i(:),m1e(:),m1i(:)
	real*8,allocatable :: c2e(:),c2i(:),s2e(:),s2i(:),m2e(:),m2i(:)
	real*8,allocatable :: gze(:),gzi(:)
	real*8,allocatable :: gke(:),gki(:)
	real*8,allocatable :: gnae(:),gnai(:)
	real*8,allocatable :: ksielgn(:),ksiilgn(:)
	real*8,allocatable :: ksieback(:),ksiiback(:)
	real*8,allocatable :: filtksielgn(:),filtksiilgn(:)
	real*8,allocatable :: filtksieback(:),filtksiiback(:)
	real*8,allocatable :: ratelgne(:),ratelgni(:)
	real*8,allocatable :: ratelgnef0(:),ratelgnif0(:)
	real*8,allocatable :: gffe(:),gffi(:)
	real*8,allocatable :: ies(:),iis(:)
	real*8,allocatable :: avervehyper(:),avervihyper(:)
	real*8,allocatable :: avervehypermodulc(:),avervihypermodulc(:)
	real*8,allocatable :: avervehypermoduls(:),avervihypermoduls(:)
	real*8,allocatable :: avervehyper2(:),avervihyper2(:)
	real*8,allocatable :: averve(:),avervi(:)
	real*8,allocatable :: averve2(:),avervi2(:)
	real*8,allocatable :: sigve(:),sigvi(:)
	real*8,allocatable :: randlgne(:),randlgni(:)
	real*8,allocatable :: zosc2e(:),zosc2i(:)
	real*8,allocatable :: zosc3e(:),zosc3i(:)
	real*8,allocatable :: zosc4e(:),zosc4i(:)
	real*8,allocatable :: randample0(:),randampli0(:)
	integer, allocatable :: mee(:,:),mei(:,:),mii(:,:),mie(:,:)   
	integer, allocatable :: ncee(:),ncei(:),ncii(:),ncie(:)
        integer, allocatable :: nspike(:),nspiki(:)
	real*8 eps,umeps,taufac,taurec,bigu  
	real*8 taufacee,taurecee,biguee
	real*8 gee,gie,geeampa,gei,gii,geenmda,ie,ii   
	real*8 ginputee,ginputei,ratelgn0,ratelgn1
	real*8 sigmaee,sigmaei,sigmaii,sigmaie   
	real*8 tauact,epsact,umact
	real*8 gieampa,gienmda
	real*8 gg(-n:n),sumg,sigma 		  
	real*8 ggee(-n:n),sumgee 	   
        real*8 rampanmda,rampanmdaie
	real*8 rescale,rescalexci,rescalinhib
	real*8 xx,yy,zz,t
	real*8 r1,r2,x1,x2
	real*8 curhyper
	integer nhyper,nqx314,qx314
	integer ng,ngee			
        integer ngei,ngie,ngii
	real*8 rnd,pi,dpi
	real*8 ample,ampli,theta,ample0,ampli0
	real*8 thetacue
	real*8 averc,avers,thetav,epsvec
	real*8 seminit,ggubfs,semadapt,semarch
	real*8 npe,npi
	integer deltit
	integer ireal,nreal
	character*99 nombfile,nombfile0,typelgn
        character*5 type_ee,type_ei,type_ie,type_ii
	character*110 r(100)
        real*8 teenmda2,tei2,tii2,tienmda2,teeampa2,tieampa2
	real*8 s1,s2,s3,dtf,dt2
	real*8 exeenmda2,exei2,exii2,exienmda2,exieampa2,exeeampa2
	real*8 dt,eps9
	real*8 spr,spt,dtf2
	real*8 snew,uu
        real*8 vpspeenmda,vpspeeampa,vpspei,vpspie,vpspii
        real*8 gle,gnaeav,gkeav,gae,gnape,phie,dgki,dgnae
        real*8 gli,gnaiav,gkiav,gai,gnapi,phii,dgke,dgnai
	real*8 tauz,taua,vna,vk,vl,va,phi
	real*8 thresh,vsyne,vsyni
	real*8 gzemin,gzemax
	real*8 gzimin,gzimax
	real*8 gbacke,gbacki,ratebacke,ratebacki
	real*8 gbacke0,gbacki0
	real*8 glgne,glgni
	real*8 contrast,alphacondcur
	real*8 vclamp
	real*8 allsigve,sigveall,averveall,averviall,veall,viall
	real*8 allsigvi,sigviall,averveall2,averviall2
	real*8 averveallt,averveall2t
	real*8 averviallt,averviall2t
	real*8 khie,khii
	real*8 coefilt0,coefilt1
	real*8 mu,f,omega,omt
	real*8 totspike,totspiki
        integer n1e,n2e,n3e
        integer n1i,n2i,n3i
	common /distrib/gg,ng
	common/cb/ tauz,taua,vna,vk,vl,va
allocate (se(0:n),si(0:n))             
allocate (xe(5,n),xi(5,n),empty(5,n))
allocate (xoute(5,n),xouti(5,n))
allocate (ue(nu),ui(nu))               
allocate (semact(0:n),simact(0:n))
allocate (seed(0:n),seem(0:n))
allocate (nspike(n),nspiki(n))
allocate (mee(nce,ne),mei(nce,ni),mii(nci,ni),mie(nci,ne))
allocate (ncee(n),ncei(n),ncii(n),ncie(n))
allocate (he(n),hi(n))    
allocate (inpee(n),inpei(n))    
allocate (inpeeclamp(n),inpeiclamp(n))    
allocate (inpie(n),inpii(n))    
allocate (inpieclamp(n),inpiiclamp(n))    
allocate (heenmda(n),hienmda(n),hei(n),hii(n))
allocate (avhee(n),avhie(n),avhei(n),avhii(n),avhetot(n),avhitot(n))
allocate (avhee2(n),avhie2(n),avhei2(n),avhii2(n),avhetot2(n),avhitot2(n))
allocate (heeampa(n),hieampa(n))
allocate (avgee(n),avgei(n),avgie(n),avgii(n))
allocate (feenmda2(0:n),fei2(0:n),fii2(0:n),fienmda2(0:n),feeampa2(0:n),fieampa2(0:n))
allocate (ke1(0:n),ki1(0:n),kp1(0:n))
allocate (ke2(0:n),ki2(0:n),kp2(0:n))
allocate (sebis(0:n),sibis(0:n),spbis(0:n))
allocate (tspike(n),tspikeold(n),tspiki(n),tspikiold(n))
allocate (tspikp(n),tspikpold(n),sumspe(n),sumsp2e(n),sumspi(n))
allocate (sumsp2i(n),sumecv2(n),sumicv2(n),sumpcv2(n))
allocate (dte(n),dteold(n),dti(n),dtiold(n),fire(n),firi(n))
allocate (fires(n),firis(n),firec(n),firic(n))
allocate (cve(n),cvi(n),cvp(n),cv2e(n),cv2i(n),cv2p(n))
allocate (cose(n),cosi(n),sine(n),sini(n))
allocate (cose2(n),cosi2(n),sine2(n),sini2(n),thetae(n),thetai(n))
allocate (thetare(n),thetari(n))
allocate (phiosc1e(n),phiosc1i(n))
allocate (phiosc2e(n),phiosc2i(n))
allocate (phiosc3e(n),phiosc3i(n))
allocate (c1e(nu),c1i(nu),s1e(nu),s1i(nu),m1e(nu),m1i(nu))
allocate (c2e(nu),c2i(nu),s2e(nu),s2i(nu),m2e(nu),m2i(nu))
allocate (ratelgne(n),ratelgni(n))
allocate (ratelgnef0(n),ratelgnif0(n))
allocate (gffe(n),gffi(n))
allocate (ksielgn(n),ksiilgn(n),ksieback(n),ksiiback(n))
allocate (filtksielgn(n),filtksiilgn(n),filtksieback(n),filtksiiback(n))
allocate (ies(n),iis(n))
allocate (avervehyper(n),avervihyper(n))
allocate (avervehypermodulc(n),avervihypermodulc(n))
allocate (avervehypermoduls(n),avervihypermoduls(n))
allocate (avervehyper2(n),avervihyper2(n))
allocate (averve(n),avervi(n))
allocate (averve2(n),avervi2(n))
allocate (sigve(n),sigvi(n))
allocate (randlgne(n),randlgni(n))
allocate (zosc2e(n),zosc2i(n))
allocate (zosc3e(n),zosc3i(n))
allocate (zosc4e(n),zosc4i(n))
allocate (randample0(n),randampli0(n))

	allocate (gze(n),gzi(n))
	allocate (gke(n),gki(n))
	allocate (gnae(n),gnai(n))

	eps9=1e-09
	epsvec=0.95
	uu=0.0001
	nreal=1

	open(2,file='nomapl42dosc.ini')
	read(2,*) semarch
	read(2,*) seminit
	read(2,*) itmax
	read(2,*) dtf
	read(2,*) deltit
	read(2,*) gee
	read(2,*) gei
	read(2,*) gii
	read(2,*) gie
	read(2,*) glgne
	read(2,*) glgni
	read(2,*) gbacke0
	read(2,*) gbacki0
	read(2,*) sigmaee
	read(2,*) sigmaei
	read(2,*) sigmaii
	read(2,*) sigmaie
        read(2,*) rampanmda
        read(2,*) rampanmdaie
        read(2,*) teeampa2
        read(2,*) teenmda2
        read(2,*) tei2
        read(2,*) tii2
        read(2,*) tieampa2
        read(2,*) tienmda2
	read(2,*) kkee
	read(2,*) kkei
	read(2,*) kkii
	read(2,*) kkie
	read(2,*) kklgn
	read(2,*) tauact
	read(2,*) it1
	read(2,*) alphacondcur
	read(2,*) contrast
	read(2,*) ratelgn0
	read(2,*) ratelgn1
	read(2,*) ratebacke
	read(2,*) ratebacki
	read(2,*) thetacue
	read(2,*) mu
	read(2,*) omega
	read(2,*) curhyper
	read(2,*) nhyper
	read(2,*) qx314
	read(2,*) nqx314
	read(2,*) typelgn
	read(2,*) ample
	read(2,*) ampli
	read(2,*) ample0
	read(2,*) ampli0
	read(2,*) npe
	read(2,*) npi
        read(2,*) n1e
        read(2,*) n2e
        read(2,*) n3e
        read(2,*) n1i
        read(2,*) n2i
        read(2,*) n3i
	read(2,*) nnn
	read(2,*) nombfile0
	close(2)


        gnaeav=100.0
	dgnae=0.0
        gkeav=40.
	dgke=0.
        gae=0.
        gnape=0.0
	gzemin=0.5
	gzemax=0.5
	gle=0.05
	phie=10.

        gnaiav=100.0
	dgnai=0.0
        gkiav=40.
	dgki=0.
        gai=0.
        gnapi=0.0
	gzimin=0.
	gzimax=0.
	gli=0.1
	phii=10.

        vna=55.
        vl=-65
        vk=-80
        va=-80.
	taua=20
	tauz=60
        thresh=-20.

	vsyni=-80
	vsyne=0
	vclamp=-40

        geeampa=gee*rampanmda/(1.+rampanmda)
        geenmda=gee/(1.+rampanmda)
        gieampa=gie*rampanmdaie/(1+rampanmdaie)
        gienmda=gie/(1+rampanmdaie)
	open(4,file='nomapl42dosc-conect-'//nombfile0)

        pi=dacos(-1.d0)
        dpi=2.d0*pi
        dtf2=dtf/2
        epsact=1/tauact
	exeenmda2=exp(-dtf/teenmda2)
	exeeampa2=exp(-dtf/teeampa2)
	exei2=exp(-dtf/tei2)
	exii2=exp(-dtf/tii2)
	exieampa2=exp(-dtf/tieampa2)
	exienmda2=exp(-dtf/tienmda2)
	coefilt0=1.-dtf/teeampa2
	coefilt1=dsqrt(dtf)/teeampa2

	umact=exp(-dtf/tauact)

        semadapt=209927.d0
        do i=1,n
        gze(i)=gzemin+(gzemax-gzemin)*ggubfs(semadapt)
        gzi(i)=gzimin+(gzimax-gzimin)*ggubfs(semadapt)
        gke(i)=gkeav*(1-2.*dgke*(ggubfs(semadapt)-0.5))
        gki(i)=gkiav*(1-2.*dgki*(ggubfs(semadapt)-0.5))
        gnae(i)=gnaeav*(1-2.*dgnae*(ggubfs(semadapt)-0.5))
        gnai(i)=gnaiav*(1-2.*dgnai*(ggubfs(semadapt)-0.5))
       	end do
	if(qx314.eq.1) then
	do i=1,nqx314
	gnae(i)=0
	gnai(i)=0
	end do
	end if

 	r(1)=nombfile0(1:nnn)//'-r10'

	open(12,file='nomapl42dosc-pr-'//r(1))
	write(12,*) 'N,K=',n,kk
	write(12,*) 'seminit=',seminit
	write(12,*) 'semarch=',semarch
	write(12,*) 'itmax=',itmax
	write(12,*) 'dt=',dtf
	write(12,*) 'geenmda=',geenmda
	write(12,*) 'geeampa=',geeampa
	write(12,*) 'gee (total)=',geenmda+geeampa
	write(12,*) 'gienmda=',gienmda
	write(12,*) 'gieampa=',gieampa
	write(12,*) 'gie (total)=',gienmda+gieampa
	write(12,*) 'gei=',gei
	write(12,*) 'gii=',gii
	write(12,*) 'sigmaee =',sigmaee
	write(12,*) 'sigmaei=',sigmaei
	write(12,*) 'sigmaii=',sigmaii
	write(12,*) 'sigmaie=',sigmaie
	write(12,*) 'rampanmda=',rampanmda
	write(12,*) 'rampanmdaie=',rampanmdaie
        write(12,*) 'ne=',ne
	write(12,*) 'ni=',ni
	write(12,*) 'gbacke0=',gbacke0
	write(12,*) 'gbacki0=',gbacki0
	write(12,*) 'it1=',it1
        write(12,*) 'glgne=',glgne
        write(12,*) 'glgni=',glgni
	write(12,*) 'alphacondcur=',alphacondcur
	write(12,*) 'contrast=',contrast
	write(12,*) 'ratelgn0,ratelgn1=',ratelgn0,ratelgn1
	write(12,*) 'ratebacke=', ratebacke
	write(12,*) 'ratebacki=', ratebacki
	write(12,*) 'thetacue=',thetacue
	write(12,*) 'mu=',mu
	write(12,*) 'f=',f
	write(12,*) 'omega=',omega
	write(12,*) 'typelgn=',typelgn
	write(12,*) 'ample=',ample
	write(12,*) 'ampli=',ampli
	write(12,*) 'npe=',npe
	write(12,*) 'npi=',npi
	write(12,*) 'ample0=',ample0
	write(12,*) 'ampli0=',ampli0
	write(12,*) 'tauact=',tauact
	write(12,*) 'teenmda2=',teenmda2
        write(12,*) 'teeampa2=',teeampa2
	write(12,*) 'tei2=',tei2
	write(12,*) 'tii2=',tii2
	write(12,*) 'tienmda2=',tienmda2
	write(12,*) 'tieampa2=',tieampa2
	write(12,*) 'nombfile=',nombfile0

	write(12,*) 'Connectivity for the simulated size:'
	write(12,*) 'kkee=',kkee
        write(12,*) 'kkei=',kkei
        write(12,*) 'kkii=',kkii
        write(12,*) 'kkie=',kkie
        write(12,*) 'kklgn=',kklgn
        write(12,*) 'gnaeav=',gnaeav
        write(12,*) 'dgnae=',dgnae
        write(12,*) 'dgke=',dgke
        write(12,*) 'gkeav=',gkeav
        write(12,*) 'gae=',gae
        write(12,*) 'gnape=',gnape
        write(12,*) 'gzemin=',gzemin
        write(12,*) 'gzemax=',gzemax
        write(12,*) 'gle=',gle
	write(12,*)  'phie=',phie
        write(12,*) 'gnaiav=',gnaiav
        write(12,*) 'dggnai=',dgnai
        write(12,*) 'dgki=',dgki
        write(12,*) 'gkiav=',gkiav
        write(12,*) 'gai=',gai
        write(12,*) 'gnapi=',gnapi
        write(12,*) 'gzimin=',gzimin
        write(12,*) 'gzimax=',gzimax
        write(12,*) 'gli=',gli
	write(12,*)  'phii=',phii
        write(12,*) 'vna=',vna
        write(12,*) 'vl=',vl
        write(12,*) 'vk=',vk
        write(12,*) 'va=',va
        write(12,*) 'vsyne=',vsyne
        write(12,*) 'vsyni=',vsyni
	write(12,*) 'vclamp=',vclamp
	write(12,*) 'curhyper=',curhyper
	write(12,*) 'nhyper=',nhyper
	write(12,*) 'qx314=',qx314
	write(12,*) 'nqx314=',nqx314
	write(12,*) 'n1e=',n1e
	write(12,*) 'n2e=',n2e
	write(12,*) 'n3e=',n3e
	write(12,*) 'n1i=',n1i
	write(12,*) 'n2i=',n2i
	write(12,*) 'n3i=',n3i
        write(12,*) 'taua=',taua
        write(12,*) 'tauz=',tauz
	write(12,*) 'seminit=',seminit
	write(12,*) 'semadapt=',semadapt

	call flush(12)


	geenmda=geenmda/sqrt((kkee))/teenmda2
	geeampa=geeampa/sqrt((kkee))/teeampa2
	gienmda=gienmda/sqrt((kkie))/tienmda2
	gieampa=gieampa/sqrt((kkie))/tieampa2
	gei=gei/sqrt((kkei))/tei2
	gii=gii/sqrt((kkii))/tii2
	gbacke=gbacke0*sqrt(kkee)
	gbacki=gbacki0*sqrt(kkee)
	glgne=glgne*sqrt(kklgn)
	glgni=glgni*sqrt(kklgn)
	mu=mu/dsqrt(2.d0)
	f=f/dsqrt(2.d0)

!----------------------------------------------------------
	do i=1,ne
 	thetae(i)=pi*float(i-ne/2)/float(ne)
 	thetare(i)=pi*(ggubfs(semarch)-0.5)
	end do
	do i=1,ni
 	thetai(i)=pi*float(i-ni/2)/float(ni)
 	thetari(i)=pi*(ggubfs(semarch)-0.5)
	end do
	cose=dcos(2.*thetae)
	sine=dsin(2.d0*thetae)
	cose2=dcos(4.d0*thetae)
	sine2=dsin(4.d0*thetae)
	cosi=dcos(2.d0*thetai)
	sini=dsin(2.d0*thetai)
	cosi2=dcos(4.d0*thetai)
	sini2=dsin(4.d0*thetai)
sigmaee=1./sigmaee
sigmaei=1./sigmaei
sigmaie=1./sigmaie
sigmaii=1./sigmaii

call crmat(semarch,ne,ne,sqne,sqne,sigmaee,kkee,0.d0,nce,'ee',ncee,mee,mei,mie,mii)
call crmat(semarch,ne,ni,sqne,sqni,sigmaei,kkei,0.d0,nce,'ei',ncei,mee,mei,mie,mii)
call crmat(semarch,ni,ne,sqni,sqne,sigmaie,kkie,0.d0,nci,'ie',ncie,mee,mei,mie,mii)
call crmat(semarch,ni,ni,sqni,sqni,sigmaii,kkii,0.d0,nci,'ii',ncii,mee,mei,mie,mii)
do i=1,n
if(typelgn.eq.'rand') then
randlgne(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
randlgni(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
end if
zosc2e(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
zosc2i(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
zosc3e(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
zosc3i(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
zosc4e(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
zosc4i(i)=1./sqrt(2.*pi)*sqrt(-2.*log(ggubfs(semarch)))
phiosc1e(i)=dpi*ggubfs(semarch)
phiosc2e(i)=dpi*ggubfs(semarch)
phiosc3e(i)=dpi*ggubfs(semarch)
phiosc1i(i)=dpi*ggubfs(semarch)
phiosc2i(i)=dpi*ggubfs(semarch)
phiosc3i(i)=dpi*ggubfs(semarch)
end do


!	r(1)=nombfile0(1:nnn)//'-r10'
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
 	open(3,file='nomapl42dosc-ux-'//r(ireal))
        open(31,file='nomapl42dosc-fire-'//r(ireal))
        open(32,file='nomapl42dosc-firi-'//r(ireal))
        open(30,file='nomapl42dosc-allfir-'//r(ireal))
!       open(36,file='nomapl42dosc-ve-'//r(ireal))
!       open(37,file='nomapl42dosc-vi-'//r(ireal))
!        open(38,file='nomapl42dosc-he-'//r(ireal))
!        open(39,file='nomapl42dosc-hi-'//r(ireal))
	open(56,file='nomapl42dosc-avinpe-'//r(ireal))
	open(57,file='nomapl42dosc-avinpi-'//r(ireal))
	open(66,file='nomapl42dosc-sdinpe-'//r(ireal))
	open(67,file='nomapl42dosc-sdinpi-'//r(ireal))
	open(58,file='nomapl42dosc-avige-'//r(ireal))
	open(59,file='nomapl42dosc-avigi-'//r(ireal))
	open(60,file='nomapl42dosc-avvehyper-'//r(ireal))
	open(80,file='nomapl42dosc-avvehypermodul-'//r(ireal))
	open(61,file='nomapl42dosc-avvihyper-'//r(ireal))
	open(81,file='nomapl42dosc-avvihypermodul-'//r(ireal))
        open(71,file='nomapl42dosc-modfire-'//r(ireal))
        open(72,file='nomapl42dosc-modfiri-'//r(ireal))
	open(96,file='nomapl42dosc-khi-'//r(ireal))
	open(27,file='nomapl42dosc-spike-'//r(ireal))
	open(28,file='nomapl42dosc-spiki-'//r(ireal))
	open(29,file='nomapl42dosc-totrate-'//r(ireal))

se(0)=0.  !NOTA: el elemento 0 de todas las poblaciones es 
si(0)=0.  !      siempre 0 y no se toca. Est�para cuando aparece 
semact(0)=0.
simact(0)=0.
averc=0.
avers=0
thetav=0.
feenmda2=0
fienmda2=0
feeampa2=0
fieampa2=0
fei2=0
fii2=0
he=0 
hi=0 
heeampa=0 
heenmda=0 
hei=0 
hii=0 
hieampa=0 
hienmda=0 
avhee=0 
avhei=0 
avhii=0 
avhie=0 
avhetot=0
avhitot=0
avhee2=0 
avhei2=0 
avhii2=0 
avhie2=0 
avhetot2=0
avhitot2=0
avgee=0 
avgei=0 
avgii=0 
avgie=0 
nspike=0
nspiki=0
fires=0
firis=0
firec=0
firic=0
tspikeold=0
tspikiold=0
sumspe=0
sumsp2e=0
sumspi=0
sumsp2i=0
sumecv2=0
sumicv2=0
dte=0
dteold=0
dti=0
dtiold=0
semact=0
simact=0
he=0
hi=0
ksielgn=0
ksieback=0
ksiilgn=0
ksiiback=0
filtksielgn=0
filtksieback=0
filtksiilgn=0
filtksiiback=0
ies=0
iis=0
avervehyper=0
avervihyper=0
avervehypermoduls=0
avervehypermodulc=0
avervihypermoduls=0
avervihypermodulc=0
avervehyper2=0
avervihyper2=0
averve=0
avervi=0
averve2=0
avervi2=0
averveall=0
averviall=0
averveall2=0
averviall2=0
do i=1,nhyper
ies(i)=curhyper
iis(i)=curhyper
end do
totspike=0.d0
totspiki=0.d0

thetacue=thetacue/180.*pi

do i=1,ne
xe(1,i)=-70.+100.*ggubfs(seminit)
xe(2,i)=ggubfs(seminit)
xe(3,i)=ggubfs(seminit)
xe(4,i)=ggubfs(seminit)
xe(5,i)=ggubfs(seminit)
end do
!call devt(xe,empty,he,ne,gle,gnae,gke,phie,gae,gnape,gze)

do i=1,ni
xi(1,i)=-70.+100.*ggubfs(seminit)
xi(2,i)=ggubfs(seminit)
xi(3,i)=ggubfs(seminit)
xi(4,i)=ggubfs(seminit)
xi(5,i)=ggubfs(seminit)
end do


!call devt(xi,empty,hi,ni,gli,gnai,gki,phii,gai,gnapi,gzi)


do i=1,n
x1=ggubfs(semarch)
x2=ggubfs(semarch)
randample0(i)=dsqrt(-2.d0*dlog(x1))*dcos(dpi*x2)/sqrt(kklgn)
randampli0(i)=dsqrt(-2.d0*dlog(x1))*dsin(dpi*x2)/sqrt(kklgn)
!NOTE THAT IN PREVIOUS VERSIONS  : randampli0(i)=dsqrt(-2.d0*dlog(x2))*dcos(dpi*x1)/sqrt(kklgn)
!SO THAT NO CORRELATION BETWEEN E AND I 
end do

ratelgne=ratelgn0*(1.+randample0)
ratelgni=ratelgn0*(1.+randampli0)

do it=1,itmax
if(it.eq.it1) then
ratelgnef0=ratelgn0*(1.+randample0)+ratelgn1*dlog(1.+contrast)/dlog(10.d0)*(ample0+randample0+ample*randlgne/sqrt(kklgn)*dcos(2.d0*(thetare-thetacue)))
ratelgnif0=ratelgn0*(1.+randampli0)+ratelgn1*dlog(1.+contrast)/dlog(10.d0)*(ampli0+randampli0+ampli*randlgni/sqrt(kklgn)*dcos(2.d0*(thetari-thetacue)))
end if

if(it.gt.it1) then
omt=float(it-it1)*dtf*omega
ratelgne=ratelgnef0+ratelgn1*dlog(1.+contrast)/dlog(10.d0)/sqrt(kklgn)*(mu*zosc2e*dcos(omt-phiosc1e)+ample*mu/2.*(zosc3e*dcos(2.*thetacue+omt-phiosc2e)+zosc4e*dcos(2.*thetare-omt-phiosc3e)))
ratelgni=ratelgnif0+ratelgn1*dlog(1.+contrast)/dlog(10.d0)/sqrt(kklgn)*(mu*zosc2i*dcos(omt-phiosc1i)+ampli*mu/2.*(zosc3i*dcos(2.*thetacue+omt-phiosc2i)+zosc4i*dcos(2.*thetare-omt-phiosc3i)))
end if

!go to 20

do i=1,n
x1=ggubfs(seminit)
x2=ggubfs(seminit)
ksieback(i)=dsqrt(-2.d0*dlog(x1))*dcos(dpi*x2)
!THE FOLLOWING LINE WAS CORRECTED ON NOV 27 MORNING : it was originally: ksiiback(i)=dsqrt(-2.d0*dlog(x2))*dcos(dpi*x1)
ksiiback(i)=dsqrt(-2.d0*dlog(x1))*dsin(dpi*x2)
end do
filtksieback=filtksieback*coefilt0+ksieback*coefilt1
filtksiiback=filtksiiback*coefilt0+ksiiback*coefilt1
write(98,*) filtksieback(1),filtksiiback(1)
do i=1,n
x1=ggubfs(seminit)
x2=ggubfs(seminit)
ksielgn(i)=dsqrt(-2.d0*dlog(x1))*dcos(dpi*x2)
!THE FOLLOWING LINE WAS CORRECTED ON NOV 27 MORNING : it was originally: ksiiback(i)=dsqrt(-2.d0*dlog(x2))*dcos(dpi*x1)
ksiilgn(i)=dsqrt(-2.d0*dlog(x1))*dsin(dpi*x2)
end do
filtksielgn=filtksielgn*coefilt0+ksielgn*coefilt1
filtksiilgn=filtksiilgn*coefilt0+ksiilgn*coefilt1


!20 continue

gffe=gbacke*(ratebacke+dsqrt(ratebacke/kkee)*filtksieback)+glgne*(ratelgne+dsqrt(ratelgne/kklgn)*filtksielgn)
gffi=gbacki*(ratebacki+dsqrt(ratebacki/kkee)*filtksiiback)+glgni*(ratelgni+dsqrt(ratelgni/kklgn)*filtksiilgn)


do k=1,ne
if(semact(k).gt.eps9) semact(k)=semact(k)*umact
end do  
do k=1,ni
if(simact(k).gt.eps9) simact(k)=simact(k)*umact
end do  

heenmda=exeenmda2*heenmda
heeampa=exeeampa2*heeampa

do k=1,ne
if(switche(k).eq.1) then
do j=1,ncee(k)
nn=mee(j,k)
heenmda(nn)=heenmda(nn)+feenmda2(k)	
heeampa(nn)=heeampa(nn)+feeampa2(k)	
end do
end if
end do
hei=exei2*hei
do k=1,ni
if(switchi(k).eq.1) then
do j=1,ncei(k)
nn=mei(j,k)
hei(nn)=hei(nn)+fei2(k)
end do
end if
end do
inpee=(geenmda*heenmda+geeampa*heeampa)*(alphacondcur*(xe(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne)) +gffe*(alphacondcur*(xe(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne))
inpei=gei*hei*(alphacondcur*(xe(1,:)-vsyni)+(1.-alphacondcur)*(vl-vsyni))
inpeeclamp=(geenmda*heenmda+geeampa*heeampa+gffe)*(alphacondcur*(vclamp-vsyne)+(1.-alphacondcur)*(vl-vsyne))
inpeiclamp=gei*hei*(alphacondcur*(vclamp-vsyni)+(1.-alphacondcur)*(vl-vsyni))
he=inpee+inpei
!he=(geenmda*heenmda+geeampa*heeampa+gffe)*(alphacondcur*(xe(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne))+gei*hei*(alphacondcur*(xe(1,:)-vsyni)+(1.-alphacondcur)*(vl-vsyni))
if(it.gt.it1) then 
avhee=avhee+inpeeclamp
avhee2=avhee2+inpeeclamp**2
avhei=avhei+inpeiclamp
avhei2=avhei2+inpeiclamp**2
avgee=avgee+alphacondcur*(geenmda*heenmda+geeampa*heeampa+gffe)
avgei=avgei+alphacondcur*gei*hei
avhetot=avhetot+inpeeclamp+inpeiclamp
avhetot2=avhetot2+(inpeeclamp+inpeiclamp)**2
end if
hienmda=exienmda2*hienmda
hieampa=exieampa2*hieampa

do k=1,ne
if(switche(k).eq.1) then
do j=1,ncie(k)
nn=mie(j,k)
hienmda(nn)=hienmda(nn)+fienmda2(k)
hieampa(nn)=hieampa(nn)+fieampa2(k)
end do
end if
end do	     

hii=exii2*hii
do k=1,ni
if(switchi(k).eq.1) then
do j=1,ncii(k)
nn=mii(j,k)
hii(nn)=hii(nn)+fii2(k) 
end do
end if
end do
inpie=(gienmda*hienmda+gieampa*hieampa)*(alphacondcur*(xi(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne)) +gffi*(alphacondcur*(xi(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne))
inpii=gii*hii*(alphacondcur*(xi(1,:)-vsyni)+ (1.-alphacondcur)*(vl-vsyni))
inpieclamp=(gienmda*hienmda+gieampa*hieampa+gffi)*(alphacondcur*(vclamp-vsyne)+(1.-alphacondcur)*(vl-vsyne))
inpiiclamp=gii*hii*(alphacondcur*(vclamp-vsyni)+ (1.-alphacondcur)*(vl-vsyni))
!hi=(gienmda*hienmda+gieampa*hieampa+gffi)*(alphacondcur*(xi(1,:)-vsyne)+(1.-alphacondcur)*(vl-vsyne))+gii*hii*(alphacondcur*(xi(1,:)-vsyni)+ (1.-alphacondcur)*(vl-vsyni))
hi=inpie+inpii
if(it.gt.it1) then 
avhie=avhie+inpieclamp
avhie2=avhie2+inpieclamp**2
avhii=avhii+inpiiclamp
avhii2=avhii2+inpiiclamp**2
avgie=avgie+alphacondcur*(gienmda*hienmda+gieampa*hieampa+gffi)
avgii=avgii+alphacondcur*gii*hii
avhitot=avhitot+(inpiiclamp+inpieclamp)
avhitot2=avhitot2+(inpiiclamp+inpieclamp)**2
end if
he=-he+ies
hi=-hi+iis
call rk4(xe,ne,dtf,xoute,he,gle,gnae,gke,phie,gae,gnape,gze)


do k=1,ne
switche(k)=0
if(xoute(1,k).gt.thresh.and.xe(1,k).lt.thresh) then
if(it.gt.40000) then 
write(27,*) k, sngl(it*dtf)
totspike=totspike+1
end if

switche(k)=1
tspike(k)=it*dtf
dte(k)=tspike(k)-tspikeold(k)
feenmda2(k)=1
feeampa2(k)=1
fieampa2(k)=1
fienmda2(k)=1
semact(k)=semact(k)+epsact
if(it.gt.it1) then
nspike(k)=nspike(k)+1
fires(k)=fires(k)+dsin(omt)
firec(k)=firec(k)+dcos(omt)
sumspe(k)=sumspe(k)+dte(k)
sumsp2e(k)=sumsp2e(k)+dte(k)**2
sumecv2(k)=sumecv2(k)+ 2*abs(dte(k)-dteold(k))/(dte(k)+dteold(k))
endif
dteold(k)=dte(k)
tspikeold(k)=tspike(k)
end if
end do

call  rk4(xi,ni,dtf,xouti,hi,gli,gnai,gki,phii,gai,gnapi,gzi)

do k=1,ni
switchi(k)=0
	if(xouti(1,k).gt.thresh.and.xi(1,k).lt.thresh) then
if(it.gt.40000) then 
write(28,*) k, sngl(it*dtf)
totspiki=totspiki+1
end if
		switchi(k)=1
   		fei2(k)=1
   		fii2(k)=1
 		simact(k)=simact(k)+epsact
                        tspiki(k)=it*dtf
                        dti(k)=tspiki(k)-tspikiold(k)
			if(it.gt.it1) then
                        	nspiki(k)=nspiki(k)+1
	firis(k)=firis(k)+dsin(omt)
	firic(k)=firic(k)+dcos(omt)
                        	sumspi(k)=sumspi(k)+dti(k)
                        	sumsp2i(k)=sumsp2i(k)+dti(k)**2
	sumicv2(k)=sumicv2(k)+ 2*abs(dti(k)-dtiold(k))/(dti(k)+dtiold(k))
			endif
                        tspikiold(k)=tspiki(k)
			dtiold(k)=dti(k)
		end if
	end do

	xe=xoute
	xi=xouti

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
	  	ue(iu)=ue(iu)+1000*semact((iu-1)*(ne/nu)+i)
  	c1e(iu)=c1e(iu)+1000*semact((iu-1)*(ne/nu)+i)*cose(i)
  	s1e(iu)=s1e(iu)+1000*semact((iu-1)*(ne/nu)+i)*sine(i)
  	c2e(iu)=c2e(iu)+1000*semact((iu-1)*(ne/nu)+i)*cose2(i)
  	s2e(iu)=s2e(iu)+1000*semact((iu-1)*(ne/nu)+i)*sine2(i)
	  end do	
	  do i=1,ni/nu
	  	ui(iu)=ui(iu)+1000*simact((iu-1)*(ni/nu)+i)
  	c1i(iu)=c1i(iu)+1000*simact((iu-1)*(ni/nu)+i)*cosi(i)
  	s1i(iu)=s1i(iu)+1000*simact((iu-1)*(ni/nu)+i)*sini(i)
  	c2i(iu)=c2i(iu)+1000*simact((iu-1)*(ni/nu)+i)*cosi2(i)
  	s2i(iu)=s2i(iu)+1000*simact((iu-1)*(ni/nu)+i)*sini2(i)
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
 write(3,*) sngl(it*dtf),sngl(ue(iu)*nu),sngl(ui(iu)*nu)
! write(36,*) sngl(it*dtf),sngl(xe(1,n1e)),sngl(xe(1,n2e)),sngl(xe(1,n3e))
!  write(38,*) sngl(it*dtf),-sngl(inpee(n1e)),-sngl(inpei(n1e)),-sngl(he(n1e))
!  write(37,*) sngl(it*dtf),sngl(xi(1,n1i)),sngl(xi(1,n2i)),sngl(xi(1,n3i))
! write(39,*) sngl(it*dtf),-sngl(inpie(n1e)),-sngl(inpii(n1e)),-sngl(hi(n1e))

if(it.gt.it1) then
do i=1,nhyper
avervehyper(i)=xe(1,i)/float(itmax-it1)+avervehyper(i)
avervihyper(i)=xi(1,i)/float(itmax-it1)+avervihyper(i)
avervehypermodulc(i)=xe(1,i)/float(itmax-it1)*dcos(omt)+avervehypermodulc(i)
avervehypermoduls(i)=xe(1,i)/float(itmax-it1)*dsin(omt)+avervehypermoduls(i)
avervihypermodulc(i)=xi(1,i)/float(itmax-it1)*dcos(omt)+avervihypermoduls(i)
avervihypermoduls(i)=xi(1,i)/float(itmax-it1)*dsin(omt)+avervihypermoduls(i)
avervehyper2(i)=xe(1,i)**2/float(itmax-it1)+avervehyper2(i)
avervihyper2(i)=xi(1,i)**2/float(itmax-it1)+avervihyper2(i)
end do
end if

if(it.gt.it1) then
veall=0
viall=0
do i=1,ne
veall=veall+xe(1,i)
end do
veall=veall/float(ne)
averveall=averveall+veall
averveall2=averveall2+veall**2
do i=1,ne
averve(i)=xe(1,i)+averve(i)
averve2(i)=xe(1,i)**2+averve2(i)
end do
do i=1,ni
viall=viall+xi(1,i)
end do
viall=viall/float(ni)
averviall=averviall+viall
averviall2=averviall2+viall**2
do i=1,ni
avervi(i)=xi(1,i)+avervi(i)
avervi2(i)=xi(1,i)**2+avervi2(i)
end do
averveallt=averveall/float(it-it1)
averviallt=averviall/float(it-it1)
averveall2t=averveall2/float(it-it1)
averviall2t=averviall2/float(it-it1)
sigveall=averveall2t-(averveallt)**2
sigviall=averviall2t-(averviallt)**2
sigve=averve2/float(it-it1)-(averve/float(it-it1))**2
sigvi=avervi2/float(it-it1)-(avervi/float(it-it1))**2
allsigve=0
allsigvi=0
do i=1,ne
allsigve=sigve(i)+allsigve
end do
do i=1,ni
allsigvi=sigvi(i)+allsigvi
end do
allsigve=allsigve/float(ne)
allsigvi=allsigvi/float(ni)
khie=(sigveall/allsigve)
khii=(sigviall/allsigvi)
if(it.gt.it1+1) write(96,*) sngl(it*dtf),sngl(sqrt(khie)),sngl(sqrt(khii))
end if



        end do

	end do !it

	close(3)

        do i=1,ne
                fire(i)=nspike(i)/float(itmax-it1)/dtf
                fires(i)=fires(i)/float(itmax-it1)/dtf
                firec(i)=firec(i)/float(itmax-it1)/dtf
                if(nspike(i).gt.2) then
                        sumspe(i)=sumspe(i)/nspike(i)
                        sumsp2e(i)=sumsp2e(i)/nspike(i)
                        cve(i)=sqrt(sumsp2e(i)-sumspe(i)**2)/sumspe(i)
			cv2e(i)=sumecv2(i)/nspike(i)
                endif
        write(71,*) i,sngl(1000.*firec(i)),sngl(1000.*fires(i))
                if((nspike(i).gt.2)) then
        write(31,*) i,sngl(1000*fire(i)),sngl(cve(i))
        write(30,*) log(sngl(1000*fire(i))+0.000001)/log(10.),sngl(cve(i))
  		else
	  write(31,*) i,sngl(1000.d0*fire(i)),0
		endif
        end do

        do i=1,ni
                firi(i)=nspiki(i)/float(itmax-it1)/dtf
                firis(i)=firis(i)/float(itmax-it1)/dtf
                firic(i)=firic(i)/float(itmax-it1)/dtf
                if(nspiki(i).gt.2) then
                        sumspi(i)=sumspi(i)/nspiki(i)
                        sumsp2i(i)=sumsp2i(i)/nspiki(i)
                        cvi(i)=sqrt(sumsp2i(i)-sumspi(i)**2)/sumspi(i)
			cv2i(i)=sumicv2(i)/nspiki(i)
                endif
        write(72,*) i,sngl(1000.*firic(i)),sngl(1000.*firis(i))
                if((nspiki(i).gt.2)) then
         write(32,*) i,sngl(1000.d0*firi(i)),sngl(cvi(i))
       write(30,*) log(sngl(1000*firi(i))+0.000001)/log(10.),sngl(cvi(i))
                       	else
	  write(32,*) i,sngl(1000.d0*firi(i)),0
		endif
        end do	
	write(29,*) totspike/dfloat(itmax-40000)/dtf*0.1d0, totspiki/dfloat(itmax-40000)/dtf*0.1
	avhee=-avhee/float(itmax-it1)
	avhei=-avhei/float(itmax-it1)
	avhii=-avhii/float(itmax-it1)
	avhie=-avhie/float(itmax-it1)
	avhetot=-avhetot/float(itmax-it1)
	avhitot=-avhitot/float(itmax-it1)
	avhee2=avhee2/float(itmax-it1)
	avhei2=avhei2/float(itmax-it1)
	avhii2=avhii2/float(itmax-it1)
	avhie2=avhie2/float(itmax-it1)
	avhetot2=avhetot2/float(itmax-it1)
	avhitot2=avhitot2/float(itmax-it1)
	avgee=avgee/float(itmax-it1)
	avgei=avgei/float(itmax-it1)
	avgii=avgii/float(itmax-it1)
	avgie=avgie/float(itmax-it1)
	do i=1,ne 
	write(56,*) i,sngl(avhee(i)),sngl(avhei(i)),sngl(avhetot(i)) 
	write(66,*) i,sqrt(sngl(avhee2(i)-avhee(i)**2)),sqrt(sngl(avhei2(i)-avhei(i)**2)),sqrt(sngl(avhetot2(i)-avhetot(i)**2)) 
	write(58,*) i,sngl(avgee(i)),sngl(avgei(i)),sngl(avgee(i)+avgei(i)) 
	end do
	do i=1,ni 
	write(57,*) i,sngl(avhie(i)),sngl(avhii(i)),sngl(avhitot(i)) 
	write(67,*) i,sqrt(sngl(avhie2(i)-avhie(i)**2)),sqrt(sngl(avhii2(i)-avhii(i)**2)),sqrt(sngl(avhitot2(i)-avhitot(i)**2)) 
	write(59,*) i,sngl(avgie(i)),sngl(avgii(i)),sngl(avgie(i)+avgii(i)) 
	end do
do i=1,nhyper
write(60,*) i,sngl(avervehyper(i)),sngl(avervehyper2(i)-avervehyper(i)**2)
write(80,*) i,sngl(avervehypermodulc(i)),sngl(avervehypermoduls(i))
write(81,*) i,sngl(avervihypermodulc(i)),sngl(avervihypermoduls(i))
write(61,*) i,sngl(avervihyper(i)),sngl(avervihyper2(i)-avervihyper(i)**2)
end do
	write(12,*) 'ireal=',ireal
	call flush(12)
		end do !ireal


	end  !---- FIN DEL PROGRAMA ----

subroutine crmat(sem,n1,n2,sqn1,sqn2,sigma,k12,g0,nm,typcon,ncmjj,mee,mei,mie,mii)
!call crmat(semarch,ne,ni,sqne,sqni,sigmaei,kkei,0.d0,nce,'ei',ncei,mee,mei,mie,mii)
	
	implicit none
	integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
	real*8 k12
	parameter(n=19600)
	parameter (kk=300)
        parameter(nc=12000,nce=12000,nci=12000)
	parameter(ne=19600,ni=19600)
	real*8 gg(-n:n) 	
	integer ng		
	integer i,j,il,ncm
	real*8  pp,rnd,ncmm,sumg
	real*8 sigma,g0,pi,dpi,p12,alpha
	real*8 prob(0:ne/2)
 	integer mee(nce,ne),mei(nce,ni),mii(nci,ni),mie(nci,ne)
        real*8 ggubfs,sem,xi1,xi2,yi1,yi2
	real*8 x1,x2,y1,y2,p12x,p12y,snorm1,snorm2,snorm
	integer ncmjj(n),ij,sqn1,sqn2,ix1,ix2,iy1,iy2,k
	character*2 typcon

 	dpi=2.d0*acos(-1.d0)
 	pi=acos(-1.d0)
!       -----------------------------------------------------------------------

	if(typcon.eq.'ee') then
        do i=1,nce
        do j=1,ne
        mee(i,j)=0
        end do
        end do
	end if

	if(typcon.eq.'ei') then
        do i=1,nce
        do j=1,ni
        mei(i,j)=0
        end do
        end do
	end if

	if(typcon.eq.'ie') then
        do i=1,nci
        do j=1,ne
        mie(i,j)=0
        end do
        end do
        end if

	if(typcon.eq.'ii') then
        do i=1,nci
        do j=1,ni
        mii(i,j)=0
        end do
        end do
        end if

	snorm1=0
	do i=1,sqn1
	x1=float(i)/float(sqn1)
	do k=-40,40
	snorm1=snorm1+dexp(-(x1-k)**2/2./sigma**2)
 	end do
	end do

         snorm2=0
        do j=1,sqn2
        x2=float(j)/float(sqn2)
        do k=-40,40
        snorm2=snorm2+dexp(-(x2-k)**2/2./sigma**2)
        end do
        end do




	 do j=1,n2
                ncm=0
                do i=1,n1
	iy1=1+(i-1)/sqn1
	ix1=i-(iy1-1)*sqn1
	iy2=1+(j-1)/sqn2
	ix2=j-(iy2-1)*sqn2
	x1=float(ix1)/float(sqn1)
	x2=float(ix2)/float(sqn2)
	y1=float(iy1)/float(sqn1)
	y2=float(iy2)/float(sqn2)
!write(98,*) i,ix1,iy1
!	write(99,*) j,ix2,iy2
	p12x=0
	p12y=0
	do k=-40,40
	p12x=p12x+dexp(-(x1-x2-k)**2/2./sigma**2)
	p12y=p12y+dexp(-(y1-y2-k)**2/2./sigma**2)
	end do
	p12=p12x*p12y/snorm1/snorm2*k12*sqn1/sqn2
!write(98,*) p12
        rnd=ggubfs(sem)
                  if (rnd.lt.p12) then
                          ncm=ncm+1
!j presynaptic; n2
        if(typcon.eq.'ee') then
	mee(ncm,j)=i
	end if
        if(typcon.eq.'ei') then
	mei(ncm,j)=i
	end if
        if(typcon.eq.'ie') then
	mie(ncm,j)=i
	end if
        if(typcon.eq.'ii') then
	mii(ncm,j)=i
	end if
                          if (ncm.eq.nm) then
                           write(*,*) 'error en crmat: ncm=nc',i,il
                           stop
                   end if
                   end if
                end do
        ncmjj(j)=ncm
                write(97,*) j,ncm
        end do

        do i=1,n2
        if(typcon.eq.'ee') write(4,*) i,mee(1,i)
        if(typcon.eq.'ei') write(4,*) i,mei(1,i)
        if(typcon.eq.'ie') write(4,*) i,mie(1,i)
        if(typcon.eq.'ii') write(4,*) i,mii(1,i)
        end do

!	open(4,file='histg.dat')
!	do i=1,nc
!		write(4,*) i,mm(5000,i)
!	end do
!	close(4)

	return
	end


	subroutine crmatcos(mm,sem,g1,n1,n2,k12,nm,theta1,theta2)
	implicit none
	integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
	real*8 k12
	parameter(n=19600)
	parameter (kk=300)
	parameter(nc=12000)
	real*8 gcos(-n2:n2),g1   	  !distrib de prob.de conexion
	integer mm(n,nc),i,j,il,ncm,nn(n)
	real*8  pp,rnd,sumg,dpi,theta1(n),theta2(n),p12
        real*8 sem,ggubfs
               
	dpi=2.d0*acos(-1.d0)

	do i=1,n
	nn(i)=0
	do j=1,nc
                mm(i,j)=0
	end do
	end do

	do i=1,n1
		ncm=0
		do j=1,n2
	p12=(1.d0+g1*dcos(theta1(i)-theta2(j)))*k12/dfloat(n2)
		   rnd=ggubfs(sem)
               	   if (rnd.lt.p12) then
			  ncm=ncm+1
		 	  mm(i,ncm)=j		
			  nn(j)=nn(j)+1
			  if (ncm.eq.nm) then
                 write(*,*) 'error en crmatcos: ncm=nc',i,il
                                stop
                          end if 
	           end if
		end do  
	end do  !i
	do i=1,n1
		write(4,*) i,mm(i,1)
	end do
!	open(4,file='histcos.dat')
!	do i=1,nc
!		write(4,*) i,mm(5000,i)
!	end do
!	close(4)

	return
	end
!=====================================================================
      subroutine  rk4(x,ntot,dt,xout,isyn,gl,gna,gk,phi,ga,gnap,gz)
      integer n,ntot
      PARAMETER (n=19600)
      real*8 x(5,n),dxdt(5,n),xout(5,n)
      real*8 xt(5,n),dxt(5,n),dxm(5,n)
      real*8 isyn(n)
      real*8 hh,h6,dt,t
      real*8 gl,gna(n),gk(n),ga,gnap,phi
	real*8 gz(n)
      HH=dt*0.5
      H6=dt/6.
      call devt(x,dxdt,isyn,ntot,gl,gna,gk,phi,ga,gnap,gz)
      call copyscale(x,dxdt,xt,ntot,hh)
!	xt=x+hh*dxdt
      call devt(xt,dxt,isyn,ntot,gl,gna,gk,phi,ga,gnap,gz)
      call copyscale(x,dxt,xt,ntot,hh)
!	xt=x+hh*dxt
      call devt(xt,dxm,isyn,ntot,gl,gna,gk,phi,ga,gnap,gz)
      call copyscale(x,dxm,xt,ntot,dt)
!	xt=x+dt*dxm
      hh=1.
      call copyscale(dxt,dxm,dxm,ntot,hh)
!	dxm=dxt+dxm
      call devt(xt,dxt,isyn,ntot,gl,gna,gk,phi,ga,gnap,gz)
      call rk4last(x,dxdt,dxt,dxm,xout,h6,ntot)
!       xout=x+h6*(dxdt+dxt+2.*dxm)
        return
        end
!===================================================================
      subroutine devt(x,f,isyn,ntot,gl,gna,gk,phi,ga,gnap,gz)
      integer n,neqmax
      parameter (n=19600)
      real*8 x(5,n),f(5,n)
      real*8 v(n),h(n),z(n)
      real*8 b(n),xn(n)
      real*8 minfv(n),ninfv(n),zinfv(n)
      real*8 ainfv(n),binfv(n)
      real*8 sainfv(n),am(n),bm(n),ah(n),bh(n)
      real*8 an(n),bn(n),taubv(n)
      real*8 taunv(n),tauhv(n),hinfv(n)
      real*8 vna,vl,vk,va,taua
      real*8 gl,gna(n),gk(n),ga,gnap
      real*8 tauz,phi
      real*8 isyn(n),t,gz(n)
      integer in,ntot
      common/cb/ tauz,taua,vna,vk,vl,va

        do in=1,ntot
        v(in)=x(1,in)
        h(in)=x(2,in)
        xn(in)=x(3,in)
        b(in)=x(4,in)
        z(in)=x(5,in)
        end do
        do in=1,ntot
        if(abs(v(in)+30.).lt.0.0001) then
        am(in)=1.
          else
        am(in)=0.1*(v(in)+30.)/(1.-exp(-0.1*(v(in)+30.)))
           end if
        bm(in)=4.*exp(-0.055555*(v(in)+55.))

        ah(in)=0.07*exp(-0.05*(v(in)+44.))
        bh(in)=1./(1.+exp(-0.1*(v(in)+14.)))

         if(abs(v(in)+34).lt.0.0001) then
            an(in)=0.1
          else
        an(in)=0.01*(v(in)+34)/(1.-exp(-0.1*(v(in)+34)))
             end if
        bn(in)=0.125*exp(-0.0125*(v(in)+44))

        minfv(in)= am(in)/(am(in)+bm(in))
        hinfv(in)= ah(in)/(ah(in)+bh(in))
        tauhv(in)=1./(ah(in)+bh(in))
        ninfv(in)= an(in)/(an(in)+bn(in))
        taunv(in)=1./(an(in)+bn(in))
        ainfv(in)=1./(1.+exp(-(v(in)+50.)/20.))
        binfv(in)=1./(1.+exp((v(in)+80.)/6.))
        taubv(in)=taua
        zinfv(in)=1./(1.+exp(-0.7*(v(in)+30.)))
        sainfv(in)=1./(1.+exp(-0.3*(v(in)+50)))
         end do

           do  in=1,ntot
     f(1,in)=-gna(in)*minfv(in)**3*h(in)*(v(in)-vna)-gl*(v(in)-vl)-gk(in)*xn(in)**4*(v(in)-vk)-ga*ainfv(in)**3*b(in)*(v(in)-va)-gnap*sainfv(in)*(v(in)-vna)-gz(in)*z(in)*(v(in)-vk)+isyn(in)
        f(2,in)= phi*(hinfv(in)-h(in))/tauhv(in)
        f(3,in)= phi*(ninfv(in)-xn(in))/taunv(in)
        f(4,in)= (binfv(in)-b(in))/taubv(in)
        f(5,in)= (zinfv(in)-z(in))/tauz
        end do
        return
        end
!==================================================================
        function ggubfs(seed)
        real*8 ggubfs
        real*8 seed
        real*8 d2p31m,d2p31
        data d2p31m /2147483647.d0/
        data d2p31  /2147483711.d0/

        seed = dmod(16807.d0*seed,d2p31m)
        ggubfs=seed/d2p31
        return
        end
!==================================================================
        subroutine copyscale(x,y,z,ntot,hh)
        integer n,i
	parameter(n=19600)
        real*8 x(5,n),y(5,n),z(5,n),hh
	z=x+hh*y
        return
        end
!==================================================================
        subroutine rk4last(x,dxdt,dxt,dxm,xout,h6,ntot)
        integer n,i,ntot
	parameter(n=19600)
        real*8 x(5,n),dxdt(5,n),dxt(5,n),xout(5,n),h6,dxm(5,n)
        xout=x+h6*(dxdt+dxt+2.*dxm)
        return
        end
