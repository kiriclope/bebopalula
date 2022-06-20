c       SPARSE NETWORK MODEL 
c       THREE POPULATIONS: se, si, sp
c       CONECTIVITY: K COUPLINGS:
c       FIRST ORDER IMPLICIT INTEGRATION
c       INTERACTION: 1 EXPONENTIAL
c       
C       E->I' Facilitating   I'->I depressing  E-E depressing
	PROGRAM INTFIRE
	implicit none
	integer n,nc,kk,nu,nce,nci
	integer ne,ni
	real*8 kkee,kkei,kkii,kkie
	real*8 kkee0,kkei0,kkii0,kkie0
	parameter(n=64000)	!poner multiplo de nu
	parameter(nu=1) 
	parameter (kk=300) 
        parameter(nc=1500,nce=1500,nci=500)
        parameter(ne=64000,ni=16000)
	integer i,j,k,iu,it,itmax,iseed1,iseed2

	integer it1,it2,it3,it4,itrans
	integer i1,i2,i3
	real*8 se(0:n),si(0:n)	!estados de las poblaciones
	real*8 ue(nu),ui(nu)	!estados de las poblaciones
	real*8 semact(0:n),simact(0:n)
	real*8 seed(0:n),seem(0:n)
	real*8 eps,umeps,taufac,taurec,bigu !valor medio temporal
	real*8 taufacee,taurecee,biguee
	real*8 gee,geeampa,gei,gii,gie,ie,ii !pesos de conexion
	real*8 g1ee,g1ei,g1ii,g1ie !pesos de conexion
	real*8 sigmaee,sigmaei,sigmaii,sigmaie !pesos de conexion
	real*8 gpec
	real*8 tauact,epsact,umact
	real*8 gee0,gampa0,gei0,gii0,gie0,gieampa
	real*8 gpi0
	real*8 gg(-n:n),sumg,sigma !distr.de prob.de conexion
	real*8 ggee(-n:n),sumgee !distr.de prob.de conexion
        real*8 rampanmda,rampanmdaie
	real*8 rescale,rescalexci,rescalinhib
	real*8 xx

	integer ng,ngee		!ancho de gg (en sitios)
        integer ngei,ngie,ngii
	integer mee(ne,nce),mei(ne,nci),mii(ni,nci) !matrices de conexion
	integer mie(ni,nce)
	real*8 f,x
	real*8 he(n),hi(n)	!campos internos
	real*8 rnd,ran3,pi,dpi,iesqk,gamma,iisqk,iisqk0
	real*8 iesqk0,iesqk1,iesqk2
	real*8 ampl,ampl0,ies,iis,theta
        real*8 ree(n),rie(n),rei(n),rii(n)
        real*8 hee(n),hie(n),hei(n),hii(n)
	real*8 heeampa(n),hieampa(n)
        integer nspike(n),nspiki(n)
	integer deltit,idir,anneal
	integer ireal,nreal
	character*80 nombfile
        character*5 type_ee,type_ei,type_ie,type_ii
	character*15 type_anneal

        real*8 tee2,tei2,tii2,tie2,teeampa2
	real*8 s1,s2,s3,dtf,dt2
	real*8 exee2,exei2,exii2,exie2
	real*8 expi2
	real*8 exeeampa2
	real*8 fee2(0:n),fei2(0:n),fii2(0:n)
	real*8 fie2(0:n)
	real*8 feeampa2(0:n),fieampa2(0:n)
	real*8 dt,eps9
	real*8 ser,set,sir,sit,spr,spt,dtf2
	real*8 ke1(0:n),ki1(0:n),kp1(0:n)
	real*8 ke2(0:n),ki2(0:n),kp2(0:n)
	real*8 sebis(0:n),sibis(0:n),spbis(0:n),snew,uu
	real*8 taue,taui

        real*8 tspike(n),tspikeold(n)
        real*8 tspiki(n),tspikiold(n)
        real*8 tspikp(n),tspikpold(n)
        real*8 sumspe(n),sumsp2e(n),sumspi(n)
	real*8 sumsp2i(n)
	real*8 sumecv2(n),sumicv2(n),sumpcv2(n)
	real*8 dte(n),dteold(n),dti(n),dtiold(n),dtp(n),dtpold(n)
        real*8 fire(n),firi(n),firp(n)
        real*8 cve(n),cvi(n),cvp(n)
        real*8 cv2e(n),cv2i(n),cv2p(n)
	real*8 coeff,coef,coef1,coef2,coefk
	real*8 vpspeenmda,vpspeeampa,vpspei,vpspie,vpspii
	integer je,ji,jp

	common /distrib/gg,ng
	common/connect/mee,mei,mie,mii

	f(x)=-x

	eps9=1e-09
	uu=0.0001
 	coeff=1.
	nreal=1

c       
	open(2,file='depeir.ini')
	read(2,*) iseed1
	read(2,*) iseed2
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
c       read(2,*) ne
c       read(2,*) ni
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
	read(2,*) iesqk1
 	read(2,*) iesqk2	
	read(2,*) theta
	read(2,*) ampl0
	read(2,*) anneal
	read(2,*) type_anneal
	read(2,*) nombfile
	close(2)

        geeampa=gee*rampanmda/(1.+rampanmda)
        gee=gee/(1.+rampanmda)
        gieampa=gie*rampanmdaie/(1+rampanmdaie)
        gie=gie/(1+rampanmdaie)
c       
 	open(3,file='ux-depr-'//nombfile)
	open(4,file='conect-depr-'//nombfile)
        open(26,file='spike-depr-'//nombfile)
        open(27,file='spiki-depr-'//nombfile)
        open(28,file='spikp-depr-'//nombfile)
        open(31,file='fire-depr-'//nombfile)
        open(32,file='firi-depr-'//nombfile)
        open(33,file='firp-depr-'//nombfile)
        open(44,file='itot-depr-'//nombfile)
        open(45,file='traces-depr-'//nombfile)
c       
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
        iesqk1=iesqk1*sqrt(1.d0*kkee)*coef
        iesqk2=iesqk2*sqrt(1.d0*kkee)*coef
c       
	exee2=exp(-dtf/tee2)
	exeeampa2=exp(-dtf/teeampa2)
	exei2=exp(-dtf/tei2)
	exii2=exp(-dtf/tii2)
	exie2=exp(-dtf/tie2)

	umeps=exp(-dtf/taufac)
	umact=exp(-dtf/tauact)

	open(12,file='pr-depr-'//nombfile)
	write(12,*) 'N,K=',n,kk
	write(12,*) 'iseed1=',iseed1
	write(12,*) 'iseed2=',iseed2
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
  	write(12,*) 'it3,it4=',it3,it4
        write(12,*) 'iesqk1=',iesqk1/sqrt(1.*kkee)
        write(12,*) 'iisqk1=',iisqk/sqrt(1.*kkii)
	write(12,*) 'theta=',theta

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
	write(12,*) 'nombfile=',nombfile

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
c----------------------------------------------------------
c	calculation PSP's
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
c       --------- Matrices de conexion -------
c	call crmatcos(mee,iseed1,g1ee,ne,ne,kkee,nce)
c	call crmatcos(mei,iseed1,g1ei,ne,ni,kkei,nci)
c	call crmatcos(mie,iseed1,g1ie,ni,ne,kkie,nce)
c	call crmatcos(mii,iseed1,g1ii,ni,ni,kkii,nci)
 	call crmat(iseed1,ne,ne,sigmaee,kkee,0.d0,nce,'ee')
	call crmat(iseed1,ne,ni,sigmaei,kkei,0.d0,nci,'ei')
	call crmat(iseed1,ni,ne,sigmaie,kkie,0.d0,nce,'ie')
	call crmat(iseed1,ni,ni,sigmaii,kkii,0.d0,nci,'ii')

	do ireal=1,nreal
c	---------- Condicion Inicial-----------

	   se(0)=0.		!NOTA: el elemento 0 de todas las poblaciones es 
	   si(0)=0.		!      siempre 0 y no se toca. Est�para cuando aparece 
	   semact(0)=0.
	   simact(0)=0.
	   seem(0)=biguee
	   seed(0)=1
	   do i=1,n
	      fee2(i)=0
	      fie2(i)=0
	      feeampa2(i)=0
	      fei2(i)=0
	      fii2(i)=0

	      he(i)=0 
	      hi(i)=0 

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

	   do i=1,n
	      rnd=ran3(iseed2)
	      se(i)=ser+rnd*(set-ser)
	      rnd=ran3(iseed2)
	      si(i)=sir+rnd*(sit-sir)
	      semact(i)=0
	      simact(i)=0
	      seed(i)=1.
	      seem(i)=biguee
	   end do

c	----------- DINAMICA ----------

c	open(3,file='ux-'//nombfile)
c	open(26,file='spike-'//nombfile)
c	open(27,file='spiki-'//nombfile)
c	open(28,file='spikp-'//nombfile)

	   idir=1
	   ampl=0

	   do it=1,itmax

	      if(anneal.eq.0) then

		 if(it.gt.it1 .and. it.lt.it2) then
		    iesqk=iesqk1
		    ampl=ampl0
		 else if(it.gt.it3 .and. it.lt.it4) then
		    iesqk=iesqk2
		    iisqk=4.*iesqk/sqrt(kkee)*sqrt(kkii)
		    ampl=0
		 else
		    iesqk=iesqk0
		    iisqk=iisqk0
		    ampl=0
		 endif

	      else
		 if(it.ge.itmax/2) idir=-1

		 if(type_anneal.eq.'2dir') 
		    iesqk=iesqk+idir*(iesqk1-iesqk0)/(itmax/2)
		    
		 if(type_anneal.eq.'2direi') then
		    iesqk=iesqk+idir*(iesqk1-iesqk0)/(itmax/2)
		    iisqk=iesqk/sqrt(kkee)*sqrt(kkii)    
		 end if
		 
		 if(type_anneal.eq.'1dir') iesqk=iesqk+(iesqk1-iesqk0)/itmax
c	if(type_anneal.eq.'1dirgpe') gpec=gpe/float(itmax)*float(it)
		 

	      endif

	      do k=1,ne
		 if(fee2(k).gt.eps9) fee2(k)=fee2(k)*exee2  
		 if(feeampa2(k).gt.eps9) feeampa2(k)=feeampa2(k)*exeeampa2  
		 if(fieampa2(k).gt.eps9) fieampa2(k)=fieampa2(k)*exeeampa2  
		 if(fie2(k).gt.eps9) fie2(k)=fie2(k)*exie2 
		 if(semact(k).gt.eps9) semact(k)=semact(k)*umact
	      end do  
	      do k=1,ni
		 if(fei2(k).gt.eps9) fei2(k)=fei2(k)*exei2 
		 if(fii2(k).gt.eps9) fii2(k)=fii2(k)*exii2  
		 if(simact(k).gt.eps9) simact(k)=simact(k)*umact
	      end do  


	      do k=1,ne
		 he(k)=0.
		 ree(k)=0
		 rei(k)=0
		 hee(k)=0
		 heeampa(k)=0
		 hei(k)=0
	      end do

	      do j=1,nce
		 do k=1,ne
		    hee(k)=hee(k)+fee2(mee(k,j))	
		 end do
	      end do
	      if(rampanmda.eq.0) go to 190
	      do j=1,nce
		 do k=1,ne
		    heeampa(k)=heeampa(k)+feeampa2(mee(k,j))	
		 end do
	      end do

 190	      continue

	      do j=1,nci
		 do k=1,ne
		    hei(k)=hei(k)+fei2(mei(k,j))
		 end do
	      end do
	      if(rampanmda.eq.0) then
		 do k=1,ne
		    he(k)=hee(k)*gee+ hei(k)*gei
		 end do
	      else       
		 do k=1,ne
		    he(k)=hee(k)*gee+ heeampa(k)*geeampa+ hei(k)*gei
		 end do
	      end if


	      do k=1,ni
		 hi(k)=0.
		 rie(k)=0
		 rii(k)=0
		 hie(k)=0
		 hieampa(k)=0
		 hii(k)=0
	      end do

	      do j=1,nce
		 do k=1,ni
		    hie(k)=hie(k)+fie2(mie(k,j))
		 end do
	      end do	     

	      do j=1,nce
		 do k=1,ni
		    hieampa(k)=hieampa(k)+fieampa2(mie(k,j))
		 end do
	      end do

	      do j=1,nci
		 do k=1,ni
		    hii(k)=hii(k)+fii2(mii(k,j)) 
		 end do
	      end do


	      do k=1,ni
		 hi(k)=hie(k)*gie+hii(k)*gii+hieampa(k)*gieampa
	      end do



	      do k=1,ne
c       ies=iesqk*(1-ampl*cos(float(k)*dpi/float(ne)-theta*dpi/360.))
		 ies=iesqk*(1+ampl-ampl*cos(float(k)*dpi/float(ne)-theta*dpi/360.))
		 sebis(k)=se(k)/(1+dtf/taue)+dtf*(he(k)+ies)/(taue+dtf)
	      end do
	      iis=iisqk
	      do k=1,ni
		 sibis(k)=si(k)/(1+dtf/taui)+dtf*(hi(k)+iis)/(taui+dtf)
	      end do


	      do k=1,ne
		 if(sebis(k).gt.set) then
		    dt2=dtf*(sebis(k)-set)/(sebis(k)-se(k))
		    snew=ser
	1		 +(sebis(k)-set)*(1+dtf/taue*(f(ser)-f(se(k)))/(sebis(k)-se(k)))
		    tspike(k)=it*dtf
		    dte(k)=tspike(k)-tspikeold(k)
		    seem(k)=seem(k)*(1.-biguee)*dexp(-dte(k)/taufacee)+biguee
		    xx=dexp(-dte(k)/taurecee)
		    seed(k)=seed(k)*(1.-seem(k))*xx+1.-xx
		    fee2(k)=fee2(k)+exp(-dt2/tee2)*seem(k)*seed(k)
		    feeampa2(k)=feeampa2(k)+exp(-dt2/teeampa2)*seem(k)*seed(k)
		    fieampa2(k)=fieampa2(k)+exp(-dt2/teeampa2)


		    fie2(k)=fie2(k)+exp(-dt2/tie2)
c       fee2(k)=fee2(k)+(1.-dt2/tee2)
c       fie2(k)=fie2(k)+(1.-dt2/tie2)
c	if(mod(k,10).eq.0) write(26,*) (it+(ireal-1)*itmax)*dtf*0.01,
		    write(26,*) (it+(ireal-1)*itmax)*dtf*0.01,
     *                   k,sebis(k)-se(k)
		    call flush(26)
		    se(k)=snew
		    semact(k)=semact(k)+epsact*exp(-dt2/tauact)
		    if(it.gt.itrans .and. it.lt.it3) then
		       nspike(k)=nspike(k)+1
		       sumspe(k)=sumspe(k)+dte(k)
		       sumsp2e(k)=sumsp2e(k)+dte(k)**2
		       sumecv2(k)=sumecv2(k)
	1		    + 2*abs(dte(k)-dteold(k))/(dte(k)+dteold(k))
		    endif
		    dteold(k)=dte(k)
		    tspikeold(k)=tspike(k)
		 else 
		    se(k)=sebis(k)
		 end if
	      end do


	      do k=1,ni
		 if(sibis(k).gt.sit) then
		    dt2=dtf*(sibis(k)-sit)/(sibis(k)-si(k))
		    snew=sir
	1		 +(sibis(k)-sit)*(1+dtf/taui*(f(sir)-f(si(k)))/(sibis(k)-si(k)))
		    fei2(k)=fei2(k)+exp(-dt2/tei2)
		    fii2(k)=fii2(k)+exp(-dt2/tii2)
c       fei2(k)=fei2(k)+(1.-dt2/tei2)
c       fii2(k)=fii2(k)+(1.-dt2/tii2)
c       if(mod(k,10).eq.0) write(27,*) (it+(ireal-1)*itmax)*dtf*0.01,
		    write(27,*) (it+(ireal-1)*itmax)*dtf*0.01,
	1		 k,sibis(k)-si(k)
		    call flush(27)
		    si(k)=snew
		    simact(k)=simact(k)+epsact*exp(-dt2/tauact)
c       sim(k)=sim(k)+eps*(1-dt2/taufac)

		    tspiki(k)=it*dtf
		    dti(k)=tspiki(k)-tspikiold(k)
		    if(it.gt.itrans .and. it.lt.it3) then
		       nspiki(k)=nspiki(k)+1
		       sumspi(k)=sumspi(k)+dti(k)
		       sumsp2i(k)=sumsp2i(k)+dti(k)**2
		       sumicv2(k)=sumicv2(k)
	1		    + 2*abs(dti(k)-dtiold(k))/(dti(k)+dtiold(k))
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
		    do i=1,ne/nu
		       ue(iu)=ue(iu)+100*semact((iu-1)*(ne/nu)+i)
		    end do	
		    do i=1,ni/nu
		       ui(iu)=ui(iu)+100*simact((iu-1)*(ni/nu)+i)
		    end do	

		    if(anneal.eq.0) then
		       write(3,*) sngl(it*dtf*0.01d0),
	1		    sngl((ue(iu)*nu)/ne),sngl((ui(iu)*nu)/ni),sngl(iesqk),sngl(iisqk)
		       write(45,*) sngl(it*dtf*0.01d0),sngl(se(1)),sngl(si(1))
		    else
		       if(type_anneal.eq.'1dirgpe') then 
			  write(3,*) sngl(gpec),
	1		       sngl((ue(iu)*nu)/ne),sngl((ui(iu)*nu)/ni)
		       else
			  write(3,*) sngl(iesqk/sqrt(1.*kkee)),sngl(iisqk/sqrt(1.*kkii)),
	1		       sngl((ue(iu)*nu)/ne),sngl((ui(iu)*nu)/ni)
		       endif
		    endif
		 end do
		 call flush(3)
		 call flush(45)

		 je=1
		 ji=1
		 jp=1
		 write(44,*)  sngl(it*dtf*0.01d0),
	1	      sngl(abs(ree(je)+iesqk+rei(je))/(ree(je)+iesqk-rei(je))),
	2	      sngl(abs(rie(ji)+rii(ji)))/(rie(ji)-rii(ji))

		 call flush(44)

	      end if
	      
	   end do		!it

	   close(3)

	   do i=1,ne
	      fire(i)=nspike(i)/float(it3-itrans)/dtf
	      if(nspike(i).gt.2) then
		 sumspe(i)=sumspe(i)/nspike(i)
		 sumsp2e(i)=sumsp2e(i)/nspike(i)
		 cve(i)=sqrt(sumsp2e(i)-sumspe(i)**2)/sumspe(i)
		 cv2e(i)=sumecv2(i)/nspike(i)
	      endif
	      if((nspike(i).gt.2)) then
		 write(31,*) sngl(i*360.d0/ne),sngl(100*fire(i)),
	1	      sngl(cve(i)),sngl(cv2e(i))
	      else
		 write(31,*) sngl(i*360.d0/ne),sngl(100.d0*fire(i)),0,0
	      endif
	   end do


	   do i=1,ni
	      firi(i)=nspiki(i)/float(it3-itrans)/dtf
	      if(nspiki(i).gt.2) then
		 sumspi(i)=sumspi(i)/nspiki(i)
		 sumsp2i(i)=sumsp2i(i)/nspiki(i)
		 cvi(i)=sqrt(sumsp2i(i)-sumspi(i)**2)/sumspi(i)
		 cv2i(i)=sumicv2(i)/nspiki(i)
	      endif
	      if((nspiki(i).gt.2)) then
                 write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),
	1	      sngl(cvi(i)),sngl(cv2i(i))
	      else
		 write(32,*) sngl(i*360.d0/ni),sngl(100.d0*firi(i)),0,0
	      endif
	   end do

	   write(12,*) 'ireal=',ireal
	   call flush(12)
	end do			!ireal


	end			!---- FIN DEL PROGRAMA ----

	subroutine crmat(iseed,n1,n2,sigma,k12,g0,nm,typcon)
	
	implicit none
	integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
	real*8 k12
	parameter(n=64000)
	parameter (kk=300)
        parameter(nc=1500,nce=1500,nci=500)
	parameter(ne=64000,ni=16000)
	real*8 gg(-n:n) 	!distr.de prob.de conexion
	integer ng		!ancho de gg (en sitios)
	integer i,j,il,iseed,ncm
	real*8  pp,rnd,ran3,ncmm,sumg
	real*8 sigma,g0,pi,dpi,p12
	real*8 theta1(n),theta2(n),deltatheta,sigmarad
	integer mee(ne,nce),mei(ne,nci),mii(ni,nci) !matrices de conexion
	integer mie(ni,nce)
	character*2 typcon
	common/connect/mee,mei,mie,mii


 	dpi=2.d0*acos(-1.d0)
 	pi=acos(-1.d0)
c       -----------------------------------------------------------------------
        do i=1,n1
	   theta1(i)=dpi*float(i)/float(n1)
        end do
        do i=1,n2
	   theta2(i)=dpi*float(i)/float(n2)
        end do

	if(typcon.eq.'ee') then
	   do i=1,n1
	      do j=1,nce
		 mee(i,j)=0
	      end do
	   end do
	end if

	if(typcon.eq.'ei') then
	   do i=1,n1
	      do j=1,nci
		 mei(i,j)=0
	      end do
	   end do
	end if

	if(typcon.eq.'ie') then
	   do i=1,n1
	      do j=1,nce
		 mie(i,j)=0
	      end do
	   end do
        end if

	if(typcon.eq.'ii') then
	   do i=1,n1
	      do j=1,nci
		 mii(i,j)=0
	      end do
	   end do
        end if

	sigmarad=sigma/360.*dpi

	sumg=0.
	do i=1,n2/2
	   sumg=sumg+g0+exp(-theta2(i)**2/2/sigmarad**2)
      	end do
	sumg=2.*sumg

c	open(4,file='ncmg.dat')

	do i=1,n1
	   ncm=0
	   do j=1,n2
	      deltatheta=dabs(theta1(i)-theta2(j))
	      if(deltatheta.gt.pi) deltatheta=dpi-deltatheta
	      p12=k12*(g0+dexp(-deltatheta**2/2/sigmarad**2))/sumg
	      rnd=ran3(iseed) 
	      if (rnd.lt.p12) then
		 ncm=ncm+1
		 if(typcon.eq.'ee') mee(i,ncm)=j	
		 if(typcon.eq.'ei') mei(i,ncm)=j	
		 if(typcon.eq.'ie') mie(i,ncm)=j	
		 if(typcon.eq.'ii') mii(i,ncm)=j	
		 if (ncm.eq.nm) then
		    write(*,*) 'error en crmat: ncm=nc',i,il
		    stop
		 end if
	      end if
	   end do  
	end do  

        do i=1,n1
	   if(typcon.eq.'ee') write(4,*) i,mee(i,1)
	   if(typcon.eq.'ei') write(4,*) i,mei(i,1)
	   if(typcon.eq.'ie') write(4,*) i,mie(i,1)
	   if(typcon.eq.'ii') write(4,*) i,mii(i,1)
        end do

c	open(4,file='histg.dat')
c	do i=1,nc
c       write(4,*) i,mm(5000,i)
c	end do
c	close(4)

	return
	end



c	---------  RUTINA QUE CREA LAS MATRICES DE CONEXION CON COSENO --------


	subroutine crmatcos(mm,iseed,g1,n1,n2,k12,nm)
	implicit none
	integer n,nc,kk,n1,n2,nce,nci,nm,ne,ni
	real*8 k12
	parameter(n=64000)
	parameter (kk=300)
	parameter(nc=1500)
	real*8 gcos(-n2:n2),g1	!distrib de prob.de conexion
	integer mm(n,nc),i,j,il,iseed,ncm,nn(n)
	real*8  pp,rnd,ran3,sumg,dpi,theta1(n),theta2(n),p12
c	---- gcos: Distribucion tipo 1+cos(x) para prob. de conexion de ie ----
	dpi=2.d0*acos(-1.d0)
c	-----------------------------------------------------------------------
	do i=1,n1
	   theta1(i)=dpi*float(i)/float(n1)
	end do
	do i=1,n2
	   theta2(i)=dpi*float(i)/float(n2)
	end do

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
	      rnd=ran3(iseed) 
	      if (rnd.lt.p12) then
		 ncm=ncm+1
		 mm(i,ncm)=j		
		 nn(j)=nn(j)+1
		 if (ncm.eq.nm) then
		    write(*,*)
	1		 'error en crmatcos: ncm=nc',i,il
		    stop
		 end if 
	      end if
	   end do  
	end do			!i
	do i=1,n1
	   write(4,*) i,mm(i,1)
	end do
c	open(4,file='histcos.dat')
c	do i=1,nc
c       write(4,*) i,mm(5000,i)
c	end do
c	close(4)

	return
	end
c--------------------------------------------------------------------
	FUNCTION ran3(idum)
	INTEGER idum
	INTEGER MBIG,MSEED,MZ
C       REAL MBIG,MSEED,MZ
	DOUBLE PRECISION ran3,FAC
	PARAMETER (MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.d0/MBIG)
C       PARAMETER (MBIG=4000000.d0,MSEED=1618033.d0,MZ=0.d0,FAC=1.d0/MBIG)
	INTEGER i,iff,ii,inext,inextp,k
	INTEGER mj,mk,ma(55)
C       REAL mj,mk,ma(55)
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
 11	   continue
	   do 13 k=1,4
	      do 12 i=1,55
		 ma(i)=ma(i)-ma(1+mod(i+30,55))
		 if(ma(i).lt.MZ)ma(i)=ma(i)+MBIG
 12	      continue
 13	   continue
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
	END
c------------------------------------------------------
	subroutine psp(g,tau,tm,vt,vpsp)
	implicit none
	real*8 g,tau,tm,vt,vpsp,tp
c       
	tp=log(tm/tau)/(1/tau-1/tm)
	vpsp=g*exp(-tp/tau)
	vpsp=20.0*vpsp/vt
c       
	return
	end
c-----------------------------------------------------
	
