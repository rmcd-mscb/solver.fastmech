Subroutine StdyStt(ich,itp,qus,zds)
	!    curxsc(ns,1)=ws elevation
	!    curxsc(ns,2)=xsec area
	!    curxsc(ns,3)=hydraulic radius =area/wet perim
	!    curxsc(ns,4)=hydraulic depth =area/top width
	!    curxsc(ns,5)=top width
	!    curxsc(ns,6)=conveyance
!	!    curxsc(ns,7)=discharge
!  xsloc(i,5)=emin  ! xsloc( ,5) set to min ch elev
	use support
	nus=1
	if(ich>1)nus=nsce(ich-1)+1
	nds=nsce(ich)
	select case(itp)
	case(1)
! itp==1, discharge-stage stepbw
 	call StepBW(1,ich,qus,zds)
	case(2)
! itp==2, stage-stage stepbw, solve for flow in channel, qus==us stage
	usstg=qus
	if(usstg==zds)then ! for usstg>zds
		do i=nus,nds
		call Farea(i,usstg)
		curxsc(i,7)=0.
		enddo
	else ! ! for usstg/=zds
		if(usstg>zds)then
			call Farea(nus,usstg)
			qcrit=sqrt(gg*curxsc(nus,2)**3/curxsc(nus,5))
			indrcn=1
			bcstg=zds
			stgfct=usstg
			isstg=nus
		else !if(usstg>zds)then
			call Farea(nds,zds)
			qcrit=sqrt(gg*curxsc(nds,2)**3/curxsc(nds,5))
			indrcn=-1
			bcstg=usstg
			stgfct=zds
			isstg=nds
		endif !(usstg>zds)then
	  qfct=qcrit*.95
	  tol=.001*curxsc(nus,4)
	  iter=0
 	  call StepBw(indrcn,ich,qfct,bcstg)
	  fo=	curxsc(isstg,1)-stgfct
	  dq=.01*qfct
	  f=1.e9
! 	  write(2,*)indrcn,isstg,usstg,zds,stgfct,fo,qfct,dq
	do while(abs(f)>tol.and.iter<100)
		iter=iter+1
		if(qfct+dq>qcrit)then
			dq=qcrit-qfct
		elseif(qfct+dq<0.)then
			dq=-qfct+tol
		endif
		qfct=qfct+dq
		call StepBw(indrcn,ich,qfct,bcstg)
!      		write(2,*)iter,qfct,f,curxsc(nus,1),dq
		f=curxsc(isstg,1)-stgfct
		if(f==fo)return
		dq=-.1*f*dq/(f-fo) !
		fo=f
	enddo
	endif !(usstg==zds)then ! for usstg>zds
	end select
end	Subroutine StdyStt
subroutine StepBW(inddrcn,nc,q,zds)
!  if inddrcn==1, usual step bw, compute ds>>us 
!  if inddrcn==-1, ds wselev>us, compute us>>ds 
use support
real kds,kus

	oflwcrtcl(nc)=1
	nfrst=1
	if(nc>1) nfrst=nsce(nc-1)+1
	nlst=nsce(nc)
	nus=nfrst
	nds=nlst
	instp=-1
	if(inddrcn<0)then
		nus=nlst
		nds=nfrst
		instp=1
	endif !(inddrcn<0)then
	sfdnom=q*abs(q)*2.
	atrm=2.*q*q/gg
	call zCrit(nds,abs(q),zc)
	if(zc>zds)then
! 		print *,'using critical downstream elev',nc,nds,zc,zds,q
! 		write(2,1001)nc,nds,zc,zds,q
		oflwcrtcl(nc)=-1
1001	format('Using critical downstream elev ',2i6,3f10.3)
	endif
	curxsc(nds,1)=max(zds,zc+.001)
	zs=curxsc(nds,1)
	call Farea(nds,zs)
	ads=curxsc(nds,2)
	kds=curxsc(nds,6)**2
	curxsc(nds,7)=q*real(inddrcn)
	do ns=nds+instp,nus,instp
		nsp=ns+1
		if(inddrcn<0)nsp=ns
		dz=.1
		call dxc(nsp,dx)
		call zCrit(ns,abs(q),zc)
		call Farea(ns,zc)
		fcrit=zs-zc+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
! 	print *,ns,fcrit,zs,dx*q*q/kds,q/sqrt(GG*curxsc(ns,2)**3/curxsc(ns,5))
		if(fcrit<=0)then  ! there is no subcritical upstream solution
			z=zc+.001
			call Farea(ns,z)
!		print *,'Supercritical flow solution'
! 		write(2,1002)nc,ns,zc,zs,q
			oflwcrtcl(nc)=-1
1002		format('Supercritical flow solution',2i6,3f10.3)
			goto 1000
		endif !(fcrit<0)then  ! there is no subcritical upstream solution
		z=max(zc+dz,curxsc(ns-instp,1))
		call Farea(ns,z)
		call dxc(nsp,dx)
!  	print *,ns,q,ffroud,z
		zprec=curxsc(ns,4)  !hd
		fo=zs-z+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
		iter=0
		f=1.e9
		do while(abs(f)>.001*zprec.and.iter<100)
			iter=iter+1
			if(z+dz<zc)then
				dz=zc+.001-z
			endif
			z=z+dz
			call Farea(ns,z)
			call dxc(nsp,dx)
			f=zs-z+dx*sfdnom/(kds+curxsc(ns,6)**2)+atrm*(1./ads-1./curxsc(ns,2))/(ads+curxsc(ns,2))
			if(f==fo)goto 1000
			dz=-f*dz/(f-fo)
			fo=f
		enddo
1000	continue
!	print *,ns,q/sqrt(GG*curxsc(ns,2)**3/curxsc(ns,5))
		zs=z
		curxsc(ns,7)=q*real(inddrcn)
		ads=curxsc(ns,2)
		kds=curxsc(ns,6)**2
	enddo  !section loop, ns
end subroutine StepBW
subroutine zCrit(is,q,zc)
use Support
	zsv=curxsc(is,1)
	qtol=.001*q
	qtol=min(qtol,.001)
	zc=xsloc(is,5)+.01 !etol
	f=-1.e9
	dz=0.1
	do while (f<0.)
		zc=zc+dz
		call Farea(is,zc)
		if(curxsc(is,2)<=0..or.curxsc(is,5)<=0.)then
			print *,' 1 area negative in zCrit',is,curxsc(is,1),curxsc(is,2),curxsc(is,5)
			stop
		endif
		f=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
	enddo 
	call Farea(is,zc)
	fo=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
	f=1.e9
	iter=0
	do while (abs(f)>qtol.and.iter<100)
		iter=iter+1
		if(zc+dz<xsloc(is,5)+.01)then  !etol
			dz=xsloc(is,5)-zc+.01 !etol
		endif
		zc=zc+dz
		call Farea(is,zc)
		if(curxsc(is,2)<=0..or.curxsc(is,5)<=0.)then
			print *,'f area negative in zCrit',is,curxsc(is,1),curxsc(is,2),curxsc(is,5)
			stop
		endif
! 		print *,is,iter,zc,f,xsloc(is,5),curxsc(is,2)
		f=rtgg*sqrt(curxsc(is,2)**3/curxsc(is,5))-q
		if(f==fo)return
		dz=-f*dz/(f-fo)
		fo=f
	enddo
end subroutine zCrit
subroutine MinCourTS(ich,tsmin)
	use support
	Call ChLims(ich,nus,nds)
	ifuck=0
	do i=nus,nds
	if(curxsc(i,2)<=0..or.curxsc(i,5)<=0.)then
	ifuck=1
	endif
	enddo
	if(ifuck==1)then
	do i=nus,nds
!	write(2,'(i4,7f10.3)')i,(curxsc(i,j),j=1,7)
	enddo
	stop
	endif
	tsmin=1.e9
	do i=nus+1,nds
		call dxc(i,dx)
		vel=abs(curxsc(i,7)/curxsc(i,2))
		c=sqrt(gg*curxsc(i,2)/curxsc(i,5))
		tsmin=min(tsmin,dx/(vel+c))
!		print *,i,dx,vel,c
	enddo
end subroutine MinCourTS
subroutine CheckFlow(ich,elus,elds)
	use support
	eltry=ellim(ich)+etol
	if( eltry<elus.or.eltry<elds)return
	call ChLims(ich,nus,nds)
	flwexst(ich)=-1
	do i=nus,nds
		call Farea(i,xsloc(i,5)+.01) !+etol)
		curxsc(i,7)=0.
	enddo
end subroutine CheckFlow
