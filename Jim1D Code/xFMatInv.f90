	Subroutine Matinv()
	use support
	real*8 rdum,cabs,cum,cnorm,cmult
	neq=2*nsc
	do i=1,neq
	cabs=-1.e9
	do j=1,ntbw
		cabs=max(cabs,abs(cm(i,j)))
	enddo
	do j=1,ntbw
		cm(i,j)=cm(i,j)/cabs
	enddo
	rhs(i)=rhs(i)/cabs
	enddo
!	do n=1,2*nsc
!		write(2,*)n,iloc(n)
!		write(2,'(5es14.5)')(cm(n,j),j=1,4),rhs(n)
!	enddo
! 	write(2,*)neq,ntbw
	!  lower diagonalization
		do i=1,neq-1
		! pivot
! 		  write(2,'(2i4,4(5es14.5,/))')i,iloc(i),(cm(i,k),k=1,ntbw),rhs(i)
			npiv=i+ntbw
			if(npiv>neq)npiv=neq
			cabs=abs(cm(i,1))
			do j=i+1,npiv
				if(iloc(j)==i)then
					if(abs(cm(j,1))>cabs)then !pivot
!					   write(2,'(2i4,4(5es14.5,/))')i,j,cm(j,1),cabs
						rdum=rhs(j)
						cabs=abs(cm(j,1))
						rhs(j)=rhs(i)
						rhs(i)=rdum
						do k=1,ntbw
							cum=cm(j,k)
							cm(j,k)=cm(i,k)
							cm(i,k)=cum
						enddo
					endif  !(abs(cm(j,i))>cabs)then !pivot
				endif !(iloc(j)==i)then
			enddo  !j=i+1,npiv
		!  normalize pivot row and move coefs to 'left'
!		  write(2,'(2i4,4(5es14.5,/))')i,iloc(i),(cm(i,k),k=1,ntbw),rhs(i)
			cnorm=cm(i,1)
			do k=2,ntbw
				cm(i,k-1)=cm(i,k)/cnorm
			enddo
			cm(i,ntbw)=0.
			rhs(i)=rhs(i)/cnorm
!		  write(2,'(2i4,4(5es14.5,/))')i,iloc(i),(cm(i,k),k=1,ntbw),rhs(i)
		!  lower diagnolization	 
			do j=i+1,npiv
				if(iloc(j)==i)then
!		   write(2,'(3i4,4(5es14.5,/))')i,j,iloc(j),(cm(j,k),k=1,ntbw),rhs(j)
						iloc(j)=iloc(j)+1
						cmult=-cm(j,1)
						do k=2,ntbw
							km=k-1
							cm(j,km)=cm(j,k)+cmult*cm(i,km)
						enddo
						rhs(j)=rhs(j)+cmult*rhs(i)
						cm(j,ntbw)=0.
!		   write(2,'(3i4,4(5es14.5,/))')i,j,iloc(j),(cm(j,k),k=1,ntbw),rhs(j)
				endif !(iloc(j)==i)then
			enddo  !j=i+1,npiv
		enddo !i=1,neq-1
	!  normalize last row
		rhs(neq)=rhs(neq)/cm(neq,1)
	!  back substitution
	do i=1,neq
!	write(2,'(i4,5es14.5)')i,(cm(i,j),j=1,4),rhs(i)
	enddo
		do i=neq-1,1,-1
			nlim=ntbw-1
			if(i+nlim>neq)nlim=neq-i
			do j=1,nlim
				rhs(i)=rhs(i)-rhs(i+j)*cm(i,j)
			enddo
		enddo !i=neq-1,1,-1
		  !write(2,'(4(5es14.5,/))')(rhs(i),i=1,neq)
	end Subroutine Matinv
