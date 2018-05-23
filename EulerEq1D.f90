	include "openNS3d_DF_AECDS.f90"
    include "DuDx_UCC.f90"
    include "mpi_1d.f90"
	include "openNS3d_DF_BOUND.f90"
	include "FluxSplit1D.f90"
	include "Output_Tecplot.f90"

	program EulerEq1D
	!implicit double precision (a-h,o-z)
	implicit none
	include 'openNS3d.h'
	!parameter (N=51)
	real(kind=OCFD_REAL_KIND) :: kappa,xbeg,xend
	!dimension Ma(N)	
	!dimension rho(N),u(N),p(N)
	!dimension a(N),S(N),x(N)
	real(kind=OCFD_REAL_KIND),allocatable :: Ma(:),rho(:),u(:),p(:)
	real(kind=OCFD_REAL_KIND),allocatable :: a(:),S(:),x(:),R(:)
	real(kind=OCFD_REAL_KIND),allocatable :: Ma_global(:),rho_global(:),u_global(:)
	real(kind=OCFD_REAL_KIND),allocatable :: a_global(:),S_global(:),p_global(:)
	real(kind=OCFD_REAL_KIND),allocatable :: xx(:),temp_x(:)
	real(kind=OCFD_REAL_KIND) :: time_start,time_end,ctime
	integer :: ka,k,i,i_global,ii,npx1,npx2,my_mod1,Nt



	call MPI_INIT(ierr)
	call MPI_COMM_SIZE(MPI_COMM_WORLD,np_size,ierr)
	call MPI_COMM_RANK(MPI_COMM_WORLD,my_id,ierr)	

	nx_global=601

	kappa=1.4d0

	xbeg=-10.d0
	xend=10.d0
	slx=xend-xbeg
	
	hx=slx/dble(nx_global-1)

	T = .025d0
	!dt=4.276d-4/1.d0
	dt=2.138d-4/6.d0
	Nt=T/dt

!---------------------------------------------------------
!call part3d()
	npx0=np_size
	npx=mod(my_id,npx0)
	nx=nx_global/npx0
	if(npx .lt. mod(nx_global,npx0)) nx=nx+1

	do k=0,npx0-1
    	ka=min(k,mod(nx_global,npx0))
    	i_offset(k)=int(nx_global/npx0)*k+ka+1
    	i_nn(k)=nx_global/npx0
    	if(k .lt. mod(nx_global,npx0)) i_nn(k)=i_nn(k)+1
	enddo
	npx1=my_mod1(npx-1,npx0)
	npx2=my_mod1(npx+1,npx0)
	ID_XM1=npx1    ! -1 proc in x-direction
	ID_XP1=npx2 
	if(npx .eq. 0) ID_XM1=MPI_PROC_NULL     ! if not periodic, 0 node donot send mesg to npx0-1 node
	if(npx .eq. npx0-1) ID_XP1=MPI_PROC_NULL
!---------------------------------------------------------
	allocate(x(1-LAP:nx+LAP),u(1-LAP:nx+LAP))
	allocate(rho(1-LAP:nx+LAP),p(1-LAP:nx+LAP))
	allocate(a(1-LAP:nx+LAP),S(1-LAP:nx+LAP))
	allocate(R(1-LAP:nx+LAP),Ma(1-LAP:nx+LAP))
!---------------------------------------------------------
!call define_grid(x)
if (my_id .eq. 0) then

    allocate(xx(1:nx_global))
    do i=1,nx_global
!        xx(i)=dble(i-1)*hx
		xx(i)=xbeg+dble(i-1)*hx
    enddo
    
    do 100 i=0,npx0-1		
		ka=i
			
		if (ka .eq. 0) then
			do ii=1,nx
					i_global=i_offset(i)-1+ii
					x(ii)=xx(i_global)
			enddo

		else
			allocate(temp_x(1:i_nn(i)))
			do ii=1,i_nn(i)
				i_global=i_offset(i)+ii-1
				temp_x(ii)=xx(i_global)
			enddo
					
	   	    call MPI_SEND(temp_x,i_nn(i),OCFD_DATA_TYPE,  &
			ka,1,MPI_COMM_WORLD,ierr) 

		    deallocate(temp_x)
		endif
100	continue

!	deallocate(xx)
	
	else
	    allocate(temp_x(1:i_nn(npx)))

		call MPI_RECV(temp_x,i_nn(npx),OCFD_DATA_TYPE,  &
        0,1,MPI_COMM_WORLD,status,ierr) 
        
		do ii=1,i_nn(npx)
			x(ii)=temp_x(ii)
		enddo	

		deallocate(temp_x)

endif
!---------------------------------------------------------

!-----------initial & boundary condition----------------
!	call test_initial(N,p,u,rho)
	do 101 i=1,nx
	if (x(i).ge.0.d0) then

!	rho(i) = 0.125d0
!	u(i)   = 0.d0
!	p(i)   = 1.d4

	rho(i) = 0.5d0
	u(i)   = 0.d0
	p(i)   = 5710.d0

	else
!	rho(i) = 1.d0
!	u(i)   = 0.d0
!	p(i)   = 1.d5

	rho(i) = 0.445d0
	u(i)   = 0.698d0
	p(i)   = 35280.d0

	end if
101	continue
!-----------------------程序主体--------------------------------

	do 500 i=1,Nt
	write(*,*) i
	call RK4(p,u,rho,kappa)
500	continue

	do 901 i=1,nx
	a(i)  = dsqrt(kappa*p(i)/rho(i))
	Ma(i) = u(i)/a(i)
	S(i)  = 717.d0*dlog(p(i))-1004.d0*dlog(rho(i))	!entropy
901	continue


	call MPI_BARRIER(MPI_COMM_WORLD,ierr)
!
!-------------------------------------------------------------
!
	allocate(Ma_global(1:nx_global),u_global(1:nx_global))
	allocate(rho_global(1:nx_global),p_global(1:nx_global))
	allocate(a_global(1:nx_global),S_global(1:nx_global))

	call write_data_global(u,u_global)
	call write_data_global(rho,rho_global)
	call write_data_global(p,p_global)
	call write_data_global(Ma,Ma_global)
	call write_data_global(a,a_global)
	call write_data_global(S,S_global)

	if  (my_id .ne. 0) then
		deallocate(u_global,rho_global,p_global,Ma_global,a_global,S_global)
	endif

!-------------------------------------------------------------
	if (my_id .eq. 0) then
		call Toplt1D(xx,u_global,'u',1,nx_global)
		call Toplt1D(xx,rho_global,'rho',3,nx_global)
		call Toplt1D(xx,p_global,'p',1,nx_global)
		call Toplt1D(xx,Ma_global,'Ma',2,nx_global)
		call Toplt1D(xx,a_global,'a',1,nx_global)
		call Toplt1D(xx,S_global,'S',1,nx_global)
	endif



	call MPI_FINALIZE(ierr)

	end program EulerEq1D

!
!
!
!===========================================================================
!
!

subroutine ComputeR(p,u,rho,kappa,R)	!采用RK算法的右端项
!	implicit double precision (a-h,o-z)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND)::kappa
	real(kind=OCFD_REAL_KIND) :: rho(1-LAP:nx+LAP),u(1-LAP:nx+LAP),p(1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: flux_plus(3,1-LAP:nx+LAP),flux_minus(3,1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: d2f_plus(3,1-LAP:nx+LAP),d2f_minus(3,1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND) :: dfdx_plus(3,1-LAP:nx+LAP),dfdx_minus(3,1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: R(3,1-LAP:nx+LAP)

	
!	call FluxSplit_StegerWarming(N,p,u,rho,kappa,
!    &	flux_plus,flux_minus)
	call FluxSplit_VanLeer(p,u,rho,kappa, &
    flux_plus,flux_minus)
!	call FluxSplit_LiouSteffen(N,p,u,rho,kappa,
!     &	flux_plus,flux_minus)
!	call FluxSplit_ZhaBilgen(N,p,u,rho,kappa,
!    &	flux_plus,flux_minus)

	do 200 j=1,3

	call check_x1d(flux_plus(j,:))
	call check_x1d(flux_minus(j,:))

	call Du2Dx_PADE4(flux_plus(j,:),d2f_plus(j,:))
	call Du2Dx_PADE4(flux_minus(j,:),d2f_minus(j,:))

	call DuDx_UCC_UpWind(flux_plus(j,:),d2f_plus(j,:),dfdx_plus(j,:))
	call DuDx_UCC_DownWind(flux_minus(j,:),d2f_minus(j,:),dfdx_minus(j,:))


	R(j,:) = -(dfdx_plus(j,:)+dfdx_minus(j,:))

200	continue

!-----------------NND-------------
!	do 300 i=1,3
!	do 300 j=2,N-1
!	if (j .eq. 2) then
!	dfdx_plus(i,j-1)=dfdx_plus(i,j)
!	dfdx_minus(i,j-1)=dfdx_minus(i,j)
!	else if (j .eq. N-1) then
!	dfdx_plus(i,j+1)=dfdx_plus(i,j)
!	dfdx_minus(i,j+1)=dfdx_minus(i,j)
!	end if
!
!	call limiter1(dfdx_plus(i,j),dfdx_plus(i,j+1),temp1)
!	call limiter1(dfdx_plus(i,j-1),dfdx_plus(i,j),temp2)
!	call limiter1(dfdx_minus(i,j),dfdx_minus(i,j+1),temp3)
!	call limiter1(dfdx_minus(i,j-1),dfdx_minus(i,j),temp4)
!	R(i,j) = -(dfdx_plus(i,j)+dfdx_minus(i,j)
!    &	+0.5d0*(temp1-temp2)+0.5d0*(-temp3+temp4))	
!300	continue
!-------------------先不考虑限制器了-----------------------	

end subroutine ComputeR











	

	!subroutine ComputeR1(N,p,u,rho,kappa,R,dx)
	!implicit double precision (a-h,o-z)
	!double precision kappa
	!dimension H_plus(3,N),H_minus(3,N)
	!dimension flux_plus(3,N),flux_minus(3,N)
	!dimension R(3,N)
 !
 !
	!call FluxSplit_VanLeer(N,p,u,rho,kappa, &
 !   flux_plus,flux_minus)
	!call ComputeH(N,flux_plus,flux_minus,dx,H_plus,H_minus)
	!call ComputeRR(N,H_plus,H_minus,R,dx)
 !
	!end subroutine ComputeR1
!
!
!==========================================================================
!
!


!----------------------------------------------------------
!      4阶RK格式,bychenjinqiang,201804
!----------------------------------------------------------
subroutine RK4(p,u,rho,kappa)	!4阶RK格式,待编
!	implicit double precision (a-h,o-z)
	include 'openNS3d.h'
    real(kind=OCFD_REAL_KIND) :: rho(1-LAP:nx+LAP),u(1-LAP:nx+LAP),p(1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: R(3,1-LAP:nx+LAP),Q(3,1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND) :: R1(3,1-LAP:nx+LAP),R2(3,1-LAP:nx+LAP),R3(3,1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND) :: Q1(3,1-LAP:nx+LAP),Q2(3,1-LAP:nx+LAP),Q3(3,1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: kappa

	
	call ComputeR(p,u,rho,kappa,R)





	do 100 i=1,nx
	Q(1,i) = rho(i)
	Q(2,i) = rho(i)*u(i)
	Q(3,i) = 0.5d0*rho(i)*u(i)*u(i)+p(i)/(kappa-1.d0)
100   continue	



	do 200 j=1,3
	Q1(j,1:nx) = Q(j,1:nx)+(dt/2.d0)*R(j,1:nx)
200   continue

	do 300 i=1,nx
	rho(i) = Q1(1,i)
	u(i)   = Q1(2,i)/rho(i)
	p(i)   = (Q1(3,i)-0.5d0*rho(i)*u(i)*u(i))*(kappa-1.d0)
300   continue	
	
	call ComputeR(p,u,rho,kappa,R1)

	do 400 j=1,3
	Q2(j,:) = Q(j,:)+0.5d0*dt*R1(j,:)
400   continue
	
	do 500 i=1,nx
	rho(i) = Q2(1,i)
	u(i)   = Q2(2,i)/rho(i)
	p(i)   = (Q2(3,i)-0.5d0*rho(i)*u(i)*u(i))*(kappa-1.d0)
500   continue

	call ComputeR(p,u,rho,kappa,R2)

	do 600 j=1,3
	Q3(j,:) = Q(j,:)+dt*R2(j,:)
600   continue
	
	do 700 i=1,nx
	rho(i) = Q3(1,i)
	u(i)   = Q3(2,i)/rho(i)
	p(i)   = (Q3(3,i)-0.5d0*rho(i)*u(i)*u(i))*(kappa-1.d0)
700   continue

	call ComputeR(p,u,rho,kappa,R3)

	do 800 j=1,3
	Q(j,:) = Q(j,:)+(dt/6.d0)*(R(j,:)+2.d0*R1(j,:)+2.d0*R2(j,:)+R3(j,:))

800   continue
	
	do 900 i=1,nx
	rho(i) = Q(1,i)
	u(i)   = Q(2,i)/rho(i)
	p(i)   = (Q(3,i)-0.5d0*rho(i)*u(i)*u(i))*(kappa-1.d0)
900   continue

	end subroutine RK4







!-------------------------------------------------------------
!数据汇总
!-------------------------------------------------------------

subroutine write_data_global(phi,phi_global)
	include 'openNS3d.h'
	integer :: i_global,ii,id_mid,ka,i,count,inum,iter
	real(kind=OCFD_REAL_KIND),dimension(1-LAP:nx+LAP) :: phi
	real(kind=OCFD_REAL_KIND),allocatable,dimension(:) :: temp_phi
	real(kind=OCFD_REAL_KIND),dimension(1:nx_global) :: phi_global

if (my_id .eq. 0) then
    
    do 100 i=0,npx0-1		
		ka=i
			
		if (ka .eq. 0) then
			do ii=1,nx
					i_global=i_offset(i)-1+ii
					phi_global(i_global)=phi(ii)
			enddo

		else
			allocate(temp_phi(1:i_nn(i)))

			call MPI_RECV(temp_phi,i_nn(i),OCFD_DATA_TYPE,  &
			ka,1,MPI_COMM_WORLD,status,ierr) 
			
			do ii=1,i_nn(i)
				i_global=i_offset(i)-1+ii
				phi_global(i_global)=temp_phi(ii)
			enddo	
	
			deallocate(temp_phi)

		endif
100	continue
	
else
	allocate(temp_phi(1:i_nn(npx)))
	do ii=1,i_nn(npx)

		temp_phi(ii)=phi(ii)

	enddo
				
	call MPI_SEND(temp_phi,i_nn(npx),OCFD_DATA_TYPE,  &
	0,1,MPI_COMM_WORLD,ierr) 

	deallocate(temp_phi)

endif
end subroutine write_data_global	



function my_mod1(i,n)
          
    integer my_mod1,i,n
    if(i.lt.0) then
        my_mod1=i+n
    else if (i.gt.n-1) then
        my_mod1=i-n
    else
        my_mod1=i
    endif
end