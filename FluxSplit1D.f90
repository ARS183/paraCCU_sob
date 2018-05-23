!==================================================================
!==============本程序提供一维Euler方程各种流通量分解结果===========	
!==================================================================


!*********************通用格式说明**********************************
!                                                                  *
!	Flux_??(N,p,u,rho,kappa,flux_plus,flux_minus)                *
!                                                                  *
!		1、p,u,rho为N维向量                                      *
!		2、flux_plus,flux_minus为3*N维矩阵					     *
!		3、kappa为比热比,通常取1.4，此处也作为输入参量			 *
!*******************************************************************	
	
	
!	subroutine FluxSplit_StegerWarming(N,p,u,rho,kappa, &
!     flux_plus,flux_minus)
!	implicit double precision (a-h,o-z)
!	double precision kappa,Ma
!	dimension p(N),u(N),rho(N)
!	dimension flux_plus(3,N),flux_minus(3,N)
!
!	do 100 i=1,N
!
!	a  = dsqrt(kappa*p(i)/rho(i))
!	Ma = u(i)/a
!
!	if (Ma.le.(-1.d0)) then
!	flux_plus(1,i)  = 0.d0
!	flux_plus(2,i)  = 0.d0
!	flux_plus(3,i)  = 0.d0
!	flux_minus(1,i) = rho(i)*u(i)
!	flux_minus(2,i) = rho(i)*u(i)*u(i)+p(i)
!	flux_minus(3,i) = 0.5d0*rho(i)*u(i)*u(i)*u(i) &
!     +kappa/(kappa-1.d0)*p(i)*u(i)
!
!	else if (Ma.le.0.d0 .and. Ma.gt.(-1.d0)) then
!	coef1 = 0.5d0*rho(i)*(u(i)+a)/kappa
!	coef2 = (kappa-1.d0)/kappa*rho(i)*u(i)
!	coef3 = 0.5d0*rho(i)*(u(i)-a)/kappa
!	flux_plus(1,i)  = coef1
!	flux_plus(2,i)  = coef1*(u(i)+a)
!	flux_plus(3,i)  = coef1*(u(i)*u(i)*0.5d0+a*a/(kappa-1.d0) &
!    +a*u(i))
!	flux_minus(1,i) = coef2+coef3
!	flux_minus(2,i) = coef2*u(i)+coef3*(u(i)-a)
!	flux_minus(3,i) = coef2*0.5d0*u(i)*u(i)+coef3*(u(i)*u(i)*0.5d0  &
!    +a*a/(kappa-1.d0)-a*u(i))
!
!	else if (Ma.le.1.d0 .and. Ma.gt.0.d0) then
!	coef1 = 0.5d0*rho(i)*(u(i)+a)/kappa
!	coef2 = (kappa-1.d0)/kappa*rho(i)*u(i)
!	coef3 = 0.5d0*rho(i)*(u(i)-a)/kappa
!	flux_plus(1,i)  = coef1+coef2
!	flux_plus(2,i)  = coef1*(u(i)+a)+coef2*u(i)
!	flux_plus(3,i)  = coef1*(u(i)*u(i)*0.5d0+a*a/(kappa-1.d0) &
!    +a*u(i))+coef2*0.5d0*u(i)*u(i)
!	flux_minus(1,i) = coef3
!	flux_minus(2,i) = coef3*(u(i)-a)
!	flux_minus(3,i) = coef3*(u(i)*u(i)*0.5d0 &
!    +a*a/(kappa-1.d0)-a*u(i))
!
!	else
!	flux_plus(1,i) = rho(i)*u(i)
!	flux_plus(2,i) = rho(i)*u(i)*u(i)+p(i)
!	flux_plus(3,i) = 0.5d0*rho(i)*u(i)*u(i)*u(i) &
!    +kappa/(kappa-1.d0)*p(i)*u(i)
!	flux_minus(1,i)  = 0.d0
!	flux_minus(2,i)  = 0.d0
!	flux_minus(3,i)  = 0.d0
!	
!	end if
!
!100	continue
!
!	end subroutine FluxSplit_StegerWarming

!--------------------------------------------------------
!----------------changed by chenjinqiang------------------ 
	subroutine FluxSplit_VanLeer(p,u,rho,kappa,flux_plus,flux_minus)
	include 'openNS3d.h'
	real(kind=OCFD_REAL_KIND)::kappa,Ma,a
	real(kind=OCFD_REAL_KIND) :: rho(1-LAP:nx+LAP),u(1-LAP:nx+LAP),p(1-LAP:nx+LAP)
	real(kind=OCFD_REAL_KIND) :: flux_plus(3,1-LAP:nx+LAP),flux_minus(3,1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND) :: dfdx_plus(3,1-LAP:nx+LAP),dfdx_minus(3,1-LAP:nx+LAP)

	
	do 100 i=1,nx

	a  = dsqrt(kappa*p(i)/rho(i))
	Ma = u(i)/a

	if (Ma.le.(-1.d0)) then
	flux_plus(1,i)  = 0.d0
	flux_plus(2,i)  = 0.d0
	flux_plus(3,i)  = 0.d0
	flux_minus(1,i) = rho(i)*u(i)
	flux_minus(2,i) = rho(i)*u(i)*u(i)+p(i)
	flux_minus(3,i) = 0.5d0*rho(i)*u(i)*u(i)*u(i)+kappa/(kappa-1.d0)*p(i)*u(i)

	else if (dabs(Ma).lt.1.d0) then
	coef1 = rho(i)*a*0.25d0*(Ma+1.d0)*(Ma+1.d0)
	coef2 = -rho(i)*a*0.25d0*(Ma-1.d0)*(Ma-1.d0)
	flux_plus(1,i)  = coef1
	flux_plus(2,i)  = coef1*((kappa-1.d0)*u(i)+2.d0*a)/kappa
	flux_plus(3,i)  = coef1*((kappa-1.d0)*u(i)+2.d0*a)**2/(2.d0*(kappa+1.d0)*(kappa-1.d0))
	flux_minus(1,i)  = coef2
	flux_minus(2,i)  = coef2*((kappa-1.d0)*u(i)-2.d0*a)/kappa
	flux_minus(3,i)  = coef2*((kappa-1.d0)*u(i)-2.d0*a)**2/(2.d0*(kappa+1.d0)*(kappa-1.d0))

	else
	flux_plus(1,i) = rho(i)*u(i)
	flux_plus(2,i) = rho(i)*u(i)*u(i)+p(i)
	flux_plus(3,i) = 0.5d0*rho(i)*u(i)*u(i)*u(i)+kappa/(kappa-1.d0)*p(i)*u(i)
	flux_minus(1,i)  = 0.d0
	flux_minus(2,i)  = 0.d0
	flux_minus(3,i)  = 0.d0

	end if

100	continue

	end subroutine FluxSplit_VanLeer

!--------------------------------------------------------

!	subroutine FluxSplit_LiouSteffen(N,p,u,rho,kappa, &
!    flux_plus,flux_minus)
!	implicit double precision (a-h,o-z)
!	double precision kappa,Ma,Ma_plus,Ma_minus
!	dimension p(N),u(N),rho(N)
!	dimension flux_plus(3,N),flux_minus(3,N)
!
!	do 100 i=1,N
!
!	a  = dsqrt(kappa*p(i)/rho(i))
!	Ma = u(i)/a
!
!	if (Ma.le.(-1.d0)) then
!	Ma_plus  = 0.d0
!	Ma_minus = Ma
!	p_plus   = 0.d0
!	p_minus  = p(i)
!
!	else if (dabs(Ma).lt.1.d0) then
!	Ma_plus  = 0.25d0*(1.d0+Ma)*(1.d0+Ma)
!	Ma_minus = -0.25d0*(Ma-1.d0)*(Ma-1.d0)
!	p_plus   = 0.5d0*(1.d0+Ma)*p(i)
!	p_minus  = 0.5d0*(1.d0-Ma)*p(i)
!
!	else
!	Ma_plus  = Ma
!	Ma_minus = 0.d0
!	p_plus   = p(i)
!	p_minus  = 0.d0 
!	
!	end if
!
!	coef1 = max(0.d0,Ma_plus+Ma_minus)
!	coef2 = min(0.d0,Ma_plus+Ma_minus)
!	ht    = 0.5d0*u(i)*u(i)+1.d0/(kappa-1.d0)*a*a
!
!	flux_plus(1,i)  = coef1*rho(i)*a
!	flux_plus(2,i)  = coef1*rho(i)*u(i)*a+p_plus
!	flux_plus(3,i)  = coef1*rho(i)*ht*a
!	flux_minus(1,i) = coef2*rho(i)*a
!	flux_minus(2,i) = coef2*rho(i)*u(i)*a+p_minus
!	flux_minus(3,i) = coef2*rho(i)*ht*a
!
!100	continue
!
!	end subroutine FluxSplit_LiouSteffen
!
!---------------------------------------------
!
	

!	subroutine FluxSplit_ZhaBilgen(N,p,u,rho,kappa, &
!    flux_plus,flux_minus)
!	implicit double precision (a-h,o-z)
!	double precision kappa,Ma
!	dimension p(N),u(N),rho(N)
!	dimension flux_plus(3,N),flux_minus(3,N)
!	
!	do 100 i=1,N
!
!	a  = dsqrt(kappa*p(i)/rho(i))
!	Ma = u(i)/a
!
!	if (Ma.le.(-1.d0)) then
!	pu_plus  = 0.d0
!	pu_minus = p(i)*u(i)
!	p_plus   = 0.d0
!	p_minus  = p(i)
!
!	else if (dabs(Ma).lt.1.d0) then
!	pu_plus  = 0.5d0*p(i)*(u(i)+a)
!	pu_minus = 0.5d0*p(i)*(u(i)-a)
!	p_plus   = 0.5d0*(1.d0+Ma)*p(i)
!	p_minus  = 0.5d0*(1.d0-Ma)*p(i)
!
!	else
!	pu_plus  = p(i)*u(i)
!	pu_minus = 0.d0
!	p_plus   = p(i)
!	p_minus  = 0.d0 
!	
!	end if
!
!	coef1 = max(0.d0,u(i))
!	coef2 = min(0.d0,u(i))
!	et    = 0.5d0*u(i)*u(i)+p(i)/rho(i)/(kappa-1.d0)
!
!	flux_plus(1,i)  = coef1*rho(i)
!	flux_plus(2,i)  = coef1*rho(i)*u(i)+p_plus
!	flux_plus(3,i)  = coef1*rho(i)*et+pu_plus
!	flux_minus(1,i) = coef2*rho(i)
!	flux_minus(2,i) = coef2*rho(i)*u(i)+p_minus
!	flux_minus(3,i) = coef2*rho(i)*et+pu_minus
!100	continue	
!
!	end subroutine FluxSplit_ZhaBilgen
