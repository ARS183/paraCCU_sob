!
!-------------------------Parallel------------------------
!

subroutine Du2Dx_PADE4(u,d2f)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

    if (npx==0) then
        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,1)
        call D2F_PADE4(u,d2f,nx,hx)

    elseif (npx==npx0-1) then
        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,2)
        call D2F_PADE4(u,d2f,nx,hx)

    else
        call OCFD_D2F_SB_PADE4(u,d2f,nx,hx,0)
        call D2F_PADE4(u,d2f,nx,hx)

    endif

!    deallocate(d2f,df)

end subroutine Du2Dx_PADE4

!
!------------------------------------------------------
!




!
!===============´®ÐÐ===================
!
subroutine Du2Dx_PADE4_serial(u,d2f)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,1)

        call OCFD_D2F_BOUND_PADE4(u,d2f,nx,hx,2)

        call D2F_PADE4(u,d2f,nx,hx)

!    deallocate(d2f,df)

end subroutine Du2Dx_PADE4_serial
!
!========================================
!




!
!-------------------------------Parallel----------------------------
!

subroutine DuDx_UCC_UpWind(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

    if (npx==0) then
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
        call DF_UCC45_P(u,d2f,df,nx,hx,1)

    elseif (npx==npx0-1) then
        call OCFD_DF_SB2_UCC45_P(u,df,d2f,nx,hx)
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
        call DF_UCC45_P(u,d2f,df,nx,hx,2)

    else
        call OCFD_DF_SB2_UCC45_P(u,df,d2f,nx,hx)
        call DF_UCC45_P(u,d2f,df,nx,hx,0)

    endif

!    deallocate(d2f,df)

end subroutine DuDx_UCC_UpWind



subroutine DuDx_UCC_DownWind(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

    if (npx==0) then

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)
        call OCFD_DF_SB2_UCC45_M(u,df,d2f,nx,hx)
        call DF_UCC45_M(u,d2f,df,nx,hx,1)

    elseif (npx==npx0-1) then

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
        call DF_UCC45_M(u,d2f,df,nx,hx,2)

    else
        call OCFD_DF_SB2_UCC45_M(u,df,d2f,nx,hx)
        call DF_UCC45_M(u,d2f,df,nx,hx,0)
        
    endif

!    deallocate(d2f,df)

end subroutine DuDx_UCC_DownWind

!
!------------------------------------------------------------
!




!
!============================Serial==========================
!
subroutine DuDx_UCC_UpWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

 
        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)

        call DF_UCC45_P(u,d2f,df,nx,hx,12)

!    deallocate(d2f,df)

end subroutine DuDx_UCC_UpWind_serial



subroutine DuDx_UCC_DownWind_serial(u,d2f,df)
    include 'openNS3d.h'		
!	real(kind=OCFD_REAL_KIND)::hx
    real(kind=OCFD_REAL_KIND)::u(1-LAP:nx+LAP),d2f(1-LAP:nx+LAP)
    real(kind=OCFD_REAL_KIND)::df(1-LAP:nx+LAP)
!	integer :: nx

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,1)

        call OCFD_DF_BOUND_UCC45(u,df,nx,hx,2)
  
        call DF_UCC45_M(u,d2f,df,nx,hx,12)


!    deallocate(d2f,df)

end subroutine DuDx_UCC_DownWind_serial

!
!============================================================
!