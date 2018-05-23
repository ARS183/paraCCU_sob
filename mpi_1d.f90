!send and recv mesg, to check array in x direction.
subroutine check_x1d(f)
	include 'openNS3d.h'
	integer i,j,k1,npx1,npx2,my_mod1
	real(kind=OCFD_REAL_KIND):: f(1-LAP:nx+LAP),  &
		    tmp_send1(1:LAP),tmp_send2(1:LAP), tmp_recv1(1:LAP),tmp_recv2(1:LAP)
	  
		
		do i=1,LAP
		k1=i
		tmp_send1(k1)=f(i)
		tmp_send2(k1)=f(nx-LAP+i)
		enddo
    
!!

		call MPI_Sendrecv(tmp_send1,LAP,OCFD_DATA_TYPE, ID_XM1, 9000, &
		    tmp_recv2, LAP,  OCFD_DATA_TYPE,ID_XP1, 9000,MPI_COMM_WORLD,Status,ierr)
		call MPI_Sendrecv(tmp_send2,LAP, OCFD_DATA_TYPE,ID_XP1,8000,   &
		    tmp_recv1,LAP, OCFD_DATA_TYPE,ID_XM1,  8000,MPI_COMM_WORLD,Status,ierr)

		if(ID_XM1 .ne. MPI_PROC_NULL) then
		 do i=1,LAP
		   k1=i
		   f(i-LAP)=tmp_recv1(k1)
		 enddo
		endif

		if(ID_XP1 .ne. MPI_PROC_NULL) then
		 do i=1,LAP
		   k1=i
		   f(nx+i)=tmp_recv2(k1)
		 enddo
        endif
        
end subroutine