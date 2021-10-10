MODULE Parallel_Finest_Grid_Smoothing_02
USE    OMP_LIB
USE    Parallel_Finest_Grid_Smoothing_01

PRIVATE

PUBLIC :: MSSize, Starting_Guess_and_Boundary_Conditions
PUBLIC :: Vanka_type_smoother_FG,Convergence_Test

CONTAINS

SUBROUTINE   MSSize    
 allocate(      U(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2))                                 ! 001
 allocate(      F(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2))                                 ! 002  
 allocate( BC_0YZ(           0:NYY_FG+2,0:NZZ_FG+2))                                 ! 003 
 allocate( BC_1YZ(           0:NYY_FG+2,0:NZZ_FG+2))                                 ! 004
 allocate( BC_X0Z(0:NXX_FG+2,           0:NZZ_FG+2))                                 ! 005
 allocate( BC_X1Z(0:NXX_FG+2,           0:NZZ_FG+2))                                 ! 006
 allocate( BC_XY0(0:NXX_FG+2,0:NYY_FG+2           ))                                 ! 007
 allocate( BC_XY1(0:NXX_FG+2,0:NYY_FG+2           ))                                 ! 008
END SUBROUTINE   MSSize  


Subroutine Starting_Guess_and_Boundary_Conditions
real       x,y,z,h
                 u = 0.d0;       f = 0.d0                        ! 001
                 h = a_**2+b_**2+c_**2                           ! 002 
           Res_000 = dexp(a_+b_+c_)*h                            ! 003
           Err_000 = dexp(a_+b_+c_)                              ! 004
      do         k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 005             
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 006            
        do   i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 007
           F(i,j,k)=-h*dexp(a_*x+b_*y+c_*z)                      ! 008
        end do                                                   ! 009
       end do                                                    ! 010
      end do                                                     ! 011
!= = = = = = = = = = = = =  X planes  = = = = = = = = = = = = =
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 012 
        do       k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 013 
      BC_0YZ(  j,k)= dexp(     b_*y+c_*z)                        ! 014 
      BC_1YZ(  j,k)= dexp(a_  +b_*y+c_*z)                        ! 015  
        end do                                                   ! 016 
       end do                                                    ! 017 
!= = = = = = = = = = = = =  Y planes  = = = = = = = = = = = = =   
      do     i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 018 
        do       k = 1,NZZ_FG+1;  z = dfloat(k-1)/dfloat(NZZ_FG) ! 019    
      BC_X0Z(i,  k)= dexp(a_*x     +c_*z)                        ! 020    
      BC_X1Z(i,  k)= dexp(a_*x+b_  +c_*z)                        ! 021               
        end do                                                   ! 022    
      end do                                                     ! 023    
!= = = = = = = = = = = = =  Z planes  = = = = = = = = = = = = =      
      do     i     = 1,NXX_FG+1;  x = dfloat(i-1)/dfloat(NXX_FG) ! 024    
       do      j   = 1,NYY_FG+1;  y = dfloat(j-1)/dfloat(NYY_FG) ! 025    
      BC_XY0(i,j  )= dexp(a_*x+b_*y     )                        ! 026
      BC_XY1(i,j  )= dexp(a_*x+b_*y+c_  )                        ! 027           
       end do                                                    ! 028
      end do                                                     ! 029
end Subroutine Starting_Guess_and_Boundary_Conditions


subroutine Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)
INTEGER    way
REAL*8     A_W,A_E,A_S,A_N,A_D,A_U,A_P,h_x,h_y_h_z 
INTEGER,   ALLOCATABLE ::  M_A_P(:,:,:)
REAL*8 ,   ALLOCATABLE ::     AG(:,:),  xG(:)   

                N_ = NX_block*NY_block*NZ_block                                     ! 001   
 allocate(   AG(N_,N_+1), xG(N_))                                                   ! 002
 allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))                               ! 003
                AG = 0.D0;         xG = 0.D0                                        ! 004 
               h_x = 1.d0/dfloat(NXX_FG)
               h_y = 1.d0/dfloat(NYY_FG)
               h_z = 1.d0/dfloat(NZZ_FG)
                
                                   m  = 0                                           ! 005
   do           k  = K_top,K_bot                                                    ! 006                                   
    do        j    = J_top,J_bot                                                    ! 007  
     do     i      = I_top,I_bot;  m  = m+1                                         ! 008
      M_A_P(i,j,k) =               m                                                ! 009    
     end do                                                                         ! 010
    end do                                                                          ! 011 
   end do                                                                           ! 012   
                                   N_ = m                                           ! 013
                                   m  = 0                                           ! 014
   do           k  = K_top,K_bot                                                    ! 015                                   
    do        j    = J_top,J_bot                                                    ! 016         
     do     i      = I_top,I_bot                                                    ! 017
                                   m  = m+1                                         ! 018   
          A_W = 0.D0; A_E = 0.D0; A_S = 0.D0; A_N = 0.D0                            ! 019
          A_D = 0.D0; A_U = 0.D0; A_P = 0.D0;                         way = 0       ! 020
! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(i==1)                                                                 then ! 021
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 022 
                             AG(m,N_+1)         = BC_0YZ(j,k)                       ! 023
      end if                                                                        ! 024
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(i==NXX_FG+1)                                                          then ! 025 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 026 
                             AG(m,N_+1)         = BC_1YZ(j,k)                       ! 027 
      end if                                                                        ! 028 
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(j==1)                                                                 then ! 029 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 030 
                             AG(m,N_+1)         = BC_X0Z(i,k)                       ! 031 
      end if                                                                        ! 032 
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(j==NYY_FG+1)                                                          then ! 033 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 034 
                             AG(m,N_+1)         = BC_X1Z(i,k)                       ! 035 
      end if                                                                        ! 036 
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
      if(k==1)                                                                 then ! 037 
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 038
                             AG(m,N_+1)         = BC_XY0(i,j)                       ! 039
      end if                                                                        ! 040
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(k==NZZ_FG+1)                                                          then ! 041
                             AG(m,M_A_P(i,j,k)) = 1.D0;               way = 1       ! 042
                             AG(m,N_+1)         = BC_XY1(i,j)                       ! 043
      end if                                                                        ! 044
      
      if(way == 0)                                                             then ! 045    
               A_W =       1.d0/h_x**2                                              ! 046
               A_E =       1.d0/h_x**2                                              ! 047
               A_S =                     1.d0/h_y**2                                ! 048
               A_N =                     1.d0/h_y**2                                ! 049
               A_U =                                   1.d0/h_z**2                  ! 050 
               A_D =                                   1.d0/h_z**2                  ! 051
               A_P =-2.d0*(1.d0/h_x**2 + 1.d0/h_y**2 + 1.d0/h_z**2)                 ! 052

      if(i-1 >= I_top)  AG(m,M_A_P(i-1,j  ,k  )) = A_W                              ! 053   
      if(i-1 >= I_top)                             A_W = 0.d0                       ! 054
	     
      if(i+1 <= I_bot)  AG(m,M_A_P(i+1,j  ,k  )) = A_E                              ! 055
      if(i+1 <= I_bot)                             A_E = 0.d0                       ! 056
	                                                                                   
      if(j-1 >= J_top)  AG(m,M_A_P(i  ,j-1,k  )) = A_S                              ! 057  
      if(j-1 >= J_top)                             A_S = 0.d0                       ! 058  
	                                                                                   
      if(j+1 <= J_bot)  AG(m,M_A_P(i  ,j+1,k  )) = A_N                              ! 059  
      if(j+1 <= J_bot)                             A_N = 0.d0                       ! 060  
	                                                                                   
      if(k-1 >= K_top)  AG(m,M_A_P(i  ,j  ,k-1)) = A_D                              ! 061  
      if(k-1 >= K_top)                             A_D = 0.d0                       ! 062  
	                                                                                   
      if(k+1 <= K_bot)  AG(m,M_A_P(i  ,j  ,k+1)) = A_U                              ! 063 
      if(k+1 <= K_bot)                             A_U = 0.d0                       ! 064  
	   
                        AG(m,M_A_P(i  ,j  ,k  )) = A_P                              ! 065
                        AG(m,N_+1)               =-F(i,j,k) &                       ! 066        
                  - A_W*U(i-1,j  ,k  ) - A_E*U(i+1,j  ,k  ) &                       ! 067
                  - A_S*U(i  ,j-1,k  ) - A_N*U(i  ,j+1,k  ) &                       ! 068
                  - A_D*U(i  ,j  ,k-1) - A_U*U(i  ,j  ,k+1)                         ! 069
      end if                                                                        ! 070
     end do                                                                         ! 071
    end do                                                                          ! 072 
   end do                                                                           ! 073  
                     call GAUSSIAN_ELIMINATION(AG,xG,N_)                            ! 074 
                                   m  = 0                                           ! 075 
   do           k  = K_top,K_bot                                                    ! 076                                    
    do        j    = J_top,J_bot                                                    ! 077
     do     i      = I_top,I_bot;  m  = m+1                                         ! 078   
          U(i,j,k) = xG(m)                                                          ! 079
     end do                                                                         ! 080
    end do                                                                          ! 081
   end do                                                                           ! 082  

 deallocate(M_A_P,AG,xG)                                                            ! 083
end subroutine Vanka_type_smoother_FG


subroutine  Convergence_Test
real*8      u_a,Error,h,D2UDx2,D2UDy2,D2UDz2,Res,h_x,h_y_h_z,x,y,z
               h_x = 1.d0/dfloat(NXX_FG)                         ! 001    
               h_y = 1.d0/dfloat(NYY_FG)                         ! 002    
               h_z = 1.d0/dfloat(NZZ_FG)                         ! 003    
            ResMAX = 0.d0                                        ! 004    
            ErrMAX = 0.d0                                        ! 005    
 do             k  = 2,NZZ_FG;  z = dfloat(k-1)/dfloat(NZZ_FG)   ! 006             
  do          j    = 2,NYY_FG;  y = dfloat(j-1)/dfloat(NYY_FG)   ! 007        
   do       i      = 2,NXX_FG;  x = dfloat(i-1)/dfloat(NXX_FG)   ! 008    
               u_a = dexp(a_*x + b_*y + c_*z)                    ! 009    
             Error = dabs(u_a - U(i,j,k))                        ! 010    
          if(Error > ErrMAX)    ErrMAX = Error                   ! 011    
                 h = U(i  ,j  ,k  )                      *2.d0   ! 012         
            D2uDx2 =(U(i-1,j  ,k  ) - h + U(i+1,j  ,k  ))/h_x**2 ! 013     
            D2uDy2 =(U(i  ,j-1,k  ) - h + U(i  ,j+1,k  ))/h_y**2 ! 014    
            D2uDz2 =(U(i  ,j  ,k-1) - h + U(i  ,j  ,k+1))/h_z**2 ! 015    
               Res = dabs(F(i,j,k) + D2uDx2 + D2uDy2 + D2uDz2)   ! 016    
            if(Res > ResMAX)    ResMAX = Res                     ! 017    
   end do                                                        ! 018    
  end do                                                         ! 019
 end do                                                          ! 020
end subroutine Convergence_Test


SUBROUTINE   GAUSSIAN_ELIMINATION(AG,xG,NG)
Implicit     Real*8(A-H,O-Z)
REAL*8,      ALLOCATABLE :: AG(:,:), xG(:)
N1         = NG+1
DO      K  = 1,NG
 K1        = K+1
 S         = AG(K,K)
 J         = K
 DO     I  = K1,NG
  R        = AG(I,K)
  IF(DABS(R).Gt.DABS(S)) then
    S      = R
    J      = I
  end if
 end do
 IF(J.Ne.K)              then
  DO    I  = K,N1
   R       = AG(K,I)
   AG(K,I) = AG(J,I)
   AG(J,I) = R
  end do
 end if
 DO     J  = K1,N1
   AG(K,J) = AG(K,J)/S
 end do
 DO   I    = K1,NG
 R         = AG(I,K)
  DO    J  = K1,N1
   AG(I,J) = AG(I,J)-AG(K,J)*R
  end do
 end do
end do
xG(NG)     = AG(NG,N1)
DO    I    = NG-1,1,-1
S          = AG(I,N1)
 DO     J  = I+1,NG
  S        = S-AG(I,J)*xG(J)
 end do
xG(I)      = S
end do
END SUBROUTINE GAUSSIAN_ELIMINATION

END MODULE Parallel_Finest_Grid_Smoothing_02


































