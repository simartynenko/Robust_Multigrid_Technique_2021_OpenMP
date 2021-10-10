!
!       Parallel Vanka-type smoother
!
!       Chapter II, Section 2
!
!       File: Parallel_Finest_Grid_Smoothing_M.f90
!
program Parallel_Finest_Grid_Smoothing
USE     OMP_LIB                                                               ! 001     
USE     Parallel_Finest_Grid_Smoothing_01                                     ! 002     
USE     Parallel_Finest_Grid_Smoothing_02                                     ! 003     
                                                                              ! 004     
open(1,file='Parallel_Finest_Grid_Smoothing.03')                              ! 005     
           write(1,1)                                                         ! 006     
           IterMAX = 2;               Nthreads=3                              ! 007
DO        Nthreads_= 1,Nthreads,Nthreads-1                                    ! 008     
 do       NColors  = 10,15                                                    ! 009
                                      NX_block=3; NY_block=3; NZ_block=3      ! 010     
            NXX_FG = NColors*Nthreads*NX_block-1                              ! 011     
            NYY_FG = NColors*Nthreads*            NY_block-1                  ! 012          
            NZZ_FG = NColors*Nthreads*                        NZ_block-1      ! 013     
           X_block = NColors*Nthreads                                         ! 014     
           Y_block = NColors*Nthreads                                         ! 015     
           Z_block = NColors*Nthreads                                         ! 016     

                a_ = 1.d0                                                     ! 017     
                b_ = 1.d0                                                     ! 018          
                c_ = 1.d0                                                     ! 019     

                     call MSSize                                              ! 020     
                     call Starting_Guess_and_Boundary_Conditions              ! 021     
         The_Time0 = omp_get_wtime()                                          ! 022     
  do        IterFG = 1,IterMAX                                                ! 023     
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       X-planes      
   do     NColors_ = 1,NColors                                                ! 024              
!$OMP PARALLEL NUM_THREADS(Nthreads_) DEFAULT(PRIVATE) &                      ! 025
!$OMP SHARED(X_block,Y_block,Z_block,NX_block,NY_block, &
!$OMP NZ_block,NColors_,NColors, &                                            ! 026
!$OMP NXX_FG,NYY_FG,NZZ_FG,U,F,BC_0YZ,BC_1YZ, &
!$OMP BC_X0Z,BC_X1Z,BC_XY0,BC_XY1)                                            ! 027
!$OMP DO                                                                      ! 028     
    do    I_block  = NColors_,X_block,NColors                                 ! 029     
                                I_top =(I_block-1)*NX_block+1                 ! 030     
                                I_bot = I_block   *NX_block                   ! 031    

     do   J_block  = 1,Y_block                                                ! 032            
                                J_top =(J_block-1)*NY_block+1                 ! 033       
                                J_bot = J_block   *NY_block                   ! 034      

      do  K_block  = 1,Z_block                                                ! 035          
                                K_top =(K_block-1)*NZ_block+1                 ! 036            
                                K_bot = K_block   *NZ_block                   ! 037       

       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)       ! 038      

      end do                                                                  ! 039   
     end do                                                                   ! 040   
    end do                                                                    ! 041   
!$OMP END DO                                                                  ! 042   
!$OMP END PARALLEL                                                            ! 043   
   end do                                                                     ! 044
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       Y-planes       
   do     NColors_ = 1,NColors                                                ! 045 
!$OMP PARALLEL NUM_THREADS(Nthreads_)  DEFAULT(PRIVATE) &                     ! 046 
!$OMP SHARED(X_block,Y_block,Z_block,NX_block,NY_block, &
!$OMP NZ_block,NColors_,NColors, &                                            ! 047
!$OMP NXX_FG,NYY_FG,NZZ_FG,U,F,BC_0YZ, &
!$OMP BC_1YZ,BC_X0Z,BC_X1Z,BC_XY0,BC_XY1)                                     ! 048
!$OMP DO                                                                      ! 049      
    do    J_block  = NColors_,Y_block,NColors                                 ! 050 
                                J_top =(J_block-1)*NY_block+1                 ! 051 
                                J_bot = J_block   *NY_block                   ! 052 

     do   I_block  = 1,X_block                                                ! 053 
                                I_top =(I_block-1)*NX_block+1                 ! 054 
                                I_bot = I_block   *NX_block                   ! 055 

      do  K_block  = 1,Z_block                                                ! 056 
                                K_top =(K_block-1)*NZ_block+1                 ! 057 
                                K_bot = K_block   *NZ_block                   ! 058 

       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)       ! 059 

      end do                                                                  ! 060 
     end do                                                                   ! 061      
    end do                                                                    ! 062
!$OMP END DO                                                                  ! 063 
!$OMP END PARALLEL                                                            ! 064 
   end do                                                                     ! 065 
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       Z-planes      
   do     NColors_ = 1,NColors                                                ! 066 
!$OMP PARALLEL NUM_THREADS(Nthreads_)  DEFAULT(PRIVATE) &                     ! 067 
!$OMP SHARED(X_block,Y_block,Z_block,NX_block,NY_block, &
!$OMP NZ_block,NColors_,NColors, &                                            ! 068
!$OMP NXX_FG,NYY_FG,NZZ_FG,U,F,BC_0YZ, &
!$OMP BC_1YZ,BC_X0Z,BC_X1Z,BC_XY0,BC_XY1)                                     ! 069
!$OMP DO                                                                      ! 070     
    do    K_block  = NColors_,Z_block,NColors                                 ! 071     
                                K_top =(K_block-1)*NZ_block+1                 ! 072         
                                K_bot = K_block   *NZ_block                   ! 073         

     do   J_block  = 1,Y_block                                                ! 074         
                                J_top =(J_block-1)*NY_block+1                 ! 075         
                                J_bot = J_block   *NY_block                   ! 076         

      do  I_block  = 1,X_block                                                ! 077         
                                I_top =(I_block-1)*NX_block+1                 ! 078     
                                I_bot = I_block   *NX_block                   ! 079     

       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)       ! 080          

      end do                                                                  ! 081     
     end do                                                                   ! 082   
    end do                                                                    ! 083   
!$OMP END DO                                                                  ! 084   
!$OMP END PARALLEL                                                            ! 085   
   end do                                                                     ! 086   
  end do                                                                      ! 087

       if(Nthreads_== 1)                                                       then     ! 088 
            TimeC1(NColors) = omp_get_wtime() - The_Time0                               ! 089 
       else                                                                             ! 090 
            TimeC2(NColors) = omp_get_wtime() - The_Time0                               ! 091 
            call Convergence_Test                                                       ! 092 
write(1,2)  Nthreads,NColors,NXX_FG-1,(NXX_FG-1)**3,TimeC1(NColors)/TimeC2(NColors),  & ! 093 
            TimeC1(NColors)/TimeC2(NColors)/dfloat(Nthreads),TimeC1(NColors),Res_000, & ! 094 
            Err_000,ResMAX,ErrMAX                                                       ! 095
       end if                                                                           ! 096
            deallocate(u,f,BC_0YZ,BC_1YZ,BC_X0Z,BC_X1Z,BC_XY0,BC_XY1)                   ! 097
 end do                                                                                 ! 098
END DO                                                                                  ! 099
close(1)                                                                                ! 100
1 format('    T    C    N       Grid       S       E              Time', &
         '         R0        E0        R         E')
2 format(i5,i5,i5,i12,3x,F5.2,3x,F5.2,'  |  ',F12.1,'  |',4(1x,D9.2))    
end program Parallel_Finest_Grid_Smoothing