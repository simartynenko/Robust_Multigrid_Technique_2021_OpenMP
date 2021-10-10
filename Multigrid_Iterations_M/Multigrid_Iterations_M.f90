!
!       Parallel RMT
!
!       Chapter IV
!
!       File: Multigrid_Iterations_M.f90    
!
program Multigrid_Iterations_M
USE     Multigrid_Iterations_1            ! 001
USE     Multigrid_Iterations_2            ! 002
USE     RMT_3D_2020_OpenMP                ! 003
USE     RMT_3D_UFG                        ! 004
USE     OMP_LIB                           ! 005

OPEN(1, file='Multigrid_iterations.03')   ! 006
        write(1,1)                        ! 007
        IterMAX  = 1;         Nthreads=3  ! 008
DO      NColors  = 10,15                  ! 009
DO      Nthreads_= 1,Nthreads,Nthreads-1
                              NX_block=3; NY_block=3; NZ_block=3  ! 011
    NXX_FG = NColors*Nthreads*NX_block-1                          ! 012
    NYY_FG = NColors*Nthreads*            NY_block-1              ! 013
    NZZ_FG = NColors*Nthreads*                        NZ_block-1  ! 014
   I_block_= NColors*Nthreads                                     ! 015
   J_block_= NColors*Nthreads                                     ! 016
   K_block_= NColors*Nthreads                                     ! 017
           
             N_coarsest_X = 5; LevelXd = 1; a_ = 1.d0             ! 018
             N_coarsest_Y = 5; LevelYd = 1; b_ = 1.d0             ! 019
             N_coarsest_Z = 5; LevelZd = 1; c_ = 1.d0             ! 020
                              
    The_number_of_smoothing_iterations = 1                        ! 021
    The_number_of_multigrid_iterations = 4                        ! 022
                                MGI_DC = 4                        ! 023
                                
         if(((NXX_FG+1)/NX_block)*NX_block== NXX_FG+1)                         then  ! 024
                                   I_block_=(NXX_FG+1)/NX_block                      ! 025
         else                                                                        ! 026
                                   I_block_=(NXX_FG+1)/NX_block + 1                  ! 027
         end if                                                                      ! 028
         if(((NYY_FG+1)/NY_block)*NY_block== NYY_FG+1)                         then  ! 029
                                   J_block_=(NYY_FG+1)/NY_block                      ! 030    
         else                                                                        ! 031    
                                   J_block_=(NYY_FG+1)/NY_block + 1                  ! 032    
         end if                                                                      ! 033    
         if(((NZZ_FG+1)/NZ_block)*NZ_block== NZZ_FG+1)                         then  ! 034    
                                   K_block_=(NZZ_FG+1)/NZ_block                      ! 035    
         else                                                                        ! 036    
                                   K_block_=(NZZ_FG+1)/NZ_block + 1                  ! 037    
         end if                                                                      ! 038    
             call MSSize                                                             ! 039    
             call U_Finest_Grid(1,1,1)                                               ! 040    
             call Coarse_Grids('+')                                                  ! 041    
             call Starting_Guess_and_Boundary_Conditions                             ! 042    
             The_Time0 = omp_get_wtime()                                             ! 043    

DO MGIteration = 1,The_number_of_multigrid_iterations                                ! 044         
!$OMP PARALLEL NUM_THREADS(Nthreads_) DEFAULT(PRIVATE)            &       ! 045
!$OMP SHARED(LeXmaxFG,NXX_FG,LevelXd,N_coarsest_X,NX_block,Xv_FG,Xf_FG, &       ! 046
!$OMP LeYmaxFG,NYY_FG,LevelYd,N_coarsest_Y,NY_block,Yv_FG,Yf_FG, &       ! 047
!$OMP LeZmaxFG,NZZ_FG,LevelZd,N_coarsest_Z,NZ_block,Zv_FG,Zf_FG, &       ! 048
!$OMP MG_v_,N_MGa_,MG_c_,N_MGi_,Istart_,                         &       ! 049
!$OMP The_number_of_smoothing_iterations, MGIteration,           &       ! 050
!$OMP The_number_of_multigrid_iterations, MGI_DC,                &       ! 051
!$OMP DCSchedule_27,LevelM,NXYZmax,The_Time0,                    &       ! 052
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,cu,hatU,RHSF,a_,b_,c_)               ! 053
             call RMT_Iterations                                          ! 054 
!$OMP END PARALLEL                                                        ! 055 

            if(MGIteration > 1)  hatU = hatU + cu; cu = 0.d0              ! 056
             
  do        IterFG = 1,IterMAX                                            ! 057
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       X-planes      
   do     NColors_ = 1,NColors                                            ! 058
!$OMP PARALLEL NUM_THREADS(Nthreads_)  DEFAULT(PRIVATE)           &       ! 059
!$OMP SHARED(I_block_,J_block_,K_block_,NX_block,NY_block,NZ_block,     &       ! 060
!$OMP NColors_,NColors, NXX_FG,NYY_FG,NZZ_FG,                    &       ! 061
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,cu,hatU,RHSF,a_,b_,c_)               ! 062
!$OMP DO                                                                  ! 063
    do    I_block  = NColors_,I_block_,NColors                            ! 064
                                I_top =(I_block-1)*NX_block+1             ! 065
                                I_bot = I_block   *NX_block               ! 068

     do   J_block  = 1,J_block_                                           ! 069
                                J_top =(J_block-1)*NY_block+1             ! 070
                                J_bot = J_block   *NY_block               ! 071

      do  K_block  = 1,K_block_                                           ! 072
                                K_top =(K_block-1)*NZ_block+1             ! 073
                                K_bot = K_block   *NZ_block               ! 074
                                                                             
       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)   ! 075

      end do                                                              ! 076
     end do                                                               ! 077
    end do                                                                ! 078
!$OMP END DO                                                              ! 079
!$OMP END PARALLEL                                                        ! 080
   end do                                                                 ! 081
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       Y-planes      
   do     NColors_ = 1,NColors                                            ! 082
!$OMP PARALLEL NUM_THREADS(Nthreads_)  DEFAULT(PRIVATE)           &       ! 083
!$OMP SHARED(I_block_,J_block_,K_block_,NX_block,NY_block,NZ_block,     &       ! 084
!$OMP NColors_,NColors, NXX_FG,NYY_FG,NZZ_FG,                    &       ! 085
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,cu,hatU,RHSF,a_,b_,c_)               ! 086
!$OMP DO                                                                  ! 087
    do    J_block  = NColors_,J_block_,NColors                            ! 088
                                J_top =(J_block-1)*NY_block+1             ! 089
                                J_bot = J_block   *NY_block               ! 090
                                                                          
     do   I_block  = 1,I_block_                                           ! 091
                                I_top =(I_block-1)*NX_block+1             ! 092
                                I_bot = I_block   *NX_block               ! 093
                                                                            
      do  K_block  = 1,K_block_                                           ! 094
                                K_top =(K_block-1)*NZ_block+1             ! 095
                                K_bot = K_block   *NZ_block               ! 096
                                                                               
       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)   ! 097 
                                                                          
      end do                                                              ! 099 
     end do                                                               ! 100 
    end do                                                                ! 101 
!$OMP END DO                                                              ! 102 
!$OMP END PARALLEL                                                        ! 103 
   end do                                                                 ! 104 
!= = = = = = = = = = = = = = = = = = = = = = = = = = = = =       Z-planes ! 105 
   do     NColors_ = 1,NColors                                            ! 106 
!$OMP PARALLEL NUM_THREADS(Nthreads_)  DEFAULT(PRIVATE)             &     ! 107
!$OMP SHARED(I_block_,J_block_,K_block_,NX_block,NY_block,NZ_block, &     ! 108
!$OMP NColors_,NColors, NXX_FG,NYY_FG,NZZ_FG,                       &     ! 109
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,cu,hatU,RHSF,a_,b_,c_)                ! 110
!$OMP DO                                                                  ! 111 
    do    K_block  = NColors_,K_block_,NColors                            ! 112 
                                K_top =(K_block-1)*NZ_block+1             ! 113 
                                K_bot = K_block   *NZ_block               ! 114
                                                                          
     do   J_block  = 1,J_block_                                           ! 115
                                J_top =(J_block-1)*NY_block+1             ! 116
                                J_bot = J_block   *NY_block               ! 117

      do  I_block  = 1,I_block_                                           ! 118
                                I_top =(I_block-1)*NX_block+1             ! 119
                                I_bot = I_block   *NX_block               ! 120
                                                                          
       call Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)   ! 121
                                                                              
      end do                                                              ! 122
     end do                                                               ! 123
    end do                                                                ! 124
!$OMP END DO                                                              ! 125
!$OMP END PARALLEL                                                        ! 126
   end do                                                                 ! 127
  end do                                                                  ! 128
END DO                                                                    ! 129
if(Nthreads_== 1)                                                              then ! 130    
     TimeC1          = omp_get_wtime() - The_Time0                                  ! 131
else                                                                                ! 132  
     TimeC2          = omp_get_wtime() - The_Time0                                  ! 133  
     call Convergence_Test                                                          ! 134
     rho =(ResMAX/Res_000)**(1./float(The_number_of_multigrid_iterations))          
write(1,2) Nthreads,NColors,NXX_FG-1,(NXX_FG-1)**3,TimeC1/TimeC2, &                 ! 135  
     TimeC1/TimeC2/dfloat(Nthreads),TimeC1,1.d0,Err_000,ResMAX/Res_000,ErrMAX,rho   ! 136  
end if                                                                              ! 137  
     call All_Deallocate                                                            ! 138  
     deallocate(cu,hatU,RHSF,B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1)                         ! 139
END DO                                                                              ! 140
END DO                                                                              ! 141
CLOSE(1)                                                                            ! 142
1 format('    T    C    N       Grid       S       E              Time', &
         '         R0        E0        R         E         AC')
2        format(i5,i5,i5,i12,3x,F5.2,3x,F5.2,'  |  ',F12.1,'  |',4(1x,D9.2),4x,F5.3)    
CONTAINS

Subroutine RMT_Iterations 
INTEGER,   ALLOCATABLE :: MG_v(:,:),  Poin_V(:,:,:),  N_MGa(:,:,:)                      ! 001
INTEGER,   ALLOCATABLE :: MG_c(:,:),  Poin_C(:,:,:),  N_MGi(:,:,:),  Istart(:,:,:)      ! 002

INTEGER,   ALLOCATABLE :: TxV_(:), TxF_(:)                                              ! 003
INTEGER,   ALLOCATABLE :: TyV_(:), TyF_(:)                                              ! 004
INTEGER,   ALLOCATABLE :: TzV_(:), TzF_(:)                                              ! 005
INTEGER    NxxD,NyyD,NzzD,LeXmax,LeYmax,LeZmax                                          ! 006

       j = 3**LevelM;                               ka =          -j                    ! 007
       i =(LevelM+1)*(NXYZmax+5+2*j);               kb = NXYZmax+2+j                    ! 008
       
           allocate(   MG_v(3,0:i       ),   MG_c(3,0:i      )                        ) ! 009
           allocate( Poin_V(3,0:LevelM,j),  N_MGa(3,0:LevelM,j)                       ) ! 010
           allocate( Poin_C(3,0:LevelM,j),  N_MGi(3,0:LevelM,j),  Istart(3,0:LevelM,j)) ! 011

           allocate( TxV_(ka:kb), TxF_(ka:kb) )                                         ! 012
           allocate( TyV_(ka:kb), TyF_(ka:kb) )                                         ! 013
           allocate( TzV_(ka:kb), TzF_(ka:kb) )                                         ! 014
           
if(MGIteration > 1)  MGI_DC = 1                                                         ! 015
           
!$OMP DO                                                                                ! 016
 do    m = 1,27                                                                         ! 017
 call Dynamic_Grids_X(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                       ! 018
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 019
 call Dynamic_Grids_Y(LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                       ! 020
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 021
 call Dynamic_Grids_Z(LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                       ! 022
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 023
 call      Mini_Cycle(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                       ! 024
                      LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                       ! 025
                      LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                       ! 026
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 027
 end do                                                                                 ! 028
!$OMP END DO                                                                            ! 029

 deallocate(TxV_, TxF_, TyV_, TyF_, TzV_, TzF_)                                         ! 030
 deallocate(MG_v, Poin_V, N_MGa, MG_c, Poin_C, N_MGi, Istart)                           ! 031
end Subroutine RMT_Iterations 

end program Multigrid_Iterations_M 