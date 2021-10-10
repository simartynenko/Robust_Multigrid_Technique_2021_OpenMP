!
!       Parallel Coarse Grid Smoothing
!
!       Chapter III, Section 4
!
!       File: Parallel_Coarse_Grid_Smoothing_M.f90
!
program Parallel_Coarse_Grid_Smoothing_M
USE     Parallel_Coarse_Grid_Smoothing_1                            ! 001
USE     Parallel_Coarse_Grid_Smoothing_2                            ! 002
USE     RMT_3D_2020_OpenMP                                          ! 003
USE     RMT_3D_UFG                                                  ! 004
USE     OMP_LIB                                                     ! 005

open(1,file='Parallel_Coarse_Grid_Smoothing.03')                    ! 006
              write(1,1)                                            ! 007     
                          Nthreads=3                                ! 008
DO   MOCKBA = 1,5                                                   ! 009
DO Nthreads_= 1,Nthreads, Nthreads-1                                ! 010  

     NXX_FG = 50*MOCKBA;  N_coarsest_X = 5; LevelXd = 1; a_ = 1.d0; NX_block = 3 ! 011 
     NYY_FG = 50*MOCKBA;  N_coarsest_Y = 5; LevelYd = 1; b_ = 1.d0; NY_block = 3 ! 012
     NZZ_FG = 50*MOCKBA;  N_coarsest_Z = 5; LevelZd = 1; c_ = 1.d0; NZ_block = 3 ! 013
           
    The_number_of_smoothing_iterations = 1                          ! 014
    The_number_of_multigrid_iterations = 1                          ! 015
                                MGI_DC = 4                          ! 016
              call MSSize                                           ! 017
              call U_Finest_Grid(1,1,1)                             ! 018    
              call Coarse_Grids('+')                                ! 019
              call Starting_Guess_and_Boundary_Conditions           ! 020   
              The_Time0 = omp_get_wtime()                           ! 021
                                                                  
!$OMP PARALLEL NUM_THREADS(Nthreads_) DEFAULT(PRIVATE)            & ! 022
!$OMP SHARED(LeXmaxFG,NXX_FG,LevelXd,N_coarsest_X,                &
!$OMP NX_block,Xv_FG,Xf_FG,                                       & ! 023
!$OMP LeYmaxFG,NYY_FG,LevelYd,N_coarsest_Y,NY_block,Yv_FG,Yf_FG,  & ! 024
!$OMP LeZmaxFG,NZZ_FG,LevelZd,N_coarsest_Z,NZ_block,Zv_FG,Zf_FG,  & ! 025
!$OMP MG_v_,N_MGa_,MG_c_,N_MGi_,Istart_,               &            ! 026
!$OMP The_number_of_smoothing_iterations,              &            ! 027
!$OMP The_number_of_multigrid_iterations, MGI_DC,      &            ! 028
!$OMP DCSchedule_27,LevelM,NXYZmax,The_Time0,          &            ! 029
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,hatU,RHSF,a_,b_,c_)             ! 020
             call RMT_Iterations                                    ! 031
!$OMP END PARALLEL                                                  ! 032
                                                                        
       if(Nthreads_== 1)                                                       then   ! 033 
             TimeC1(MOCKBA) = omp_get_wtime() - The_Time0                             ! 034 
       else                                                                           ! 035 
             TimeC2(MOCKBA) = omp_get_wtime() - The_Time0                             ! 036 
             call Convergence_Test                                                    ! 037 
write(1,2)   Nthreads,NColors,NXX_FG-1,(NXX_FG-1)**3,TimeC1(MOCKBA)/TimeC2(MOCKBA), & ! 038 
             TimeC1(MOCKBA)/TimeC2(MOCKBA)/dfloat(Nthreads),TimeC1(MOCKBA),Res_000, & ! 039 
             Err_000,ResMAX,ErrMAX                                                    ! 040
       end if                                                                         ! 041
             call All_Deallocate                                                      ! 042 
             deallocate(hatU,RHSF,B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1)                      ! 043 
END DO                                                                                ! 044   
END DO                                                                                ! 045
close(1)                                                                              ! 046 
1 format('    T    C    N       Grid       S       E              Time', &
         '         R0        E0        R         E')
2 format(i5,i5,i5,i12,3x,F5.2,3x,F5.2,'  |  ',F12.1,'  |',4(1x,D9.2))    
CONTAINS

Subroutine RMT_Iterations
INTEGER,   ALLOCATABLE :: MG_v(:,:),  Poin_V(:,:,:),  N_MGa(:,:,:)                       ! 001 
INTEGER,   ALLOCATABLE :: MG_c(:,:),  Poin_C(:,:,:),  N_MGi(:,:,:),  Istart(:,:,:)       ! 002 

INTEGER,   ALLOCATABLE :: TxV_(:), TxF_(:)                                               ! 003
INTEGER,   ALLOCATABLE :: TyV_(:), TyF_(:)                                               ! 004
INTEGER,   ALLOCATABLE :: TzV_(:), TzF_(:)                                               ! 005 
INTEGER    NxxD,NyyD,NzzD,LeXmax,LeYmax,LeZmax                                           ! 006 

       j = 3**LevelM;                               ka =          -j                     ! 007 
       i =(LevelM+1)*(NXYZmax+5+2*j);               kb = NXYZmax+2+j                     ! 008  
       
           allocate(   MG_v(3,0:i       ),    MG_c(3,0:i      )                       )  ! 009 
           allocate( Poin_V(3,0:LevelM,j),  N_MGa(3,0:LevelM,j)                       )  ! 010 
           allocate( Poin_C(3,0:LevelM,j),  N_MGi(3,0:LevelM,j),  Istart(3,0:LevelM,j))  ! 011 

           allocate( TxV_(ka:kb), TxF_(ka:kb) )                                          ! 012 
           allocate( TyV_(ka:kb), TyF_(ka:kb) )                                          ! 013 
           allocate( TzV_(ka:kb), TzF_(ka:kb) )                                          ! 014
           
!$OMP DO                                                                                 ! 015
 do    m = 1,27                                                                          ! 016 
 call Dynamic_Grids_X(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                        ! 017 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                  ! 017 
 call Dynamic_Grids_Y(LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                        ! 018  
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                  ! 018 
 call Dynamic_Grids_Z(LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                        ! 019 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                  ! 019   
 call      Mini_Cycle(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                        ! 020 
                      LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                        ! 020 
                      LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                        ! 020 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                  ! 020
 end do                                                                                  ! 021 
!$OMP END DO                                                                             ! 022 

 deallocate(TxV_, TxF_, TyV_, TyF_, TzV_, TzF_)                                          ! 023 
 deallocate(MG_v, Poin_V, N_MGa, MG_c, Poin_C, N_MGi, Istart)                            ! 024 
end Subroutine RMT_Iterations

end program Parallel_Coarse_Grid_Smoothing_M
 