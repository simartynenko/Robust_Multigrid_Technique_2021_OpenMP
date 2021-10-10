!
! All multigrid iterations
!
! File: Integrals_M.f90    
!
program Integrals_M
USE     Integrals_1         ! 001
USE     Integrals_2         ! 002
USE     RMT_3D_2020_OpenMP  ! 003
USE     RMT_3D_UFG          ! 004
USE     OMP_LIB             ! 005

        Nthreads_= 3        ! 006
   
        NXX_FG = 67; N_coarsest_X = 5; LevelXd = 1                  ! 007
        NYY_FG = 67; N_coarsest_Y = 5; LevelYd = 1                  ! 008
        NZZ_FG = 67; N_coarsest_Z = 5; LevelZd = 1                  ! 009
                                
           call MSSize                                              ! 010
           call U_Finest_Grid(1,1,1)                                ! 011    
           call Coarse_Grids('+')                                   ! 012
           call Starting_Guess_and_Boundary_Conditions              ! 013   
             
!$OMP PARALLEL NUM_THREADS(Nthreads_) DEFAULT(PRIVATE)            & ! 014
!$OMP SHARED(LeXmaxFG,NXX_FG,LevelXd,N_coarsest_X,NX_block, &
!$OMP Xv_FG,Xf_FG, &                                                ! 015
!$OMP LeYmaxFG,NYY_FG,LevelYd,N_coarsest_Y,NY_block,Yv_FG,Yf_FG, & ! 016
!$OMP LeZmaxFG,NZZ_FG,LevelZd,N_coarsest_Z,NZ_block,Zv_FG,Zf_FG, & ! 017
!$OMP MG_v_,N_MGa_,MG_c_,N_MGi_,Istart_,              &            ! 018
!$OMP The_number_of_smoothing_iterations, MGIteration,&            ! 019
!$OMP The_number_of_multigrid_iterations, MGI_DC,     &            ! 020
!$OMP DCSchedule_27,LevelM,NXYZmax,The_Time0,         &            ! 021
!$OMP B_X0,B_X1,B_Y0,B_Y1,B_Z0,B_Z1,hatU,RHSF,a_,b_,c_)            ! 022
           call RMT_Iterations                                      ! 023
!$OMP END PARALLEL                                                  ! 024
           call All_Deallocate                                      ! 025  
           deallocate(RHSF)                                         ! 026 
CONTAINS

Subroutine RMT_Iterations 
INTEGER,   ALLOCATABLE :: MG_v(:,:),  Poin_V(:,:,:),  N_MGa(:,:,:)                      ! 001 
INTEGER,   ALLOCATABLE :: MG_c(:,:),  Poin_C(:,:,:),  N_MGi(:,:,:),  Istart(:,:,:)      ! 002 

INTEGER,   ALLOCATABLE :: TxV_(:), TxF_(:)                                              ! 003
INTEGER,   ALLOCATABLE :: TyV_(:), TyF_(:)                                              ! 004
INTEGER,   ALLOCATABLE :: TzV_(:), TzF_(:)                                              ! 005 
INTEGER    NxxD,NyyD,NzzD,LeXmax,LeYmax,LeZmax,Mmax                                     ! 006 

       j = 3**LevelM;                               ka =          -j                    ! 007 
       i =(LevelM+1)*(NXYZmax+5+2*j);               kb = NXYZmax+2+j                    ! 008  
       
           allocate(   MG_v(3,0:i       ),    MG_c(3,0:i      )                       ) ! 009 
           allocate( Poin_V(3,0:LevelM,j),  N_MGa(3,0:LevelM,j)                       ) ! 010 
           allocate( Poin_C(3,0:LevelM,j),  N_MGi(3,0:LevelM,j),  Istart(3,0:LevelM,j)) ! 011 

           allocate( TxV_(ka:kb), TxF_(ka:kb) )                                         ! 012 
           allocate( TyV_(ka:kb), TyF_(ka:kb) )                                         ! 013 
           allocate( TzV_(ka:kb), TzF_(ka:kb) )                                         ! 014
           
                               Mmax = 27                                                ! 015              
if(LevelXd*LevelYd*LevelZd==0) Mmax = 1                                                 ! 016           
 
!$OMP DO                                                                                ! 017           
 do    m = 1,Mmax                                                                       ! 018 
 call Dynamic_Grids_X(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                       ! 019 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 019 
 call Dynamic_Grids_Y(LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                       ! 020  
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 020 
 call Dynamic_Grids_Z(LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                       ! 021 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 021   
 call      Mini_Cycle(LeXmax,DCSchedule_27(m,1),NxxD,TxV_,TxF_, &                       ! 022 
                      LeYmax,DCSchedule_27(m,2),NyyD,TyV_,TyF_, &                       ! 022 
                      LeZmax,DCSchedule_27(m,3),NzzD,TzV_,TzF_, &                       ! 022 
                      MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)                 ! 022 
 end do                                                                                 ! 023 
!$OMP END DO                                                                            ! 024 
 deallocate(TxV_, TxF_, TyV_, TyF_, TzV_, TzF_)                                         ! 025 
 deallocate(MG_v, Poin_V, N_MGa, MG_c, Poin_C, N_MGi, Istart)                           ! 026 
end Subroutine RMT_Iterations 

end program Integrals_M
 