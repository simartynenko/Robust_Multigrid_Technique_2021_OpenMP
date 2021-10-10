MODULE  Parallel_Coarse_Grid_Smoothing_1
integer NX_block,NY_block,NZ_block,MGI_DC                   ! 001
integer The_number_of_smoothing_iterations                  ! 002
integer The_number_of_multigrid_iterations                  ! 003
real*8  The_Time0, TimeC1(100), TimeC2(100)                 ! 004      
real*8  a_,b_,c_,Res_000,Err_000,ResMAX,ErrMAX              ! 005

real*8, allocatable, save :: hatU(:,:,:),  RHSF(:,:,:)      ! 006
real*8, allocatable, save :: B_X0(:,:),    B_X1(:,:)        ! 007 
real*8, allocatable, save :: B_Y0(:,:),    B_Y1(:,:)        ! 008
real*8, allocatable, save :: B_Z0(:,:),    B_Z1(:,:)        ! 009 

END MODULE Parallel_Coarse_Grid_Smoothing_1                                                       