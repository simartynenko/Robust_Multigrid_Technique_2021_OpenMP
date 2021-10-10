MODULE  Multigrid_Iterations_1
integer The_number_of_smoothing_iterations                                ! 001
integer The_number_of_multigrid_iterations, MGI_DC                        ! 002
INTEGER NX_block,NY_block,NZ_block,NColors_,NColors, MGIteration          ! 003
real*8  a_,b_,c_,Res_000,Err_000,ResMAX,ErrMAX, The_Time0, TimeC1, TimeC2 ! 004

real*8, allocatable, save :: hatU(:,:,:),  RHSF(:,:,:),  cu(:,:,:)        ! 005
real*8, allocatable, save :: B_X0(:,:),    B_X1(:,:)                      ! 006
real*8, allocatable, save :: B_Y0(:,:),    B_Y1(:,:)                      ! 007
real*8, allocatable, save :: B_Z0(:,:),    B_Z1(:,:)                      ! 008

END MODULE Multigrid_Iterations_1                                              