MODULE  Parallel_Finest_Grid_Smoothing_01
INTEGER NColors_,NColors,NXX_FG,NYY_FG,NZZ_FG              ! 001      
INTEGER X_block,Y_block,Z_block,NX_block,NY_block,NZ_block ! 002      
REAL*8  a_,b_,c_,Res_000,Err_000,ResMAX,ErrMAX             ! 003      
REAL*8  The_Time0, TimeC1(100), TimeC2(100)                ! 004      

REAL*8, ALLOCATABLE ::      U(:,:,:),      F(:,:,:)        ! 005      
REAL*8, ALLOCATABLE :: BC_0YZ(:,:)  , BC_1YZ(:,:)          ! 006      
REAL*8, ALLOCATABLE :: BC_X0Z(:,:)  , BC_X1Z(:,:)          ! 007      
REAL*8, ALLOCATABLE :: BC_XY0(:,:)  , BC_XY1(:,:)          ! 008      

END MODULE Parallel_Finest_Grid_Smoothing_01