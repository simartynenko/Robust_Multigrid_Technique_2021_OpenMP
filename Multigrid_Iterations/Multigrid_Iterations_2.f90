MODULE  Multigrid_Iterations_2
USE     Multigrid_Iterations_1
USE     RMT_3D_2020_OpenMP
    
PRIVATE    
PUBLIC ::     Matrix_and_Vector_in_Mini_Cycle, Mini_Cycle,Convergence_Test
PUBLIC ::     Starting_Guess_and_Boundary_Conditions,Vanka_type_smoother_FG

CONTAINS 
    
subroutine Starting_Guess_and_Boundary_Conditions
real*8     h
allocate(    cu(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2)                                          )
allocate(  hatU(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2),  RHSF(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2) )
allocate(  B_X0(           0:NYY_FG+2,0:NZZ_FG+2),  B_X1(           0:NYY_FG+2,0:NZZ_FG+2) )
allocate(  B_Y0(0:NXX_FG+2,           0:NZZ_FG+2),  B_Y1(0:NXX_FG+2,           0:NZZ_FG+2) )
allocate(  B_Z0(0:NXX_FG+2,0:NYY_FG+2           ),  B_Z1(0:NXX_FG+2,0:NYY_FG+2           ) )
                                                  
               cu = 0.d0;   hatU = 0.d0;   RHSF = 0.d0               ! 001
             B_X0 = 0.d0;   B_X1 = 0.d0                              ! 002
             B_Y0 = 0.d0;   B_Y1 = 0.d0                              ! 003
             B_Z0 = 0.d0;   B_Z1 = 0.d0                              ! 004
! 
! Right Hand Side Function
!            
                h  = a_**2+b_**2+c_**2                               ! 005
          Res_000  = dexp(a_+b_+c_)*h                                ! 006
          Err_000  = dexp(a_+b_+c_)                                  ! 007
   do           k  = 1,NZZ_FG+1                                      ! 008
    do        j    = 1,NYY_FG+1                                      ! 009
     do     i      = 1,NXX_FG+1                                      ! 010
       RHSF(i,j,k) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))*h ! 011
     end do                                                          ! 012
    end do                                                           ! 013
   end do                                                            ! 014
!
! Boundary Conditions
!
    do          k  = 1,NZZ_FG+1                                      ! 015
     do       j    = 1,NYY_FG+1                                      ! 016
         B_X0(j,k) = dexp(              b_*Yv_FG(j) + c_*Zv_FG(k))   ! 017
         B_X1(j,k) = dexp(a_          + b_*Yv_FG(j) + c_*Zv_FG(k))   ! 018
     end do                                                          ! 019
    end do                                                           ! 020
    
   do           k  = 1,NZZ_FG+1                                      ! 021
     do       i    = 1,NXX_FG+1                                      ! 022
         B_Y0(i,k) = dexp(a_*Xv_FG(i)               + c_*Zv_FG(k))   ! 023
         B_Y1(i,k) = dexp(a_*Xv_FG(i) + b_          + c_*Zv_FG(k))   ! 024
     end do                                                          ! 025
   end do                                                            ! 026
   
   do           j  = 1,NYY_FG+1                                      ! 027
    do        i    = 1,NXX_FG+1                                      ! 028
         B_Z0(i,j) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j)              )   ! 029
         B_Z1(i,j) = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_         )   ! 030
    end do                                                           ! 031
   end do                                                            ! 032
end subroutine Starting_Guess_and_Boundary_Conditions


subroutine    Mini_Cycle(LeXmax,GridsXd,NxxD,TxV_,TxF_, &
                         LeYmax,GridsYd,NyyD,TyV_,TyF_, &
                         LeZmax,GridsZd,NzzD,TzV_,TzF_, &
                         MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)
                                                                                          
INTEGER,      ALLOCATABLE ::  TxV(:),       TxF(:),       TxV_(:),        TxF_(:)      ! 001
INTEGER,      ALLOCATABLE ::  TyV(:),       TyF(:),       TyV_(:),        TyF_(:)      ! 002
INTEGER,      ALLOCATABLE ::  TzV(:),       TzF(:),       TzV_(:),        TzF_(:)      ! 003
INTEGER,      ALLOCATABLE :: MG_v(:,:),  Poin_V(:,:,:),  N_MGa(:,:,:)                  ! 004
INTEGER,      ALLOCATABLE :: MG_c(:,:),  Poin_C(:,:,:),  N_MGi(:,:,:),  Istart(:,:,:)  ! 005
REAL*8 ,      ALLOCATABLE :: BCX0(:,:),    BCX1(:,:),       J3(:,:,:)                  ! 006
REAL*8 ,      ALLOCATABLE :: BCY0(:,:),    BCY1(:,:),      J3H(:,:,:)                  ! 007
REAL*8 ,      ALLOCATABLE :: BCZ0(:,:),    BCZ1(:,:)                                   ! 008
REAL*8 ,      ALLOCATABLE :: BUX0(:,:),    BUX1(:,:)                                   ! 009
REAL*8 ,      ALLOCATABLE :: BUY0(:,:),    BUY1(:,:)                                   ! 010
REAL*8 ,      ALLOCATABLE :: BUZ0(:,:),    BUZ1(:,:)                                   ! 011

real*8        h,h_x,h_y,h_z,D2uDx2,D2uDy2,D2uDz2                                       ! 012
integer       LevelX,LevelY,LevelZ,LevelCC                                             ! 013
integer       GridsXd,GridsYd,GridsZd,GridsX,GridsY,GridsZ                             ! 014
integer       Nx,NxA,topX,botX,Ny,NyA,topY,botY,Nz,NzA,topZ,botZ                       ! 015
                                                                                           
              II =(3**LevelXd-1)/2                                                     ! 016
              JJ =(3**LevelYd-1)/2                                                     ! 017
              KK =(3**LevelZd-1)/2                                                     ! 018
    allocate( J3(-II:NxxD+2+II,-JJ:NyyD+2+JJ,-KK:NzzD+2+KK) )                          ! 019
                                                                                           
              ka =          -3**LevelM                                                 ! 020
              kb = NXYZmax+2+3**LevelM                                                 ! 021
    allocate( TxV(ka:kb), TxF(ka:kb) )                                                 ! 022
    allocate( TyV(ka:kb), TyF(ka:kb) )                                                 ! 023
    allocate( TzV(ka:kb), TzF(ka:kb) )                                                 ! 024
                                                                                           
    allocate( BCX0(           0:NYY_FG+2,0:NZZ_FG+2),  BCX1(           0:NYY_FG+2,0:NZZ_FG+2) )
    allocate( BCY0(0:NXX_FG+2,           0:NZZ_FG+2),  BCY1(0:NXX_FG+2,           0:NZZ_FG+2) )
    allocate( BCZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BCZ1(0:NXX_FG+2,0:NYY_FG+2           ) )
    allocate( BUX0(           0:NYY_FG+2,0:NZZ_FG+2),  BUX1(           0:NYY_FG+2,0:NZZ_FG+2) )
    allocate( BUY0(0:NXX_FG+2,           0:NZZ_FG+2),  BUY1(0:NXX_FG+2,           0:NZZ_FG+2) )
    allocate( BUZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BUZ1(0:NXX_FG+2,0:NYY_FG+2           ) )
    allocate(  J3H(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2)                                          )
    
do             MGI = 1,MGI_DC                                                          ! 032
                                                                                       
              BCX0 = 0.d0; BCX1 = 0.d0;   BUX0 = B_X0; BUX1 = B_X1                     ! 033
              BCY0 = 0.d0; BCY1 = 0.d0;   BUY0 = B_Y0; BUY1 = B_Y1                     ! 034
              BCZ0 = 0.d0; BCZ1 = 0.d0;   BUZ0 = B_Z0; BUZ1 = B_Z1                     ! 035
                                                                                       
 if(MGIteration==1)                                                            then    ! 036
     if(MGI==1)                                                                then    
              BCX0 = B_X0; BCX1 = B_X1;   BUX0 = 0.d0; BUX1 = 0.d0                     ! 037
              BCY0 = B_Y0; BCY1 = B_Y1;   BUY0 = 0.d0; BUY1 = 0.d0                     ! 038
              BCZ0 = B_Z0; BCZ1 = B_Z1;   BUZ0 = 0.d0; BUZ1 = 0.d0                     ! 039
     end if                                                                            
               J3H = RHSF                                                              ! 040
 else                                                                                  
               J3H = 0.d0                                                              ! 041
              h_x  = 1.d0/dfloat(NXX_FG)                                               ! 042
              h_y  = 1.d0/dfloat(NYY_FG)                                               ! 043
              h_z  = 1.d0/dfloat(NZZ_FG)                                               ! 044
     do         k  = 2,NZZ_FG                                                          ! 045
      do      j    = 2,NYY_FG                                                          ! 046
       do   i      = 2,NXX_FG                                                          ! 047
                 h = hatU(i  ,j  ,k  )                         *2.d0                   ! 048
            D2uDx2 =(hatU(i-1,j  ,k  ) - h + hatU(i+1,j  ,k  ))/h_x**2                 ! 049
            D2uDy2 =(hatU(i  ,j-1,k  ) - h + hatU(i  ,j+1,k  ))/h_y**2                 ! 050
            D2uDz2 =(hatU(i  ,j  ,k-1) - h + hatU(i  ,j  ,k+1))/h_z**2                 ! 051
        J3H(i,j,k) = RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2                            ! 052
       end do                                                                          ! 053
      end do                                                                           ! 054
     end do                                                                            ! 055
 end if                                                                                ! 056
                                                                                       
 do         LevelCC= max(LeXmax,LeYmax,LeZmax        ),0,-1                            ! 057
            LevelX = min(LeXmax              ,LevelCC)                                 ! 058
            LevelY = min(       LeYmax       ,LevelCC)                                 ! 059
            LevelZ = min(              LeZmax,LevelCC)                                  
                                                                                       ! 060
       call Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_,        &    ! 061
                                            NyyD,LevelY,TyV,TyF,TyV_,TyF_,        &    ! 062
                                            NzzD,LevelZ,TzV,TzF,TzV_,TzF_,        &    ! 063
                                            BUX0,BUX1,BUY0,BUY1,BUZ0,BUZ1,J3,J3H, MGI) ! 064
  do        GridsX = 1,3**LevelX                                                       ! 065
       call Index_Mapping_X('-',LevelX,GridsX,Nx,NxA,topX,botX,TxV,TxF,TxV_,TxF_, &    ! 066
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)      ! 067
   do       GridsY = 1,3**LevelY                                                       ! 068
       call Index_Mapping_Y('-',LevelY,GridsY,Ny,NyA,topY,botY,TyV,TyF,TyV_,TyF_, &    ! 069
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)      ! 070
    do      GridsZ = 1,3**LevelZ                                                       ! 071
       call Index_Mapping_Z('-',LevelZ,GridsZ,Nz,NzA,topZ,botZ,TzV,TzF,TzV_,TzF_, &    ! 072
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)      ! 073
                                                                                           
       call Vanka_type_Smoother(LevelX,LeXmax,TxV,TxF,TxV_,TxF_,Nx,               &    ! 074
                                LevelY,LeYmax,TyV,TyF,TyV_,TyF_,Ny,               &    ! 075
                                LevelZ,LeZmax,TzV,TzF,TzV_,TzF_,Nz,               &    ! 076
                                BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3)                      ! 077
    end do                                                                             ! 078
   end do                                                                              ! 079
  end do                                                                               ! 080
 end do                                                                                ! 081
 if(MGIteration  ==  1)                                                        then    ! 082
  do            k  = 1,NzzD+1;                           KK  = TzV_(k)                 ! 083
   do        j     = 1,NyyD+1;                        JJ     = TyV_(j)                 ! 084
    do    i        = 1,NxxD+1;                     II        = TxV_(i)                 ! 085
    hatU(II,JJ,KK) = hatU(II,JJ,KK)           + cu(II,JJ,KK)                           ! 086
                                                cu(II,JJ,KK) = 0.d0                    ! 087
    end do                                                                             ! 088
   end do                                                                              ! 089
  end do                                                                               ! 090
 end if                                                                                ! 091
end do                                                                                 ! 092
deallocate(BCX0, BCX1, BCY0, BCY1, BCZ0, BCZ1, BUX0, BUX1, BUY0, BUY1, BUZ0, BUZ1, J3, J3H)
end subroutine  Mini_Cycle                                                                        


subroutine    Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_, &
                                              NyyD,LevelY,TyV,TyF,TyV_,TyF_, &
                                              NzzD,LevelZ,TzV,TzF,TzV_,TzF_, &
                                              BUX0,BUX1,BUY0,BUY1,BUZ0,BUZ1,J3,J3H, MGI) 

INTEGER,      ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)      ! 001
INTEGER,      ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)      ! 002
INTEGER,      ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:)      ! 003

REAL*8 ,      ALLOCATABLE :: BUX0(:,:), BUX1(:,:),   F3X(:),     J3H(:,:,:)  ! 004
REAL*8 ,      ALLOCATABLE :: BUY0(:,:), BUY1(:,:),   F3Y(:),      J3(:,:,:)  ! 005
REAL*8 ,      ALLOCATABLE :: BUZ0(:,:), BUZ1(:,:),   F3Z(:)                  ! 006

integer       IV_m,IV_0,IV_p,JV_m,JV_0,JV_p,KV_m,KV_0,KV_p, MGI                ! 007
real*8        AX_m,AX_0,AX_p,Sx,AY_m,AY_0,AY_p,Sy,AZ_m,AZ_0,AZ_p,Sz,xi,h___    ! 008
real*8        h_x,h_y,h_z,h_x_,h_y_,h_z_,h_Lx,h_Ly,h_Lz,d2udx2,d2udy2,d2udz2,h ! 009

              ka   =          -3**LevelM                                     ! 010
              kb   = NXYZmax+2+3**LevelM                                     ! 011
    allocate( F3X(ka:kb), F3Y(ka:kb), F3Z(ka:kb) )                           ! 012
!
!- - - - -    X-Integral        - - - - -
!
               Ih  =(3**LevelXd-1)/2                                         ! 013
                N_ = NxxD                                                    ! 014
 if(LevelXd==0) N_ = NXX_FG                                                  ! 015
                ia =    1-Ih;  if(TxV_(ia)+Ih <        1) ia = ia + 1        ! 016
                ib = N_+1+Ih;  if(TxV_(ib)-Ih > NXX_FG+1) ib = ib - 1        ! 017
            F3X    = 0.d0                                                    ! 018
           do   i  = ia,ib                                                   ! 019
            F3X(i) = Xf_FG(min(NXX_FG+1,TxF_(i))) - Xf_FG(max(0,TxF_(i-1)))  ! 020
           end do                                                            ! 021
!                                                                           
!- - - - -    Y-Integral        - - - - -
!
               Jh  =(3**LevelYd-1)/2                                         ! 022
                N_ = NyyD                                                    ! 023
 if(LevelYd==0) N_ = NYY_FG                                                  ! 024
                ja =    1-Jh;  if(TyV_(ja)+Jh <        1) ja = ja + 1        ! 025
                jb = N_+1+Jh;  if(TyV_(jb)-Jh > NYY_FG+1) jb = jb - 1        ! 026
            F3Y    = 0.d0                                                    ! 027
           do   j  = ja,jb                                                   ! 028
            F3Y(j) = Yf_FG(min(NYY_FG+1,TyF_(j))) - Yf_FG(max(0,TyF_(j-1)))  ! 029
           end do                                                            ! 030
!                                                                                 
!- - - - -    Z-Integral        - - - - -
!
               Kh  =(3**LevelZd-1)/2                                         ! 031
                N_ = NzzD                                                    ! 032
 if(LevelZd==0) N_ = NZZ_FG                                                  ! 033
                ka =    1-Kh;  if(TzV_(ka)+Kh <        1) ka = ka + 1        ! 034
                kb = N_+1+Kh;  if(TzV_(kb)-Kh > NZZ_FG+1) kb = kb - 1        ! 035
            F3Z    = 0.d0                                                    ! 036
           do   k  = ka,kb                                                   ! 037
            F3Z(k) = Zf_FG(min(NZZ_FG+1,TzF_(k))) - Zf_FG(max(0,TzF_(k-1)))  ! 038
           end do                                                            ! 039

                J3 = 0.d0                                                                 ! 040
      do        k  = 1,NzzD+1                                                             ! 041
       do     j    = 1,NyyD+1                                                             ! 042
        do  i      = 1,NxxD+1                                                             ! 043
                                                                                          ! 044
      do        k_ = TzV_(k)-Kh,TzV_(k)+Kh                                                ! 045
      if((1.le.k_).and.(k_.le.NZZ_FG+1))                                       then       ! 046
                                              h_z_ = Zf_FG(k_) - Zf_FG(k_-1)              ! 047
       do     j_   = TyV_(j)-Jh,TyV_(j)+Jh                                                ! 048
       if((1.le.j_).and.(j_.le.NYY_FG+1))                                      then       ! 049
                                              h_y_ = Yf_FG(j_) - Yf_FG(j_-1)              ! 050
        do  i_     = TxV_(i)-Ih,TxV_(i)+Ih                                                ! 051
        if((1.le.i_).and.(i_.le.NXX_FG+1))                                     then       ! 052
                                              h_x_ = Xf_FG(i_) - Xf_FG(i_-1)              ! 053
                                                                                          ! 054
         J3(i,j,k)= J3(i,j,k) + h_x_*h_y_*h_z_*J3H(i_,j_,k_)                              ! 055
                                                                                          ! 056
        end if                                                                            ! 057
        end do                                                                            ! 058
       end if                                                                             ! 059
       end do                                                                             ! 060
      end if                                                                              ! 061
      end do                                                                              ! 062
                                                                                          ! 063
     end do                                                                               ! 064
    end do                                                                                ! 065
   end do                                                                                 ! 066
    if(MGIteration== 1)                                                        then       ! 067
              h_x  = 1.d0/dfloat(NXX_FG)
              h_y  = 1.d0/dfloat(NYY_FG)
              h_z  = 1.d0/dfloat(NZZ_FG)
              h_Lx = 1.d0/(h_x*dfloat(3**LevelXd))**2                                     ! 068
              h_Ly = 1.d0/(h_y*dfloat(3**LevelYd))**2                                     ! 069
              h_Lz = 1.d0/(h_z*dfloat(3**LevelZd))**2                                     ! 070
              h___ = h_x*dfloat(3**LevelXd)*h_y*dfloat(3**LevelYd)*h_z*dfloat(3**LevelZd) ! 071
      do        k  = 1,NzzD+1; KV_m = TzV_(k-1); KV_0 = TzV_(k); KV_p = TzV_(k+1)         ! 072
       do     j    = 1,NyyD+1; JV_m = TyV_(j-1); JV_0 = TyV_(j); JV_p = TyV_(j+1)         ! 073
        do  i      = 1,NxxD+1; IV_m = TxV_(i-1); IV_0 = TxV_(i); IV_p = TxV_(i+1)         ! 074
                                                                                            
        if((IV_0==       1).or.(JV_0==       1).or.(KV_0==       1).or. &                 ! 075
           (IV_0==NXX_FG+1).or.(JV_0==NYY_FG+1).or.(KV_0==NZZ_FG+1)) then                 ! 076
                                                                                            
        J3(i,j,k)  = 0.d0                                                                 ! 077
                                                                                            
        else                                                                              ! 078
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  X-Direction                           
         if( i==1)                                                   then                 ! 079
                xi =      Xv_FG(IV_0) *dsqrt(h_Lx)                                        ! 080
              AX_m = 0.d0                                                                 ! 081
              AX_o = 2.d0          /xi                                                    ! 082
              AX_p = 2.d0/(xi+1.d0)                                                       ! 083
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 084
                   *             BUX0(     JV_0,KV_0)                                     ! 085
            d2udx2 =(      -AX_o*hatU(IV_0,JV_0,KV_0)          &                          ! 086
                           +AX_p*hatU(IV_p,JV_0,KV_0) ) * h_Lx + Sx                       ! 087
         end if                                                                           ! 088
         if((1<i).and.(i<NxxD+1))                                    then                 ! 089
            d2udx2 = 1.d0*(      hatU(IV_m,JV_0,KV_0)          &                          ! 090
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 091
                           +     hatU(IV_p,JV_0,KV_0) ) * h_Lx                            ! 092
         end if                                                                           ! 093
         if(           i==NxxD+1)                                    then                 ! 094
                xi =(1.d0-Xv_FG(IV_0))*dsqrt(h_Lx)                                        ! 095
              AX_p = 0.d0                                                                 ! 096
              AX_o = 2.d0          /xi                                                    ! 097
              AX_m = 2.d0/(xi+1.d0)                                                       ! 098
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 099
                   *             BUX1(     JV_0,KV_0)                                     ! 100
            d2udx2 =(       AX_m*hatU(IV_m,JV_0,KV_0)          &                          ! 101
                           -AX_o*hatU(IV_0,JV_0,KV_0) ) * h_Lx + Sx                       ! 102
         end if                                                                           ! 103
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  Y-Direction                         
                                                                                          ! 104
         if( j==1)                                                   then                 ! 105
                xi =      Yv_FG(JV_0) *dsqrt(h_Ly)                                        ! 106
              AY_m = 0.d0                                                                 ! 107
              AY_o = 2.d0          /xi                                                    ! 108
              AY_p = 2.d0/(xi+1.d0)                                                       ! 109
              Sy   = 2.d0/(xi+1.d0)/xi                  * h_Ly &                          ! 110
                   *             BUY0(IV_0,     KV_0)                                     ! 111
            d2udy2 =(      -AY_o*hatU(IV_0,JV_0,KV_0)          &                          ! 112
                           +AY_p*hatU(IV_0,JV_p,KV_0) ) * h_Ly + Sy                       ! 113
         end if                                                                           ! 114
         if((1<j).and.(j<NyyD+1))                                    then                 ! 115
            d2udy2 =(            hatU(IV_0,JV_m,KV_0)          &                          ! 116
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 117
                           +     hatU(IV_0,JV_p,KV_0) ) * h_Ly                            ! 118
         end if                                                                           ! 119
         if(           j==NyyD+1)                                    then                 ! 120
                xi =(1.d0-Yv_FG(JV_0))*dsqrt(h_Ly)                                        ! 121
              AY_p = 0.d0                                                                 ! 122
              AY_o = 2.d0          /xi                                                    ! 123
              AY_m = 2.d0/(xi+1.d0)                                                       ! 124
              Sy   = 2.d0/(xi+1.d0)/xi                  * h_Ly &                          ! 125
                   *             BUY1(IV_0,     KV_0)                                     ! 126
            d2udy2 =(       AY_m*hatU(IV_0,JV_m,KV_0)          &                          ! 127
                           -AY_o*hatU(IV_0,JV_0,KV_0) ) * h_Ly + Sy                       ! 128
         end if                                                                           ! 129
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  Z-Direction
         if( k==1)                                                   then                 ! 130
                xi =      Zv_FG(KV_0) *dsqrt(h_Lz)                                        ! 131
              AZ_m = 0.d0                                                                 ! 132
              AZ_o = 2.d0          /xi                                                    ! 133
              AZ_p = 2.d0/(xi+1.d0)                                                       ! 134
              Sz   = 2.d0/(xi+1.d0)/xi                  * h_Lz &                          ! 135
                   *             BUZ0(IV_0,JV_0     )                                     ! 136
            d2udz2 =(      -AZ_o*hatU(IV_0,JV_0,KV_0)          &                          ! 137
                           +AZ_p*hatU(IV_0,JV_0,KV_p) ) * h_Lz + Sz                       ! 138
         end if                                                                           ! 139
         if((1<k).and.(k<NzzD+1))                                    then                 ! 140
            d2udz2 =(            hatU(IV_0,JV_0,KV_m)          &                          ! 141
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 142
                           +     hatU(IV_0,JV_0,KV_p) ) * h_Lz                            ! 143
         end if                                                                           ! 144
         if(           k==NzzD+1)                                    then                 ! 145
                xi =(1.d0-Zv_FG(KV_0))*dsqrt(h_Lz)                                        ! 146
              AZ_p = 0.d0                                                                 ! 147
              AZ_o = 2.d0          /xi                                                    ! 148
              AZ_m = 2.d0/(xi+1.d0)                                                       ! 149
              Sz   = 2.d0/(xi+1.d0)/xi                  * h_Lz &                          ! 150
                   *             BUZ1(IV_0,JV_0     )                                     ! 151
            d2udz2 =(       AZ_m*hatU(IV_0,JV_0,KV_m)          &                          ! 152
                           -AZ_o*hatU(IV_0,JV_0,KV_0) ) * h_Lz + Sz                       ! 153
         end if                                                                           ! 154
        J3(i,j,k)  = J3(i,j,k) - (d2udx2 + d2udy2 + d2udz2)*h___                          ! 155
        end if                                                                            ! 156
        end do                                                                            ! 157
       end do                                                                             ! 158
      end do                                                                              ! 159
    end if                                                                                ! 160
    call RHSI3D(LevelX,LevelY,LevelZ,NxxD,NyyD,NzzD,F3X,F3Y,F3Z,J3,ia,ib,ja,jb,ka,kb)     ! 161
    deallocate(F3X,F3Y,F3Z)                                                               ! 162
end subroutine Matrix_and_Vector_in_Mini_Cycle  

                                              
subroutine  Vanka_type_Smoother(LevelX,LeXmax,TxV,TxF,TxV_,TxF_,Nx, &       
                                LevelY,LeYmax,TyV,TyF,TyV_,TyF_,Ny, &       
                                LevelZ,LeZmax,TzV,TzF,TzV_,TzF_,Nz, &             
                                BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3)            
implicit    real*8 (a-h,o-z)                                                 
integer     X_block,Y_block,Z_block                                                   ! 001
real*8      h_x,h_y,h_z,h_Lx,h_Ly,h_Lz                                                ! 002
REAL*8 ,    ALLOCATABLE :: BCX0(:,:), BCX1(:,:)                                       ! 003
REAL*8 ,    ALLOCATABLE :: BCY0(:,:), BCY1(:,:)                                       ! 004
REAL*8 ,    ALLOCATABLE :: BCZ0(:,:), BCZ1(:,:)                                       ! 005
REAL*8 ,    ALLOCATABLE ::   J3(:,:,:), AG(:,:),    xG(:)                             ! 006
INTEGER,    ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)                 ! 007
INTEGER,    ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)                 ! 008
INTEGER,    ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:),   M_A_P(:,:,:) ! 009
        h_x   = 1.d0/dfloat(NXX_FG)                                                   ! 010
        h_y   = 1.d0/dfloat(NYY_FG)                                                   ! 011
        h_z   = 1.d0/dfloat(NZZ_FG)                                                   ! 012
        h_Lx  = 1.d0/(h_x*dfloat(3**(LevelX + LevelXd)))**2                           ! 013
        h_Ly  = 1.d0/(h_y*dfloat(3**(LevelY + LevelYd)))**2                           ! 014
        h_Lz  = 1.d0/(h_z*dfloat(3**(LevelZ + LevelZd)))**2                           ! 015

    NX_block_ =                   NX_block                                            ! 016
    NY_block_ =          NY_block                                                     ! 017
    NZ_block_ = NZ_block                                                              ! 018
    NSIL      = The_number_of_smoothing_iterations                                    ! 019
    NSIL_     = The_number_of_smoothing_iterations                                    ! 020
    
            i = NX_block_*NY_block_*NZ_block_;  allocate(AG(i,i+1),xG(i))             ! 021
    
 if(NX_block_ > Nx+1) NX_block_ = Nx+1                                                ! 022
 if(NY_block_ > Ny+1) NY_block_ = Ny+1                                                ! 023
 if(NZ_block_ > Nz+1) NZ_block_ = Nz+1                                                ! 024

if(((Nx+1)/NX_block_)*NX_block_== Nx+1)                                        then   ! 025
                       X_block  =(Nx+1)/NX_block_                                     ! 026
else                                                                                  ! 027
                       X_block  =(Nx+1)/NX_block_ + 1                                 ! 028
end if                                                                                ! 029
if(((Ny+1)/NY_block_)*NY_block_== Ny+1)                                        then   ! 030
                       Y_block  =(Ny+1)/NY_block_                                     ! 031
else                                                                                  ! 032
                       Y_block  =(Ny+1)/NY_block_ + 1                                 ! 033
end if                                                                                ! 034
if(((Nz+1)/NZ_block_)*NZ_block_== Nz+1)                                        then   ! 035
                       Z_block  =(Nz+1)/NZ_block_                                     ! 036
else                                                                                  ! 037
                       Z_block  =(Nz+1)/NZ_block_ + 1                                 ! 038
end if                                                                                ! 039

if((levelX==LeXmax).and.(levelY==LeYmax).and.(levelZ==LeZmax))       NSIL_ = 2*NSIL   ! 040
if(X_block*Y_block*Z_block==1)                                                 then
  do I_block  = 1,X_block                                                             ! 040
                           I_top =(I_block-1)*NX_block_+1                             ! 041
                           I_bot = I_block   *NX_block_                               ! 042
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 043
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 044

  do J_block  = 1,Y_block                                                             ! 045
                           J_top =(J_block-1)*NY_block_+1                             ! 046
                           J_bot = J_block   *NY_block_                               ! 047
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 048
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 049

  do K_block  = 1,Z_block                                                             ! 050
                           K_top =(K_block-1)*NZ_block_+1                             ! 051
                           K_bot = K_block   *NZ_block_                               ! 052
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 053
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 054
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 055
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 056
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 057
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3,AG,xG)                    ! 058
  end do                                                                              ! 059
  end do                                                                              ! 060
  end do                                                                              ! 061
else                                                                                  ! 062
 do Iteration = 1,NSIL_                                                               ! 063
  do I_block  = 1,X_block                                                             ! 064
                           I_top =(I_block-1)*NX_block_+1                             ! 065
                           I_bot = I_block   *NX_block_                               ! 066
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 067
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 068
                                                                                      ! 069
  do J_block  = 1,Y_block                                                             ! 070
                           J_top =(J_block-1)*NY_block_+1                             ! 071
                           J_bot = J_block   *NY_block_                               ! 072
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 073
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 074
                                                                                      ! 075
  do K_block  = 1,Z_block                                                                  
                           K_top =(K_block-1)*NZ_block_+1                             ! 076
                           K_bot = K_block   *NZ_block_                               ! 077
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 078
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 079
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 080
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 081
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 082
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3,AG,xG)                    ! 083
  end do                                                                              ! 084
  end do                                                                              ! 085
  end do                                                                              ! 086
                                                                                        
  do K_block  = 1,Z_block                                                             ! 087
                           K_top =(K_block-1)*NZ_block_+1                             ! 088
                           K_bot = K_block   *NZ_block_                               ! 089
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 090
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 091
                                                                                      
  do I_block  = 1,X_block                                                             ! 092
                           I_top =(I_block-1)*NX_block_+1                             ! 093
                           I_bot = I_block   *NX_block_                               ! 094
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 095
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 096
                                                                                        
  do J_block  = 1,Y_block                                                             ! 097
                           J_top =(J_block-1)*NY_block_+1                             ! 098
                           J_bot = J_block   *NY_block_                               ! 099
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 100
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 101
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 102
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 103
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 104
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3,AG,xG)                    ! 105
  end do                                                                              ! 106
  end do                                                                              ! 107
  end do                                                                              ! 108
                                                                                           
  do J_block  = 1,Y_block                                                             ! 109
                           J_top =(J_block-1)*NY_block_+1                             ! 110
                           J_bot = J_block   *NY_block_                               ! 111
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 112
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 113
                                                                                        
  do K_block  = 1,Z_block                                                             ! 114
                           K_top =(K_block-1)*NZ_block_+1                             ! 115
                           K_bot = K_block   *NZ_block_                               ! 116
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 117
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 118
                                                                                        
  do I_block  = 1,X_block                                                             ! 119
                           I_top =(I_block-1)*NX_block_+1                             ! 120
                           I_bot = I_block   *NX_block_                               ! 121
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 122
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 123
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 124
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 125
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 126
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3,AG,xG)                    ! 127
  end do                                                                              ! 128
  end do                                                                              ! 129
  end do                                                                              ! 130

 end do                                                                               ! 131
end if                                                                                ! 132
deallocate(AG,xG)                                                                     ! 133
end subroutine  Vanka_type_Smoother                                                        
                                
                                
SUBROUTINE  Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx,  &
                            TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly,  &     
                            TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz,  &
                            BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,J3,AG,xG)
integer     M_A_P_,way                                                                    ! 001
real*8      h_x,h_y,h_z,h_Lx,h_Ly,h_Lz,h                                                  ! 002
real*8      H_W,H_E,H_S,H_N,H_U,H_D,xi                                                    ! 003
real*8      A_W,A_E,A_S,A_N,A_U,A_D,A_Px,A_Py,A_Pz,S_x,S_y,S_z                            ! 004
REAL*8 ,    ALLOCATABLE ::   J3(:,:,:), AG(:,:),    xG(:)                                 ! 005
INTEGER,    ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)                     ! 006
INTEGER,    ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)                     ! 007
INTEGER,    ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:),   M_A_P(:,:,:)     ! 008
REAL*8 ,    ALLOCATABLE :: BCX0(:,:), BCX1(:,:)                                           ! 009
REAL*8 ,    ALLOCATABLE :: BCY0(:,:), BCY1(:,:)                                           ! 010
REAL*8 ,    ALLOCATABLE :: BCZ0(:,:), BCZ1(:,:)                                           ! 011
                                                                                               
       AG = 0.D0; xG = 0.D0;  allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))        ! 012
                                                                                               
        m = 0                                                                             ! 013
   do   k = K_top,K_bot                                                                   ! 014
    do  j = J_top,J_bot                                                                   ! 015
     do i = I_top,I_bot                                                                   ! 016
        m = m+1                                                                           ! 017
                                       M_A_P(i,j,k) = m                                   ! 018
     end do                                                                               ! 019
    end do                                                                                ! 020
   end do                                                                                 ! 021
       N_ = m                                                                             ! 022
                                                                                               
        m = 0                                                                             ! 023
   do   k = K_top,K_bot;  KK = TzV_(TzV(k));    KKK = TzV(k)                              ! 024
    do  j = J_top,J_bot;  JJ = TyV_(TyV(j));    JJJ = TyV(j)                              ! 025
     do i = I_top,I_bot;  II = TxV_(TxV(i));    III = TxV(i)                              ! 026
        m = m+1                                                                           ! 027
      way = 0;                               M_A_P_ = M_A_P(i,j,k)                        ! 028
! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((i==   1).and.(II==       1))                                         then       ! 029
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 030
                                       AG(m,N_+1)   = BCX0(   JJ,KK)                      ! 031
      end if                                                                              ! 032
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((i==Nx+1).and.(II==NXX_FG+1))                                         then       ! 033
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 034
                                       AG(m,N_+1)   = BCX1(   JJ,KK)                      ! 035
      end if                                                                              ! 036
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((j==   1).and.(JJ==       1))                                         then       ! 037
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 038
                                       AG(m,N_+1)   = BCY0(II,   KK)                      ! 039
      end if                                                                              ! 040
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((j==Ny+1).and.(JJ==NYY_FG+1))                                         then       ! 041
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 042
                                       AG(m,N_+1)   = BCY1(II,   KK)                      ! 043
      end if                                                                              ! 044
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((k==   1).and.(KK==       1))                                         then       ! 045
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 046
                                       AG(m,N_+1)   = BCZ0(II,JJ   )                      ! 047
      end if                                                                              ! 048
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -            
      if((k==Nz+1).and.(KK==NZZ_FG+1))                                         then       ! 049
                                       AG(m,M_A_P_) = 1.D0;           way = 1             ! 050
                                       AG(m,N_+1)   = BCZ1(II,JJ   )                      ! 051
      end if                                                                              ! 052

      if(way == 0)                                                             then       ! 053
     H_W = 0.d0; H_E = 0.d0; S_x = 0.d0                                                   ! 054
   if(i == 1)                                                                  then       ! 055
     xi  =          Xv_FG(II)             *dsqrt(h_Lx)                                    ! 056
    A_W  =  0.d0                                                                          ! 057
    A_Px =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lx                                     ! 058
    A_E  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lx; H_E = A_E*cu(TxV_(TxV(i+1)),JJ,KK) ! 059
    S_x  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lx      *   BCX0(               JJ,KK) ! 060
   end if                                                                                 ! 061
   if((1<i).and.(i<Nx+1))                                                      then       ! 062
    A_W  =                                       h_Lx; H_W = A_W*cu(TxV_(TxV(i-1)),JJ,KK) ! 063
    A_E  =                                       h_Lx; H_E = A_E*cu(TxV_(TxV(i+1)),JJ,KK) ! 064
    A_Px =-A_W-A_E                                                                        ! 065
   end if                                                                                 ! 066
   if(i == Nx+1)                                                               then       ! 067
      xi = (1.d0 -  Xv_FG(II))            *dsqrt(h_Lx)                                    ! 068
    A_W  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lx; H_W = A_W*cu(TxV_(TxV(i-1)),JJ,KK) ! 069
    A_Px =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lx                                     ! 070
    A_E  =  0.d0                                                                          ! 071
    S_x  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lx      *   BCX1(               JJ,KK) ! 072
   end if                                                                                 ! 073
     H_S = 0.d0; H_N = 0.d0; S_y = 0.d0                                                   ! 074
   if(j == 1)                                                                  then       ! 075
     xi  =          Yv_FG(JJ)             *dsqrt(h_Ly)                                    ! 076
    A_S  =  0.d0                                                                          ! 077
    A_Py =-(1.d0 - (xi-1.d0)          /xi)* 2.d0*h_Ly                                     ! 078
    A_N  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Ly; H_N = A_N*cu(II,TyV_(TyV(j+1)),KK) ! 079
    S_y  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Ly      *   BCY0(II,               KK) ! 080
   end if                                                                                 ! 081
   if((1<j).and.(j<Ny+1))                                                      then       ! 082
    A_S  =                                       h_Ly; H_S = A_S*cu(II,TyV_(TyV(j-1)),KK) ! 083
    A_N  =                                       h_Ly; H_N = A_N*cu(II,TyV_(TyV(j+1)),KK) ! 084
    A_Py =-A_S-A_N                                                                        ! 085
   end if                                                                                 ! 086
   if(j == Ny+1)                                                               then       ! 087
     xi  = (1.d0 -  Yv_FG(JJ))            *dsqrt(h_Ly)                                    ! 088
    A_S  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Ly; H_S = A_S*cu(II,TyV_(TyV(j-1)),KK) ! 089
    A_Py =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Ly                                     ! 090
    A_N  =  0.d0                                                                          ! 091
    S_y  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Ly      *   BCY1(II,               KK) ! 092
   end if                                                                                 ! 093
     H_D = 0.d0; H_U = 0.d0; S_z = 0.d0                                                   ! 094
   if(k == 1)                                                                  then       ! 095
     xi  =          Zv_FG(KK)             *dsqrt(h_Lz)                                    ! 096
    A_D  =  0.d0                                                                          ! 097
    A_Pz =-(1.d0 - (xi-1.d0)          /xi)* 2.d0*h_Lz                                     ! 098
    A_U  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lz; H_U = A_U*cu(II,JJ,TzV_(TzV(k+1))) ! 099
    S_z  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lz      *   BCZ0(II,JJ               ) ! 100
   end if                                                                                 ! 101
   if((1<k).and.(k<Nz+1))                                                      then       ! 102
    A_D  =                                       h_Lz; H_D = A_D*cu(II,JJ,TzV_(TzV(k-1))) ! 103
    A_U  =                                       h_Lz; H_U = A_U*cu(II,JJ,TzV_(TzV(k+1))) ! 104
    A_Pz =-A_D-A_U                                                                        ! 105
   end if                                                                                 ! 106
   if(k == Nz+1)                                                               then       ! 107
     xi  = (1.d0 -  Zv_FG(KK))            *dsqrt(h_Lz)                                    ! 108
    A_D  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lz; H_D = A_D*cu(II,JJ,TzV_(TzV(k-1))) ! 109
    A_Pz =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lz                                     ! 110
    A_U  =  0.d0                                                                          ! 111
    S_z  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lz      *   BCZ1(II,JJ               ) ! 112
   end if                                                                                 ! 113
                                                                                               
       if(i-1 >= I_top)     AG(m,M_A_P(i-1,j  ,k  )) = A_W                                ! 114
       if(i-1 >= I_top)                                H_W = 0.D0                         ! 115
                                                                                               
       if(i+1 <= I_bot)     AG(m,M_A_P(i+1,j  ,k  )) = A_E                                ! 116
       if(i+1 <= I_bot)                                H_E = 0.D0                         ! 117
                                                                                               
       if(j-1 >= J_top)     AG(m,M_A_P(i  ,j-1,k  )) = A_S                                ! 118
       if(j-1 >= J_top)                                H_S = 0.D0                         ! 119
                                                                                               
       if(j+1 <= J_bot)     AG(m,M_A_P(i  ,j+1,k  )) = A_N                                ! 120
       if(j+1 <= J_bot)                                H_N = 0.D0                         ! 121
                                                                                               
       if(k-1 >= K_top)     AG(m,M_A_P(i  ,j  ,k-1)) = A_D                                ! 122
       if(k-1 >= K_top)                                H_D = 0.D0                         ! 123
                                                                                               
       if(k+1 <= K_bot)     AG(m,M_A_P(i  ,j  ,k+1)) = A_U                                ! 124
       if(k+1 <= K_bot)                                H_U = 0.D0                         ! 125
                            AG(m,M_A_P(i  ,j  ,k  )) = A_Px + A_Py + A_Pz                 ! 126
                            AG(m,N_+1)               = J3(III,JJJ,KKK) &                  ! 127
                       - S_x - S_y - S_z - H_W - H_E - H_S - H_N - H_D - H_U              ! 128
      end if                                                                              ! 129
     end do                                                                               ! 130
    end do                                                                                ! 131
   end do                                                                                 ! 132
                            call GAUSSIAN_ELIMINATION(AG,xG,N_)                           ! 133
        m = 0                                                                             ! 134
   do   k = K_top,K_bot                                                                   ! 135
    do  j = J_top,J_bot                                                                   ! 136
     do i = I_top,I_bot                                                                   ! 137
        m = m+1                                                                           ! 138
                            cu(TxV_(TxV(i)),TyV_(TyV(j)),TzV_(TzV(k))) = xG(m)            ! 139
     end do                                                                               ! 140
    end do                                                                                ! 141
   end do                                                                                 ! 142
                            deallocate(M_A_P)                                             ! 143
END SUBROUTINE   Vanka_iteration                                                                   
                                                                                           
                                                                                          
subroutine Vanka_type_smoother_FG(I_top,I_bot,J_top,J_bot,K_top,K_bot)                    
INTEGER    way,N_,                I_top,I_bot,J_top,J_bot,K_top,K_bot               ! 001            
REAL*8     A_W,A_E,A_S,A_N,A_D,A_U,A_P,h_x,h_y,h_z                                  ! 002            
INTEGER,   ALLOCATABLE ::  M_A_P(:,:,:)                                             ! 003            
REAL*8 ,   ALLOCATABLE ::     AG(:,:),   xG(:)                                      ! 004        
                                                                                                    
                N_ = NX_block*NY_block*NZ_block                                     ! 005           
 allocate(AG(N_,N_+1), xG(N_))                                                      ! 006           
 allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))                               ! 007
                AG = 0.D0;    xG = 0.D0                                             ! 008
               h_x = 1.d0/dfloat(NXX_FG)                                            ! 009
               h_y = 1.d0/dfloat(NYY_FG)                                            ! 010
               h_z = 1.d0/dfloat(NZZ_FG)                                            ! 011
                                                                                         
                                   m  = 0                                           ! 012
   do           k  = K_top,K_bot                                                    ! 013
    do        j    = J_top,J_bot                                                    ! 014
     do     i      = I_top,I_bot;  m  = m+1                                         ! 015
      M_A_P(i,j,k) =               m                                                ! 016
     end do                                                                         ! 018 
    end do                                                                          ! 019 
   end do                                                                           ! 020 
                                   N_ = m                                           ! 021 
                                   m  = 0                                           ! 022 
   do           k  = K_top,K_bot                                                    ! 023 
    do        j    = J_top,J_bot                                                    ! 024 
     do     i      = I_top,I_bot;  m  = m+1                                         ! 025 
               A_W = 0.D0;        A_E = 0.D0;  A_S = 0.D0;  A_N = 0.D0              ! 026 
               A_D = 0.D0;        A_U = 0.D0;  A_P = 0.D0;            way = 0       ! 027 
! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 028 
      if(i==1)                                                                 then ! 029 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 030 
                                AG(m,N_+1)         = B_X0(j,k)                      ! 031 
      end if                                                                        ! 032 
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 033 
      if(i==NXX_FG+1)                                                          then ! 034 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 035 
                                AG(m,N_+1)         = B_X1(j,k)                      ! 036 
      end if                                                                        ! 037 
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 038 
      if(j==1)                                                                 then ! 039 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 040 
                                AG(m,N_+1)         = B_Y0(i,k)                      ! 041 
      end if                                                                        ! 042 
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 043 
      if(j==NYY_FG+1)                                                          then ! 044 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 045 
                                AG(m,N_+1)         = B_Y1(i,k)                      ! 046 
      end if                                                                        ! 047 
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 048 
      if(k==1)                                                                 then ! 049 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 050 
                                AG(m,N_+1)         = B_Z0(i,j)                      ! 051 
      end if                                                                        ! 052 
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - ! 053 
      if(k==NZZ_FG+1)                                                          then ! 054 
                                AG(m,M_A_P(i,j,k)) = 1.D0;            way = 1       ! 055 
                                AG(m,N_+1)         = B_Z1(i,j)                      ! 056 
      end if                                                                        ! 057
                                                                                     
      if(way == 0)                                                             then ! 058    
          A_W =       1.d0/h_x**2                                                   ! 059    
          A_E =       1.d0/h_x**2                                                   ! 060    
          A_S =                     1.d0/h_y**2                                     ! 061    
          A_N =                     1.d0/h_y**2                                     ! 062    
          A_U =                                   1.d0/h_z**2                       ! 063    
          A_D =                                   1.d0/h_z**2                       ! 064    
          A_P =-2.d0*(1.d0/h_x**2 + 1.d0/h_y**2 + 1.d0/h_z**2)                      ! 065
                                                                                      
      if(i-1 >= I_top)  AG(m,M_A_P(i-1,j  ,k  )) = A_W                              ! 066
      if(i-1 >= I_top)                             A_W = 0.d0                       ! 067
	                                                                            
      if(i+1 <= I_bot)  AG(m,M_A_P(i+1,j  ,k  )) = A_E                              ! 068
      if(i+1 <= I_bot)                             A_E = 0.d0                       ! 069
	                                                                            
      if(j-1 >= J_top)  AG(m,M_A_P(i  ,j-1,k  )) = A_S                              ! 070
      if(j-1 >= J_top)                             A_S = 0.d0                       ! 071
	                                                                               
      if(j+1 <= J_bot)  AG(m,M_A_P(i  ,j+1,k  )) = A_N                              ! 072
      if(j+1 <= J_bot)                             A_N = 0.d0                       ! 073
	                                                                            
      if(k-1 >= K_top)  AG(m,M_A_P(i  ,j  ,k-1)) = A_D                              ! 074
      if(k-1 >= K_top)                             A_D = 0.d0                       ! 075
	                                                                            
      if(k+1 <= K_bot)  AG(m,M_A_P(i  ,j  ,k+1)) = A_U                              ! 076
      if(k+1 <= K_bot)                             A_U = 0.d0                       ! 077
	                                                                              
                        AG(m,M_A_P(i  ,j  ,k  )) = A_P                              ! 078      
                        AG(m,N_+1)               =     RHSF(i  ,j  ,k  ) &          ! 079      
                        - A_W*hatU(i-1,j  ,k  )  - A_E*hatU(i+1,j  ,k  ) &          ! 080      
                        - A_S*hatU(i  ,j-1,k  )  - A_N*hatU(i  ,j+1,k  ) &          ! 081      
                        - A_D*hatU(i  ,j  ,k-1)  - A_U*hatU(i  ,j  ,k+1)            ! 082 
      end if                                                                        ! 083      
     end do                                                                         ! 084      
    end do                                                                          ! 085      
   end do                                                                           ! 086      
                     call GAUSSIAN_ELIMINATION(AG,xG,N_)                            ! 087  
                                   m  = 0                                           ! 096                      
   do           k  = K_top,K_bot                                                    ! 088      
    do        j    = J_top,J_bot                                                    ! 089      
     do     i      = I_top,I_bot;  m = m+1                                          ! 090      
       hatU(i,j,k) = xG(m)                                                          ! 091      
     end do                                                                         ! 092      
    end do                                                                          ! 093      
   end do                                                                           ! 094      
 deallocate(M_A_P,AG,xG)                                                            ! 095
end subroutine Vanka_type_smoother_FG                                              


subroutine  Convergence_Test                                                
real*8      u_a,Error,h,D2UDx2,D2UDy2,D2UDz2,Res,h_x,h_y_h_z           ! 001
               h_x = 1.d0/dfloat(NXX_FG)                               ! 002
               h_y = 1.d0/dfloat(NYY_FG)                               ! 003
               h_z = 1.d0/dfloat(NZZ_FG)                               ! 004
            ResMAX = 0.d0                                              ! 005
            ErrMAX = 0.d0                                              ! 006
 do             k  = 2,NZZ_FG                                          ! 007
  do          j    = 2,NYY_FG                                          ! 008
   do       i      = 2,NXX_FG                                          ! 009
               u_a = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))     ! 010
             Error = dabs(u_a - hatU(i,j,k))                           ! 011
          if(Error > ErrMAX)    ErrMAX = Error                         ! 012
                 h = hatU(i  ,j  ,k  )                         *2.d0   ! 013
            D2uDx2 =(hatU(i-1,j  ,k  ) - h + hatU(i+1,j  ,k  ))/h_x**2 ! 014
            D2uDy2 =(hatU(i  ,j-1,k  ) - h + hatU(i  ,j+1,k  ))/h_y**2 ! 015
            D2uDz2 =(hatU(i  ,j  ,k-1) - h + hatU(i  ,j  ,k+1))/h_z**2 ! 016
               Res = dabs(RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2)      ! 017
            if(Res > ResMAX)    ResMAX = Res                           ! 018
   end do                                                              ! 019
  end do                                                               ! 020
 end do                                                                ! 021
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

END MODULE Multigrid_Iterations_2