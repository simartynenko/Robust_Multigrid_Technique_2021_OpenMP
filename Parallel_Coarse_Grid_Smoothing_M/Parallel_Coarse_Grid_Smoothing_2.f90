MODULE  Parallel_Coarse_Grid_Smoothing_2 
USE     Parallel_Coarse_Grid_Smoothing_1 
USE     RMT_3D_2020_OpenMP
    
PRIVATE    
PUBLIC ::     Starting_Guess_and_Boundary_Conditions
PUBLIC ::     Matrix_and_Vector_in_Mini_Cycle, Mini_Cycle,Convergence_Test

CONTAINS 
    
subroutine Starting_Guess_and_Boundary_Conditions
real*8     h
allocate(  hatU(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2),  RHSF(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2) ) 
allocate(  B_X0(           0:NYY_FG+2,0:NZZ_FG+2),  B_X1(           0:NYY_FG+2,0:NZZ_FG+2) ) 
allocate(  B_Y0(0:NXX_FG+2,           0:NZZ_FG+2),  B_Y1(0:NXX_FG+2,           0:NZZ_FG+2) ) 
allocate(  B_Z0(0:NXX_FG+2,0:NYY_FG+2           ),  B_Z1(0:NXX_FG+2,0:NYY_FG+2           ) ) 
                                                  
              hatU = 0.d0;   RHSF = 0.d0                             ! 001                            
              B_X0 = 0.d0;   B_X1 = 0.d0                             ! 002
              B_Y0 = 0.d0;   B_Y1 = 0.d0                             ! 003
              B_Z0 = 0.d0;   B_Z1 = 0.d0                             ! 004
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


subroutine    Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_,      &
                                              NyyD,LevelY,TyV,TyF,TyV_,TyF_,      &
                                              NzzD,LevelZ,TzV,TzF,TzV_,TzF_, J3 , &
                                              BUX0, BUX1, BUY0, BUY1, BUZ0, BUZ1) 

REAL*8,       ALLOCATABLE ::  F3X(:),    F3Y(:),     F3Z(:),      J3(:,:,:)           ! 001
INTEGER,      ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)               ! 002
INTEGER,      ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)               ! 003
INTEGER,      ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:)               ! 004

REAL*8 ,      ALLOCATABLE :: BUX0(:,:), BUX1(:,:)                                     ! 005
REAL*8 ,      ALLOCATABLE :: BUY0(:,:), BUY1(:,:)                                     ! 006
REAL*8 ,      ALLOCATABLE :: BUZ0(:,:), BUZ1(:,:)                                     ! 007

integer       IV_m,IV_0,IV_p,JV_m,JV_0,JV_p,KV_m,KV_0,KV_p                            ! 008
real*8        AX_m,AX_0,AX_p,Sx,AY_m,AY_0,AY_p,Sy,AZ_m,AZ_0,AZ_p,Sz                   ! 009
real*8        h_x_,h_y_,h_z_,h_Lx,h_Ly,h_Lz,d2udx2,d2udy2,d2udz2,xi,h___              ! 010

              ka   =          -3**LevelM                                              ! 011   
              kb   = NXYZmax+2+3**LevelM                                              ! 012   
    allocate( F3X(ka:kb), F3Y(ka:kb), F3Z(ka:kb) )                                    ! 013   
!
!- - - - -    X-Integral        - - - - -
!
               Ih  =(3**LevelXd-1)/2                                                  ! 014   
                N_ = NxxD                                                             ! 015
 if(LevelXd==0) N_ = NXX_FG                                                           ! 016
                ia =    1-Ih;  if(TxV_(ia)+Ih <        1) ia = ia + 1                 ! 017
                ib = N_+1+Ih;  if(TxV_(ib)-Ih > NXX_FG+1) ib = ib - 1                 ! 018
            F3X    = 0.d0                                                             ! 019
           do   i  = ia,ib                                                            ! 020
            F3X(i) = Xf_FG(min(NXX_FG+1,TxF_(i))) - Xf_FG(max(0,TxF_(i-1)))           ! 021
           end do                                                                     ! 022
!                                                                                   
!- - - - -    Y-Integral        - - - - -
!
               Jh  =(3**LevelYd-1)/2                                                  ! 023        
                N_ = NyyD                                                             ! 024    
 if(LevelYd==0) N_ = NYY_FG                                                           ! 025    
                ja =    1-Jh;  if(TyV_(ja)+Jh <        1) ja = ja + 1                 ! 026    
                jb = N_+1+Jh;  if(TyV_(jb)-Jh > NYY_FG+1) jb = jb - 1                 ! 027    
            F3Y    = 0.d0                                                             ! 028    
           do   j  = ja,jb                                                            ! 030    
            F3Y(j) = Yf_FG(min(NYY_FG+1,TyF_(j))) - Yf_FG(max(0,TyF_(j-1)))           ! 031    
           end do                                                                     ! 032    
!                                                                                           
!- - - - -    Z-Integral        - - - - -
!
               Kh  =(3**LevelZd-1)/2                                                  ! 033       
                N_ = NzzD                                                             ! 034
 if(LevelZd==0) N_ = NZZ_FG                                                           ! 035
                ka =    1-Kh;  if(TzV_(ka)+Kh <        1) ka = ka + 1                 ! 036
                kb = N_+1+Kh;  if(TzV_(kb)-Kh > NZZ_FG+1) kb = kb - 1                 ! 037
            F3Z    = 0.d0                                                             ! 038
           do   k  = ka,kb                                                            ! 039
            F3Z(k) = Zf_FG(min(NZZ_FG+1,TzF_(k))) - Zf_FG(max(0,TzF_(k-1)))           ! 040
           end do                                                                     ! 041
          
                J3 = 0.d0                                                             ! 042
      do        k  = 1,NzzD+1                                                         ! 043
       do     j    = 1,NyyD+1                                                         ! 044
        do  i      = 1,NxxD+1                                                         ! 045
                                                                                      ! 046
      do        k_ = TzV_(k)-Kh,TzV_(k)+Kh                                            ! 047   
      if((1.le.k_).and.(k_.le.NZZ_FG+1))                                       then   ! 048  
                                              h_z_ = Zf_FG(k_) - Zf_FG(k_-1)          ! 049               
       do     j_   = TyV_(j)-Jh,TyV_(j)+Jh                                            ! 050   
       if((1.le.j_).and.(j_.le.NYY_FG+1))                                      then   ! 051   
                                              h_y_ = Yf_FG(j_) - Yf_FG(j_-1)          ! 052   
        do  i_     = TxV_(i)-Ih,TxV_(i)+Ih                                            ! 053   
        if((1.le.i_).and.(i_.le.NXX_FG+1))                                     then   ! 054   
                                              h_x_ = Xf_FG(i_) - Xf_FG(i_-1)          ! 055   
                                                                                           
         J3(i,j,k)= J3(i,j,k) + h_x_*h_y_*h_z_*RHSF(i_,j_,k_)                         ! 056  

        end if                                                                        ! 057  
        end do                                                                        ! 058   
       end if                                                                         ! 059   
       end do                                                                         ! 060   
      end if                                                                          ! 061   
      end do                                                                          ! 062   
                                                                                 
     end do                                                                           ! 063   
    end do                                                                            ! 064   
   end do                                                                             ! 065   

              h_x  = 1.d0/dfloat(NXX_FG)                                                  ! 066   
              h_y  = 1.d0/dfloat(NYY_FG)                                                  ! 067
              h_z  = 1.d0/dfloat(NZZ_FG)                                                  ! 068
              h_Lx = 1.d0/(h_x*dfloat(3**LevelXd))**2                                     ! 069        
              h_Ly = 1.d0/(h_y*dfloat(3**LevelYd))**2                                     ! 070 
              h_Lz = 1.d0/(h_z*dfloat(3**LevelZd))**2                                     ! 071 
              h___ = h_x*dfloat(3**LevelXd)*h_y*dfloat(3**LevelYd)*h_z*dfloat(3**LevelZd) ! 072 
                
      do        k  = 1,NzzD+1; KV_m = TzV_(k-1); KV_0 = TzV_(k); KV_p = TzV_(k+1)         ! 073   
       do     j    = 1,NyyD+1; JV_m = TyV_(j-1); JV_0 = TyV_(j); JV_p = TyV_(j+1)         ! 074 
        do  i      = 1,NxxD+1; IV_m = TxV_(i-1); IV_0 = TxV_(i); IV_p = TxV_(i+1)         ! 075     

        if((IV_0==       1).or.(JV_0==       1).or.(KV_0==       1).or. &                 ! 076 
           (IV_0==NXX_FG+1).or.(JV_0==NYY_FG+1).or.(KV_0==NZZ_FG+1)) then                 ! 077 
            
        J3(i,j,k)  = 0.d0                                                                 ! 078     

        else                                                                              ! 079 
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  X-Direction                   
         if( i==1)                                                   then                 ! 080 
                xi =      Xv_FG(IV_0) *dsqrt(h_Lx)                                        ! 081 
              AX_m = 0.d0                                                                 ! 082 
              AX_o = 2.d0          /xi                                                    ! 083 
              AX_p = 2.d0/(xi+1.d0)                                                       ! 084 
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 085 
                   *             BUX0(     JV_0,KV_0)                                     ! 086 
            d2udx2 =(      -AX_o*hatU(IV_0,JV_0,KV_0)          &                          ! 087     
                           +AX_p*hatU(IV_p,JV_0,KV_0) ) * h_Lx + Sx                       ! 088      
         end if                                                                           ! 089
         if((1<i).and.(i<NxxD+1))                                    then                 ! 090 
            d2udx2 = 1.d0*(      hatU(IV_m,JV_0,KV_0)          &                          ! 091 
                           -2.d0*hatU(IV_0,JV_0,KV_0)          &                          ! 092 
                           +     hatU(IV_p,JV_0,KV_0) ) * h_Lx                            ! 093 
         end if                                                                           ! 094 
         if(           i==NxxD+1)                                    then                 ! 095 
                xi =(1.d0-Xv_FG(IV_0))*dsqrt(h_Lx)                                        ! 096 
              AX_p = 0.d0                                                                 ! 097 
              AX_o = 2.d0          /xi                                                    ! 098 
              AX_m = 2.d0/(xi+1.d0)                                                       ! 099  
              Sx   = 2.d0/(xi+1.d0)/xi                  * h_Lx &                          ! 100  
                   *             BUX1(     JV_0,KV_0)                                     ! 101  
            d2udx2 =(       AX_m*hatU(IV_m,JV_0,KV_0)          &                          ! 102  
                           -AX_o*hatU(IV_0,JV_0,KV_0) ) * h_Lx + Sx                       ! 103  
         end if                                                                           ! 104  
!- - - - - - - - - - - - - - - - - - - - - - - - - - -  Y-Direction                         
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
              Sz   = 2.d0/(xi+1.d0)/ xi                 * h_Lz &                          ! 150
                   *             BUZ1(IV_0,JV_0     )                                     ! 151
            d2udz2 =(       AZ_m*hatU(IV_0,JV_0,KV_m)          &                          ! 152  
                           -AZ_o*hatU(IV_0,JV_0,KV_0) ) * h_Lz + Sz                       ! 153  
         end if                                                                           ! 154  
        J3(i,j,k)  = J3(i,j,k) - (d2udx2 + d2udy2 + d2udz2)*h___                          ! 155        
        end if                                                                            ! 157  
        end do                                                                            ! 158  
       end do                                                                             ! 159  
      end do                                                                              ! 160
      call RHSI3D(LevelX,LevelY,LevelZ,NxxD,NyyD,NzzD,F3X,F3Y,F3Z,J3,ia,ib,ja,jb,ka,kb)   ! 161
 deallocate(F3X,F3Y,F3Z)                                                                  ! 162
end subroutine Matrix_and_Vector_in_Mini_Cycle  
    
    
subroutine    Mini_Cycle(LeXmax,GridsXd,NxxD,TxV_,TxF_, &
                         LeYmax,GridsYd,NyyD,TyV_,TyF_, &
                         LeZmax,GridsZd,NzzD,TzV_,TzF_, &
                         MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)

INTEGER,      ALLOCATABLE :: MG_v(:,:),  Poin_V(:,:,:),  N_MGa(:,:,:)                 ! 001
INTEGER,      ALLOCATABLE :: MG_c(:,:),  Poin_C(:,:,:),  N_MGi(:,:,:),  Istart(:,:,:) ! 002
INTEGER,      ALLOCATABLE ::  TxV(:),       TxF(:),       TxV_(:),        TxF_(:)     ! 003
INTEGER,      ALLOCATABLE ::  TyV(:),       TyF(:),       TyV_(:),        TyF_(:)     ! 004
INTEGER,      ALLOCATABLE ::  TzV(:),       TzF(:),       TzV_(:),        TzF_(:)     ! 005
REAL*8 ,      ALLOCATABLE ::   J3(:,:,:),    cu(:,:,:)                                ! 006
REAL*8 ,      ALLOCATABLE :: BCX0(:,:),    BCX1(:,:)                                  ! 007
REAL*8 ,      ALLOCATABLE :: BCY0(:,:),    BCY1(:,:)                                  ! 008
REAL*8 ,      ALLOCATABLE :: BCZ0(:,:),    BCZ1(:,:)                                  ! 009
REAL*8 ,      ALLOCATABLE :: BUX0(:,:),    BUX1(:,:)                                  ! 010
REAL*8 ,      ALLOCATABLE :: BUY0(:,:),    BUY1(:,:)                                  ! 011
REAL*8 ,      ALLOCATABLE :: BUZ0(:,:),    BUZ1(:,:)                                  ! 012

INTEGER       GridsXd,GridsYd,GridsZd,LevelX,LevelY,LevelZ                            ! 013
INTEGER       GridsX ,GridsY ,GridsZ ,LevelCC                                         ! 014
INTEGER       Nx,NxA,topX,botX,Ny,NyA,topY,botY,Nz,NzA,topZ,botZ                      ! 015
REAL*8        x_a,x_b,y_a,y_b,z_a,z_b,Exact_In,Error,ErrMAXG                          ! 016    

              II =(3**LevelXd-1)/2                                                    ! 017 
              JJ =(3**LevelYd-1)/2                                                    ! 018 
              KK =(3**LevelZd-1)/2                                                    ! 019 
    allocate( J3(-II:NxxD+2+II,-JJ:NyyD+2+JJ,-KK:NzzD+2+KK) )                         ! 020               
    allocate( cu(-II:NxxD+2+II,-JJ:NyyD+2+JJ,-KK:NzzD+2+KK) )                         ! 021                   

              ka =          -3**LevelM                                                ! 022
              kb = NXYZmax+2+3**LevelM                                                ! 023
    allocate( TxV(ka:kb), TxF(ka:kb) )                                                ! 024
    allocate( TyV(ka:kb), TyF(ka:kb) )                                                ! 025
    allocate( TzV(ka:kb), TzF(ka:kb) )                                                ! 026
    
    allocate( BCX0(           0:NYY_FG+2,0:NZZ_FG+2),  BCX1(           0:NYY_FG+2,0:NZZ_FG+2) ) 
    allocate( BCY0(0:NXX_FG+2,           0:NZZ_FG+2),  BCY1(0:NXX_FG+2,           0:NZZ_FG+2) )
    allocate( BCZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BCZ1(0:NXX_FG+2,0:NYY_FG+2           ) ) 
    allocate( BUX0(           0:NYY_FG+2,0:NZZ_FG+2),  BUX1(           0:NYY_FG+2,0:NZZ_FG+2) ) 
    allocate( BUY0(0:NXX_FG+2,           0:NZZ_FG+2),  BUY1(0:NXX_FG+2,           0:NZZ_FG+2) )
    allocate( BUZ0(0:NXX_FG+2,0:NYY_FG+2           ),  BUZ1(0:NXX_FG+2,0:NYY_FG+2           ) )
 
              BCX0 = B_X0; BCX1 = B_X1;   BUX0 = 0.d0; BUX1 = 0.d0                    ! 027      
              BCY0 = B_Y0; BCY1 = B_Y1;   BUY0 = 0.d0; BUY1 = 0.d0                    ! 028
              BCZ0 = B_Z0; BCZ1 = B_Z1;   BUZ0 = 0.d0; BUZ1 = 0.d0                    ! 029
                                                                   
do             MGI = 1,MGI_DC                                                         ! 030 
 do         LevelCC= max(LeXmax,LeYmax,LeZmax        ),0,-1                           ! 031 
     
   
     
            LevelX = min(LeXmax              ,LevelCC)                                ! 032
            LevelY = min(       LeYmax       ,LevelCC)                                ! 033
            LevelZ = min(              LeZmax,LevelCC)                                ! 034
      
       call Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_,      &     ! 035
                                            NyyD,LevelY,TyV,TyF,TyV_,TyF_,      &     ! 035
                                            NzzD,LevelZ,TzV,TzF,TzV_,TzF_,J3,   &     ! 035
                                            BUX0, BUX1, BUY0, BUY1, BUZ0, BUZ1)       ! 035
           
  do        GridsX = 1,3**LevelX                                                      ! 036
       call Index_Mapping_X('-',LevelX,GridsX,Nx,NxA,topX,botX,TxV,TxF,TxV_,TxF_, &   ! 037
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 037
   do       GridsY = 1,3**LevelY                                                      ! 038
       call Index_Mapping_Y('-',LevelY,GridsY,Ny,NyA,topY,botY,TyV,TyF,TyV_,TyF_, &   ! 039
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 039
    do      GridsZ = 1,3**LevelZ                                                      ! 040
       call Index_Mapping_Z('-',LevelZ,GridsZ,Nz,NzA,topZ,botZ,TzV,TzF,TzV_,TzF_, &   ! 041
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 041

       call Vanka_type_Smoother(LevelX,LeXmax,TxV,TxF,TxV_,TxF_,Nx, &                 ! 042
                                LevelY,LeYmax,TyV,TyF,TyV_,TyF_,Ny, &                 ! 042
                                LevelZ,LeZmax,TzV,TzF,TzV_,TzF_,Nz, &                 ! 042
                                BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)                  ! 042
       
    end do                                                                            ! 043      
   end do                                                                             ! 044         
  end do                                                                              ! 045
 end do                                                                               ! 046

 do             k  = 1,NzzD+1;                 KK  = TzV_(k)          ! 047    
  do         j     = 1,NyyD+1;              JJ     = TyV_(j)          ! 048          
   do     i        = 1,NxxD+1;           II        = TxV_(i)          ! 049   
    hatU(II,JJ,KK) = hatU(II,JJ,KK) + cu(i,j,k)                       ! 050   
                                      cu(i,j,k)    = 0.d0             ! 051   
   end do                                                             ! 052   
  end do                                                              ! 053   
 end do                                                               ! 054  
 
    do         k   = 1,NzzD+1;                          KK  = TzV_(k) ! 055                                                                       
     do      j     = 1,NyyD+1;                       JJ     = TyV_(j) ! 056   
       BUX0(JJ,KK) = BUX0(JJ,KK) + BCX0(JJ,KK); BCX0(JJ,KK) = 0.d0    ! 057   
       BUX1(JJ,KK) = BUX1(JJ,KK) + BCX1(JJ,KK); BCX1(JJ,KK) = 0.d0    ! 058   
     end do                                                           ! 059   
    end do                                                            ! 060   
    
    do         k   = 1,NzzD+1;                          KK  = TzV_(k) ! 061
     do      i     = 1,NxxD+1;                       II     = TxV_(i) ! 062   
       BUY0(II,KK) = BUY0(II,KK) + BCY0(II,KK); BCY0(II,KK) = 0.d0    ! 063    
       BUY1(II,KK) = BUY1(II,KK) + BCY1(II,KK); BCY1(II,KK) = 0.d0    ! 064    
     end do                                                           ! 065 
    end do                                                            ! 066
     
   do          j   = 1,NyyD+1;                          JJ  = TyV_(j) ! 067 
    do       i     = 1,NxxD+1;                       II     = TxV_(i) ! 068 
       BUZ0(II,JJ) = BUZ0(II,JJ) + BCZ0(II,JJ); BCZ0(II,JJ) = 0.d0    ! 069 
       BUZ1(II,JJ) = BUZ1(II,JJ) + BCZ1(II,JJ); BCZ1(II,JJ) = 0.d0    ! 070 
    end do                                                            ! 071
   end do                                                             ! 072      
 
end do                                                                ! 073
deallocate(J3,cu)
deallocate(BCX0, BCX1, BCY0, BCY1, BCZ0, BCZ1, BUX0, BUX1, BUY0, BUY1, BUZ0, BUZ1)
end subroutine  Mini_Cycle 
                             
                              
subroutine  Convergence_Test
real*8      u_a,Error,h,D2UDx2,D2UDy2,D2UDz2,Res,h_x,h_y_h_z            ! 001 
               h_x = 1.d0/dfloat(NXX_FG)                                ! 002 
               h_y = 1.d0/dfloat(NYY_FG)                                ! 003 
               h_z = 1.d0/dfloat(NZZ_FG)                                ! 004 
            ResMAX = 0.d0                                               ! 005          
            ErrMAX = 0.d0                                               ! 006    
 do             k  = 2,NZZ_FG                                           ! 007    
  do          j    = 2,NYY_FG                                           ! 008             
   do       i      = 2,NXX_FG                                           ! 009    
               u_a = dexp(a_*Xv_FG(i) + b_*Yv_FG(j) + c_*Zv_FG(k))      ! 010    
             Error = dabs(u_a - hatU(i,j,k))                            ! 011             
          if(Error > ErrMAX)    ErrMAX = Error                          ! 012    
                 h = hatU(i  ,j  ,k  )                         *2.d0    ! 013    
            D2uDx2 =(hatU(i-1,j  ,k  ) - h + hatU(i+1,j  ,k  ))/h_x**2  ! 014                
            D2uDy2 =(hatU(i  ,j-1,k  ) - h + hatU(i  ,j+1,k  ))/h_y**2  ! 015                
            D2uDz2 =(hatU(i  ,j  ,k-1) - h + hatU(i  ,j  ,k+1))/h_z**2  ! 016                
               Res = dabs(RHSF(i,j,k) - D2uDx2 - D2uDy2 - D2uDz2)       ! 017       
            if(Res > ResMAX)    ResMAX = Res                            ! 018       
   end do                                                               ! 019       
  end do                                                                ! 020       
 end do                                                                 ! 021
end subroutine Convergence_Test                                         
                                                                        
subroutine  Vanka_type_Smoother(LevelX,LeXmax,TxV,TxF,TxV_,TxF_,Nx, &   
                                LevelY,LeYmax,TyV,TyF,TyV_,TyF_,Ny, &   
                                LevelZ,LeZmax,TzV,TzF,TzV_,TzF_,Nz, &   
                                BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3) 
implicit    real*8 (a-h,o-z)
integer     X_block,Y_block,Z_block                                                   ! 001
real*8      h_x,h_y,h_z,h_Lx,h_Ly,h_Lz                                                ! 002

REAL*8 ,    ALLOCATABLE ::   J3(:,:,:),             cu(:,:,:)                         ! 003
INTEGER,    ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)                 ! 004
INTEGER,    ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)                 ! 005
INTEGER,    ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:)                 ! 006
REAL*8 ,    ALLOCATABLE :: BCX0(:,:), BCX1(:,:)                                       ! 007
REAL*8 ,    ALLOCATABLE :: BCY0(:,:), BCY1(:,:)                                       ! 008
REAL*8 ,    ALLOCATABLE :: BCZ0(:,:), BCZ1(:,:)                                       ! 009

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
    
 if(NX_block_ > Nx+1) NX_block_ = Nx+1                                                ! 021  
 if(NY_block_ > Ny+1) NY_block_ = Ny+1                                                ! 022
 if(NZ_block_ > Nz+1) NZ_block_ = Nz+1                                                ! 023  

if(((Nx+1)/NX_block_)*NX_block_== Nx+1)                                        then   ! 024  
                       X_block  =(Nx+1)/NX_block_                                     ! 025
else                                                                                  ! 026
                       X_block  =(Nx+1)/NX_block_ + 1                                 ! 027
end if                                                                                ! 028
if(((Ny+1)/NY_block_)*NY_block_== Ny+1)                                        then   ! 029
                       Y_block  =(Ny+1)/NY_block_                                     ! 030     
else                                                                                  ! 031
                       Y_block  =(Ny+1)/NY_block_ + 1                                 ! 032   
end if                                                                                ! 033
if(((Nz+1)/NZ_block_)*NZ_block_== Nz+1)                                        then   ! 034
                       Z_block  =(Nz+1)/NZ_block_                                     ! 035  
else                                                                                  ! 036
                       Z_block  =(Nz+1)/NZ_block_ + 1                                 ! 037 
end if                                                                                ! 038 

if((levelX==LeXmax).and.(levelY==LeYmax).and.(levelZ==LeZmax))       NSIL_ = 2*NSIL   ! 039
if(X_block*Y_block*Z_block==1)                                                 then   ! 040
  do I_block  = 1,X_block                                                             ! 041
                           I_top =(I_block-1)*NX_block_+1                             ! 042 
                           I_bot = I_block   *NX_block_                               ! 043  
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 044
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 045  

  do J_block  = 1,Y_block                                                             ! 046
                           J_top =(J_block-1)*NY_block_+1                             ! 047 
                           J_bot = J_block   *NY_block_                               ! 048
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 049  
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 050  

  do K_block  = 1,Z_block                                                             ! 051
                           K_top =(K_block-1)*NZ_block_+1                             ! 052
                           K_bot = K_block   *NZ_block_                               ! 053   
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 054
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 055   

      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 056
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 056
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 056
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)                       ! 056
  end do                                                                              ! 057     
  end do                                                                              ! 058     
  end do                                                                              ! 059         
else                                                                                  ! 060
 do Iteration = 1,NSIL_                                                               ! 061
  do I_block  = 1,X_block                                                             ! 062
                           I_top =(I_block-1)*NX_block_+1                             ! 063 
                           I_bot = I_block   *NX_block_                               ! 064 
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 065
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 066  

  do J_block  = 1,Y_block                                                             ! 067
                           J_top =(J_block-1)*NY_block_+1                             ! 068 
                           J_bot = J_block   *NY_block_                               ! 069
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 070
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 071  

  do K_block  = 1,Z_block                                                             ! 072
                           K_top =(K_block-1)*NZ_block_+1                             ! 073
                           K_bot = K_block   *NZ_block_                               ! 074   
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 075
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 076   

      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 077
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 077
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 077
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)                       ! 078
  end do                                                                              ! 079     
  end do                                                                              ! 080     
  end do                                                                              ! 081         
  
  do K_block  = 1,Z_block                                                             ! 082
                           K_top =(K_block-1)*NZ_block_+1                             ! 083
                           K_bot = K_block   *NZ_block_                               ! 084   
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 085
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 086   
 
  do I_block  = 1,X_block                                                             ! 087
                           I_top =(I_block-1)*NX_block_+1                             ! 088 
                           I_bot = I_block   *NX_block_                               ! 089  
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 090
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 091  

  do J_block  = 1,Y_block                                                             ! 092  
                           J_top =(J_block-1)*NY_block_+1                             ! 093 
                           J_bot = J_block   *NY_block_                               ! 094
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 095
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 096  
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 097
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 097
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 097
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)                       ! 097
  end do                                                                              ! 098     
  end do                                                                              ! 099     
  end do                                                                              ! 100
  
  do J_block  = 1,Y_block                                                             ! 101
                           J_top =(J_block-1)*NY_block_+1                             ! 102 
                           J_bot = J_block   *NY_block_                               ! 103
  if(J_block ==   Y_block) J_top = Ny+2      -NY_block_                               ! 104
  if(J_block ==   Y_block) J_bot = Ny+1                                               ! 105  
  
  do K_block  = 1,Z_block                                                             ! 106
                           K_top =(K_block-1)*NZ_block_+1                             ! 107
                           K_bot = K_block   *NZ_block_                               ! 108   
  if(K_block ==   Z_block) K_top = Nz+2      -NZ_block_                               ! 109
  if(K_block ==   Z_block) K_bot = Nz+1                                               ! 110   
 
  do I_block  = 1,X_block                                                             ! 111
                           I_top =(I_block-1)*NX_block_+1                             ! 112 
                           I_bot = I_block   *NX_block_                               ! 113  
  if(I_block ==   X_block) I_top = Nx+2      -NX_block_                               ! 114
  if(I_block ==   X_block) I_bot = Nx+1                                               ! 115  
      call Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx, &               ! 116
                           TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly, &               ! 116
                           TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz, &               ! 116
                           BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)                       ! 116
  end do                                                                              ! 117     
  end do                                                                              ! 118     
  end do                                                                              ! 119              
  
 end do                                                                               ! 120 
end if                                                                                ! 121
end subroutine  Vanka_type_Smoother                                

                                
SUBROUTINE  Vanka_iteration(TxV,TxF,TxV_,TxF_,Nx,I_top,I_bot,h_x,h_Lx,  &
                            TyV,TyF,TyV_,TyF_,Ny,J_top,J_bot,h_y,h_Ly,  &     
                            TzV,TzF,TzV_,TzF_,Nz,K_top,K_bot,h_z,h_Lz,  &
                            BCX0,BCX1,BCY0,BCY1,BCZ0,BCZ1,cu,J3)
integer     M_A_P_,way
real*8      A_W,A_E,A_S,A_N,A_U,A_D,A_Px,A_Py,A_Pz,S_x,S_y,S_z                        ! 001   
real*8      h_x,h_y,h_z,h_Lx,h_Ly,h_Lz,h, H_W,H_E,H_S,H_N,H_U,H_D,xi                  ! 002
REAL*8,     ALLOCATABLE ::   AG(:,:),   xG(:)                                         ! 003
REAL*8 ,    ALLOCATABLE ::   J3(:,:,:),             cu(:,:,:)                         ! 004
INTEGER,    ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)                 ! 005
INTEGER,    ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)                 ! 006
INTEGER,    ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:),   M_A_P(:,:,:) ! 007
REAL*8 ,    ALLOCATABLE :: BCX0(:,:), BCX1(:,:)                                       ! 008
REAL*8 ,    ALLOCATABLE :: BCY0(:,:), BCY1(:,:)                                       ! 009
REAL*8 ,    ALLOCATABLE :: BCZ0(:,:), BCZ1(:,:)                                       ! 010

        i = NX_block*NY_block*NZ_block;  allocate(AG(i,i+1),xG(i))                    ! 011
       AG = 0.D0; xG = 0.D0;  allocate(M_A_P(I_top:I_bot,J_top:J_bot,K_top:K_bot))    ! 012

        m = 0                                                                         ! 013
   do   k = K_top,K_bot                                                               ! 014         
    do  j = J_top,J_bot                                                               ! 015         
     do i = I_top,I_bot                                                               ! 016
        m = m+1                                                                       ! 017 
                                       M_A_P(i,j,k) = m                               ! 018
     end do                                                                           ! 019  
    end do                                                                            ! 020
   end do                                                                             ! 021
       N_ = m                                                                         ! 022

        m = 0                                                                         ! 023
   do   k = K_top,K_bot;  KK = TzV_(TzV(k));    KKK = TzV(k)                          ! 024
    do  j = J_top,J_bot;  JJ = TyV_(TyV(j));    JJJ = TyV(j)                          ! 025
     do i = I_top,I_bot;  II = TxV_(TxV(i));    III = TxV(i)                          ! 026
        m = m+1                                                                       ! 027  
      way = 0;                               M_A_P_ = M_A_P(i,j,k)                    ! 028 

! Plane x = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==   1).and.(II==       1))                                         then   ! 029
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 030
                                       AG(m,N_+1)   = BCX0(   JJ,KK)                  ! 031   
      end if                                                                          ! 032
! Plane x = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((i==Nx+1).and.(II==NXX_FG+1))                                         then   ! 033
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 034
                                       AG(m,N_+1)   = BCX1(   JJ,KK)                  ! 035
      end if                                                                          ! 036
! Plane y = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==   1).and.(JJ==       1))                                         then   ! 037
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 038
                                       AG(m,N_+1)   = BCY0(II,   KK)                  ! 039
      end if                                                                          ! 040
! Plane y = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((j==Ny+1).and.(JJ==NYY_FG+1))                                         then   ! 041
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 042
                                       AG(m,N_+1)   = BCY1(II,   KK)                  ! 043
      end if                                                                          ! 044
! Plane z = 0 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==   1).and.(KK==       1))                                         then   ! 045
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 046
                                       AG(m,N_+1)   = BCZ0(II,JJ   )                  ! 047
      end if                                                                          ! 048
! Plane z = 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if((k==Nz+1).and.(KK==NZZ_FG+1))                                         then   ! 049
                                       AG(m,M_A_P_) = 1.D0;           way = 1         ! 050
                                       AG(m,N_+1)   = BCZ1(II,JJ   )                  ! 051 
      end if                                                                          ! 052

      if(way == 0)                                                             then   ! 053  

     H_W = 0.d0; H_E = 0.d0; S_x = 0.d0                                               ! 054       
   if(i == 1)                                                                  then   ! 055
     xi  =          Xv_FG(II)             *dsqrt(h_Lx)                                ! 056
    A_W  =  0.d0                                                                      ! 057
    A_Px =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lx                                 ! 058
    A_E  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lx; H_E = A_E*cu(TxV(i+1),JJJ,KKK) ! 059
    S_x  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lx      *   BCX0(          JJ, KK) ! 060
   end if                                                                             ! 061 
   if((1<i).and.(i<Nx+1))                                                      then   ! 062 
    A_W  =                                       h_Lx; H_W = A_W*cu(TxV(i-1),JJJ,KKK) ! 063 
    A_E  =                                       h_Lx; H_E = A_E*cu(TxV(i+1),JJJ,KKK) ! 064 
    A_Px =-A_W-A_E                                                                    ! 065 
   end if                                                                             ! 066 
   if(i == Nx+1)                                                               then   ! 067      
      xi = (1.d0 -  Xv_FG(II))            *dsqrt(h_Lx)                                ! 068 
    A_W  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lx; H_W = A_W*cu(TxV(i-1),JJJ,KKK) ! 069 
    A_Px =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lx                                 ! 070 
    A_E  =  0.d0                                                                      ! 071 
    S_x  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lx      *   BCX1(          JJ, KK) ! 072       
   end if                                                                             ! 073
                                                                                       
     H_S = 0.d0; H_N = 0.d0; S_y = 0.d0                                               ! 074 
   if(j == 1)                                                                  then   ! 075 
     xi  =          Yv_FG(JJ)             *dsqrt(h_Ly)                                ! 076 
    A_S  =  0.d0                                                                      ! 077         
    A_Py =-(1.d0 - (xi-1.d0)          /xi)* 2.d0*h_Ly                                 ! 077   
    A_N  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Ly; H_N = A_N*cu(III,TyV(j+1),KKK) ! 077   
    S_y  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Ly      *   BCY0( II,          KK) ! 078   
   end if                                                                             ! 079   
   if((1<j).and.(j<Ny+1))                                                      then   ! 080   
    A_S  =                                       h_Ly; H_S = A_S*cu(III,TyV(j-1),KKK) ! 081   
    A_N  =                                       h_Ly; H_N = A_N*cu(III,TyV(j+1),KKK) ! 082              
    A_Py =-A_S-A_N                                                                    ! 083      
   end if                                                                             ! 084      
   if(j == Ny+1)                                                               then   ! 085      
     xi  = (1.d0 -  Yv_FG(JJ))            *dsqrt(h_Ly)                                ! 086      
    A_S  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Ly; H_S = A_S*cu(III,TyV(j-1),KKK) ! 087               
    A_Py =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Ly                                 ! 088      
    A_N  =  0.d0                                                                      ! 089      
    S_y  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Ly      *   BCY1( II,          KK) ! 090      
   end if                                                                             ! 091      
                                                                                                 
     H_D = 0.d0; H_U = 0.d0; S_z = 0.d0                                               ! 092      
   if(k == 1)                                                                  then   ! 093      
     xi  =          Zv_FG(KK)             *dsqrt(h_Lz)                                ! 094      
    A_D  =  0.d0                                                                      ! 095      
    A_Pz =-(1.d0 - (xi-1.d0)          /xi)* 2.d0*h_Lz                                 ! 096      
    A_U  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lz; H_U = A_U*cu(III,JJJ,TzV(k+1)) ! 097      
    S_z  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lz      *   BCZ0( II, JJ         ) ! 097      
   end if                                                                             ! 097      
   if((1<k).and.(k<Nz+1))                                                      then   ! 097      
    A_D  =                                       h_Lz; H_D = A_D*cu(III,JJJ,TzV(k-1)) ! 098      
    A_U  =                                       h_Lz; H_U = A_U*cu(III,JJJ,TzV(k+1)) ! 099      
    A_Pz =-A_D-A_U                                                                    ! 100      
   end if                                                                             ! 101            
   if(k == Nz+1)                                                               then   ! 102      
     xi  = (1.d0 -  Zv_FG(KK))            *dsqrt(h_Lz)                                ! 103      
    A_D  = (1.d0 - (xi-1.d0)/(xi+1.d0)   )      *h_Lz; H_D = A_D*cu(III,JJJ,TzV(k-1)) ! 104      
    A_Pz =-(1.d0 - (xi-1.d0)/          xi)* 2.d0*h_Lz                                 ! 105      
    A_U  =  0.d0                                                                      ! 106                
    S_z  =  1.d0            /(xi+1.d0)/xi * 2.d0*h_Lz      *   BCZ1( II, JJ         ) ! 107      
   end if                                                                             ! 108      
                                                                                            
       if(i-1 >= I_top)     AG(m,M_A_P(i-1,j  ,k  )) = A_W                            ! 109      
       if(i-1 >= I_top)                                H_W = 0.D0                     ! 110 
                                                                                           
       if(i+1 <= I_bot)     AG(m,M_A_P(i+1,j  ,k  )) = A_E                            ! 111                                 
       if(i+1 <= I_bot)                                H_E = 0.D0                     ! 112 
                                                                                       
       if(j-1 >= J_top)     AG(m,M_A_P(i  ,j-1,k  )) = A_S                            ! 113  
       if(j-1 >= J_top)                                H_S = 0.D0                     ! 114 
                                                                                      
       if(j+1 <= J_bot)     AG(m,M_A_P(i  ,j+1,k  )) = A_N                            ! 115   
       if(j+1 <= J_bot)                                H_N = 0.D0                     ! 116   
                                                                                      
       if(k-1 >= K_top)     AG(m,M_A_P(i  ,j  ,k-1)) = A_D                            ! 116     
       if(k-1 >= K_top)                                H_D = 0.D0                     ! 116     
                                                                                      
       if(k+1 <= K_bot)     AG(m,M_A_P(i  ,j  ,k+1)) = A_U                            ! 116     
       if(k+1 <= K_bot)                                H_U = 0.D0                     ! 117     
                            AG(m,M_A_P(i  ,j  ,k  )) = A_Px + A_Py + A_Pz             ! 118     
                            AG(m,N_+1)               = J3(III,JJJ,KKK) &              ! 119     
                          - S_x - S_y - S_z - H_W - H_E - H_S - H_N - H_D - H_U       ! 119          
      end if                                                                          ! 120     
     end do                                                                           ! 121     
    end do                                                                            ! 122               
   end do                                                                             ! 123               
                            call GAUSSIAN_ELIMINATION(AG,xG,N_)                       ! 124               
        m = 0                                                                         ! 125               
   do   k = K_top,K_bot                                                               ! 126     
    do  j = J_top,J_bot                                                               ! 127              
     do i = I_top,I_bot                                                               ! 128     
        m = m+1                                                                       ! 129     
                            cu(TxV(i),TyV(j),TzV(k)) = xG(m)                          ! 130     
     end do                                                                           ! 131     
    end do                                                                            ! 132     
   end do                                                                             ! 133     
                            deallocate(AG,xG,M_A_P)                                   ! 134     
END SUBROUTINE   Vanka_iteration

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

END MODULE Parallel_Coarse_Grid_Smoothing_2 