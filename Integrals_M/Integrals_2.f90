MODULE  Integrals_2
USE     Integrals_1
USE     OMP_LIB           
USE     RMT_3D_2020_OpenMP
    
PRIVATE    
PUBLIC ::     Starting_Guess_and_Boundary_Conditions
PUBLIC ::     Matrix_and_Vector_in_Mini_Cycle, Mini_Cycle

CONTAINS 
    
subroutine Starting_Guess_and_Boundary_Conditions
 allocate( RHSF(0:NXX_FG+2,0:NYY_FG+2,0:NZZ_FG+2) ); RHSF = 0.d0  ! 001
do           k  = 1,NZZ_FG+1                                      ! 002          
 do        j    = 1,NYY_FG+1                                      ! 003
  do     i      = 1,NXX_FG+1                                      ! 004
    RHSF(i,j,k) = dexp(Xv_FG(i) + Yv_FG(j) + Zv_FG(k))            ! 005 
  end do                                                          ! 006
 end do                                                           ! 007
end do                                                            ! 008
end subroutine Starting_Guess_and_Boundary_Conditions


subroutine    Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_,   &
                                              NyyD,LevelY,TyV,TyF,TyV_,TyF_,   &
                                              NzzD,LevelZ,TzV,TzF,TzV_,TzF_, J3)

REAL*8,       ALLOCATABLE ::  F3X(:),    F3Y(:),     F3Z(:),      J3(:,:,:)           ! 001   
INTEGER,      ALLOCATABLE ::  TxV(:),    TxF(:),    TxV_(:),    TxF_(:)               ! 002   
INTEGER,      ALLOCATABLE ::  TyV(:),    TyF(:),    TyV_(:),    TyF_(:)               ! 003   
INTEGER,      ALLOCATABLE ::  TzV(:),    TzF(:),    TzV_(:),    TzF_(:)               ! 004   

real*8        h_x_,h_y_,h_z_                                                          ! 005   

              ka   =          -3**LevelM                                              ! 006   
              kb   = NXYZmax+2+3**LevelM                                              ! 007   
    allocate( F3X(ka:kb), F3Y(ka:kb), F3Z(ka:kb) )                                    ! 008   
!
!- - - - -    X-Integral        - - - - -
!
               Ih  =(3**LevelXd-1)/2                                                  ! 009   
                N_ = NxxD                                                             ! 010
 if(LevelXd==0) N_ = NXX_FG                                                           ! 011
                ia =    1-Ih;  if(TxV_(ia)+Ih <        1) ia = ia + 1                 ! 012
                ib = N_+1+Ih;  if(TxV_(ib)-Ih > NXX_FG+1) ib = ib - 1                 ! 013
            F3X    = 0.d0                                                             ! 014
           do   i  = ia,ib                                                            ! 015
            F3X(i) = Xf_FG(min(NXX_FG+1,TxF_(i))) - Xf_FG(max(0,TxF_(i-1)))           ! 016
           end do                                                                     ! 017
!                                                                                   
!- - - - -    Y-Integral        - - - - -
!
               Jh  =(3**LevelYd-1)/2                                                  ! 018        
                N_ = NyyD                                                             ! 019    
 if(LevelYd==0) N_ = NYY_FG                                                           ! 020    
                ja =    1-Jh;  if(TyV_(ja)+Jh <        1) ja = ja + 1                 ! 021    
                jb = N_+1+Jh;  if(TyV_(jb)-Jh > NYY_FG+1) jb = jb - 1                 ! 022    
            F3Y    = 0.d0                                                             ! 023    
           do   j  = ja,jb                                                            ! 024    
            F3Y(j) = Yf_FG(min(NYY_FG+1,TyF_(j))) - Yf_FG(max(0,TyF_(j-1)))           ! 025    
           end do                                                                     ! 026    
!                                                                                           
!- - - - -    Z-Integral        - - - - -
!
               Kh  =(3**LevelZd-1)/2                                                  ! 027       
                N_ = NzzD                                                             ! 028
 if(LevelZd==0) N_ = NZZ_FG                                                           ! 029
                ka =    1-Kh;  if(TzV_(ka)+Kh <        1) ka = ka + 1                 ! 030
                kb = N_+1+Kh;  if(TzV_(kb)-Kh > NZZ_FG+1) kb = kb - 1                 ! 031
            F3Z    = 0.d0                                                             ! 032
           do   k  = ka,kb                                                            ! 033
            F3Z(k) = Zf_FG(min(NZZ_FG+1,TzF_(k))) - Zf_FG(max(0,TzF_(k-1)))           ! 034
           end do                                                                     ! 035
           
                J3 = 0.d0                                                             ! 036
                II =(3**LevelXd-1)/2                                                  ! 037
                JJ =(3**LevelYd-1)/2                                                  ! 038
                KK =(3**LevelZd-1)/2                                                  ! 039
   do           k  =-KK+1,NzzD+1+KK                                                   ! 040
    do        j    =-JJ+1,NyyD+1+JJ                                                   ! 041
     do     i      =-II+1,NxxD+1+II                                                   ! 042
      do        k_ = TzV_(k)-Kh,TzV_(k)+Kh;   h_z_ = Zf_FG(k_) - Zf_FG(k_-1)          ! 043
      if(0<k_.and.k_<NZZ_FG+2)                                                 then   ! 044 
       do     j_   = TyV_(j)-Jh,TyV_(j)+Jh;   h_y_ = Yf_FG(j_) - Yf_FG(j_-1)          ! 045
       if(0<j_.and.j_<NYY_FG+2)                                                then   ! 046            
        do  i_     = TxV_(i)-Ih,TxV_(i)+Ih;   h_x_ = Xf_FG(i_) - Xf_FG(i_-1)          ! 047
        if(0<i_.and.i_<NXX_FG+2)                                               then   ! 048                        
         J3(i,j,k) = J3(i,j,k) + h_x_*h_y_*h_z_*RHSF(i_,j_,k_)                        ! 049  
        end if                                                                        ! 050
        end do                                                                        ! 051   
       end if                                                                         ! 052
       end do                                                                         ! 053   
      end if                                                                          ! 054
      end do                                                                          ! 055   
     end do                                                                           ! 056   
    end do                                                                            ! 057   
   end do                                                                             ! 058 
   call RHSI3D(LevelX,LevelY,LevelZ,NxxD,NyyD,NzzD,F3X,F3Y,F3Z,J3,ia,ib,ja,jb,ka,kb)  ! 059
   deallocate(F3X,F3Y,F3Z)                                                            ! 060
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
REAL*8 ,      ALLOCATABLE ::   J3(:,:,:)                                              ! 006

INTEGER       GridsX ,GridsY ,GridsZ ,way                                             ! 008
INTEGER       GridsXd,GridsYd,GridsZd,LevelX,LevelY,LevelZ,LevelCC                    ! 007
INTEGER       Nx,NxA,topX,botX,Ny,NyA,topY,botY,Nz,NzA,topZ,botZ                      ! 009
REAL*8        x_a,x_b,y_a,y_b,z_a,z_b,Exact_In,Error,ErrMAX                           ! 010  
character     strX*2,strY*2,strZ*2,Name*37                                            ! 011  
              II =(3**LevelXd-1)/2                                                    ! 012 
              JJ =(3**LevelYd-1)/2                                                    ! 013 
              KK =(3**LevelZd-1)/2                                                    ! 014 
    allocate( J3(-II:NxxD+2+II,-JJ:NyyD+2+JJ,-KK:NzzD+2+KK) )                         ! 015               

              ka =          -3**LevelM                                                ! 016
              kb = NXYZmax+2+3**LevelM                                                ! 017
    allocate( TxV(ka:kb), TxF(ka:kb) )                                                ! 018
    allocate( TyV(ka:kb), TyF(ka:kb) )                                                ! 019
    allocate( TzV(ka:kb), TzF(ka:kb) )                                                ! 020
                                                write(strX,'(i2)'   ) GridsXd         ! 021
                                                   if(strX(1:1)==' ')   strX(1:1)='0' ! 022
                                                   if(strX(2:2)==' ')   strX(2:2)='0' ! 023
                                                write(strY,'(i2)'   ) GridsYd         ! 024
                                                   if(strY(1:1)==' ')   strY(1:1)='0' ! 025
                                                   if(strY(2:2)==' ')   strY(2:2)='0' ! 026
                                                write(strZ,'(i2)'   ) GridsZd         ! 027
                                                   if(strZ(1:1)==' ')   strZ(1:1)='0' ! 028
                                                   if(strZ(2:2)==' ')   strZ(2:2)='0' ! 029
Name = '_GridsXd_'//strX(1:2)//' GridsYd_'//strY(1:2)//' GridsZd_'//strZ(1:2)//'.txt' ! 030 
               way = omp_get_thread_num()+10                                          ! 031 
open(way,file=Name)                                                                   ! 032  
 do         LevelCC= max(LeXmax,LeYmax,LeZmax        ),0,-1                           ! 033 
            LevelX = min(LeXmax              ,LevelCC)                                ! 034
            LevelY = min(       LeYmax       ,LevelCC)                                ! 035
            LevelZ = min(              LeZmax,LevelCC)                                ! 036
            
write(way,"(//4x,'LevelX =',i2,4x,'LevelY =',i2,4x,'LevelZ =',i2/)") LevelX,LevelY,LevelZ ! 037       
write(way,"('    GridsX   Nx   GridsY   Ny   GridsZ   Nz     ErrMAX')")               ! 038      
       call Matrix_and_Vector_in_Mini_Cycle(NxxD,LevelX,TxV,TxF,TxV_,TxF_,  &         ! 039
                                            NyyD,LevelY,TyV,TyF,TyV_,TyF_,  &         ! 040
                                            NzzD,LevelZ,TzV,TzF,TzV_,TzF_,J3)         ! 041
  do        GridsX = 1,3**LevelX                                                      ! 042
       call Index_Mapping_X('-',LevelX,GridsX,Nx,NxA,topX,botX,TxV,TxF,TxV_,TxF_, &   ! 043
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 044
   do       GridsY = 1,3**LevelY                                                      ! 045
       call Index_Mapping_Y('-',LevelY,GridsY,Ny,NyA,topY,botY,TyV,TyF,TyV_,TyF_, &   ! 046
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 047
    do      GridsZ = 1,3**LevelZ                                                      ! 048
       call Index_Mapping_Z('-',LevelZ,GridsZ,Nz,NzA,topZ,botZ,TzV,TzF,TzV_,TzF_, &   ! 049
                                MG_v, MG_c, Poin_V, N_MGa, Poin_C, N_MGi, Istart)     ! 050
       ErrMAX = 0.D0                                                                  ! 051   
     do     i = 1,Nx+1;       III = TxV(i  )                                          ! 052
          x_a = Xf_FG(max(       0, TxF_(TxF(i-1))))                                  ! 053     
          x_b = Xf_FG(min(NXX_FG+1, TxF_(TxF(i  ))))                                  ! 054     
      do    j = 1,Ny+1;       JJJ = TyV(j  )                                          ! 055
          y_a = Yf_FG(max(       0, TyF_(TyF(j-1))))                                  ! 056        
          y_b = Yf_FG(min(NYY_FG+1, TyF_(TyF(j  ))))                                  ! 057        
       do   k = 1,Nz+1;       KKK = TzV(k  )                                          ! 058  
          z_a = Zf_FG(max(       0, TzF_(TzF(k-1))))                                  ! 059
          z_b = Zf_FG(min(NZZ_FG+1, TzF_(TzF(k  ))))                                  ! 060
     Exact_In =(dexp(x_b)-dexp(x_a))/(x_b-x_a) &                                      ! 061
              *(dexp(y_b)-dexp(y_a))/(y_b-y_a) &                                      ! 062
              *(dexp(z_b)-dexp(z_a))/(z_b-z_a)                                        ! 063
        Error = dabs(Exact_In - J3(III,JJJ,KKK))                                      ! 064
     if(Error > ErrMAX) ErrMAX = Error                                                ! 065
       end do                                                                         ! 066 
      end do                                                                          ! 067 
     end do                                                                           ! 068
    write(way,"(1x,3(5x,i2,4x,i3),2x,D11.4)") GridsX,Nx,GridsY,Ny,GridsZ,Nz,ErrMAX    ! 069
    end do                                                                            ! 070
   end do                                                                             ! 071
  end do                                                                              ! 072 
 end do                                                                               ! 073   
close(way)                                                                            ! 074 
deallocate(J3)                                                                        ! 075
end subroutine  Mini_Cycle 

END MODULE Integrals_2