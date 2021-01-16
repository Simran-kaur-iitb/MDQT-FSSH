program Aj
  implicit none
  integer i,j,k,l
  double precision pi
  double precision A0, B0
  double precision mass, omega_1, omega_2, k_1, k_2, k_3
  double precision x1,x2,delta_x, q
  integer, parameter :: n=31
  double precision:: kin_en(n,n),pot(n,n)
  double precision:: mat(n,n)
  integer nold, m_values
  double precision, allocatable:: work(:), iwork(:), isuppz(:)
  double precision  eigen_value(n), eigen_vect(n,n) 

  call setup_parameters
  call main
!-------------------------------------------------------!
  contains
!-------------------------------------------------------!
 subroutine setup_parameters
   implicit none
   pi=dacos(-1.d0)                            !!----in atomic units---!!
   A0=-33367.6*4.55d-6*(1.d0/1.889)**2
   B0=133501.2*4.55d-6*(1.d0/1.889)**4
   mass=20*1822.d0
   omega_1=200*4.55d-6
   omega_2=287*4.55d-6
   k_1=1.5d-4*4.55d-6*(1.d0/1.889)**2
   k_2=-5.d-4*4.55d-6*(1.d0/1.889)**3
   k_3=-7.d-3*4.55d-6*(1.d0/1.889)**2
 end subroutine
!-------------------------------------------------------!
 subroutine main
   implicit none
   call initial_conditions
   call compute_kin_en(kin_en)
   x2=0.18*1.889d0
   do i=1,31
     x1=-0.5d0+(i-1)*10.d0/300.d0   !!in Angstroms
     x1=x1*1.889d0                  !! converted to a.u.
     write(12,*) i, x1
     call compute_potential(i,x1,x2,pot)
     mat=pot+kin_en
     !call diag(mat,n,eigen_value,eigen_vect,m_values)
     write(24,*) x1, mat
  enddo
  
 end subroutine
!-------------------------------------------------------!
 subroutine initial_conditions
   implicit none
   nold=0
   pot=0.d0
   kin_en=0.d0
   delta_x=(2.d0-(-2.d0))/32.d0
 end subroutine
!-------------------------------------------------------!
 subroutine compute_kin_en(kin_en)
   implicit none
   integer i, j, a ,b
   double precision prefactor
   double precision, intent(out) :: kin_en(:,:)
   prefactor=1.d0/(2.d0*mass*delta_x)
   do i=1,31
     do j=1,31
       if(i==j)kin_en(i,j)=prefactor*(-1)**(i-j)*((pi*pi)/3.d0) 
       if(i.ne.j)kin_en(i,j)=prefactor*(-1)**(i-j)*(2.d0/(i-j)*(i-j))
       write(10,*) i,j, kin_en(i,j) 
     enddo
   enddo
 end subroutine
!-------------------------------------------------------!
 subroutine compute_potential(i,x1,x2,pot)
   implicit none
   integer j
   double precision v0,v1,vc
   integer, intent(in)::i
   double precision, intent(in) ::x1,x2
   double precision, intent(out):: pot(:,:)
   do j=1,31
     q=-2.d0+delta_x*j
     write(23,*) q
     q=1.889d0*q
     v0=A0*q**2+B0*q**4
     v1=0.5d0*mass*(omega_1**2*x1**2+omega_2**2*x2**2) + k_3*x1*x2
     vc=k_1*q*x1 + k_2*q**2*x2
     if(i==j)pot(i,j)=v0+v1+vc
     if(i.ne.j)pot(i,j)=0.d0
     write(14,*) i,j,x1,pot(i,j)
   enddo 
 end subroutine
!-------------------------------------------------------!
subroutine diag(mat,n,eigen_value,eigen_vect,m_values)
  !! Diaganalizing matrix using dsyevr. First m_values eigen values and eigenvectors computed.
  !! The module's common variables should contain:

  !! Initialize nold=0

  !! nold makes sure that everytime value of n changes, work and iwork are re-allocated for optimal performance.
  !! mat is destroyed after use.

  implicit none
  integer,intent(in) :: n,m_values
  real*8,intent(out) :: eigen_value(n),eigen_vect(n,m_values)
  real*8,intent(inout) :: mat(n,n)
  real*8 vl,vu,abstol
  integer il,iu,info,m,AllocateStatus
  integer lwork,liwork

  vl=0.d0;vu=0.d0   !! not referenced
  il=1;iu=m_values
  abstol=0.d0
  info=0

  if(nold.ne.n .or. .not.allocated(work) .or. .not.allocated(iwork) .or. .not.allocated(isuppz)) then
  !if(nold.ne.n) then
    lwork=-1;liwork=-1
    if(allocated(isuppz))deallocate(isuppz)
    if(allocated(work))deallocate(work)
    if(allocated(iwork))deallocate(iwork)
    allocate(isuppz(2*m_values),work(n),iwork(n))
    call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
    lwork=nint(work(1)); liwork=iwork(1)
    deallocate(work,iwork)
     allocate(work(lwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    allocate(iwork(liwork),STAT=AllocateStatus)
    if(allocatestatus.ne.0) write(6,*)"problem in diag, allocation"
    nold=n
  endif

  lwork=size(work)
  liwork=size(iwork)

  call dsyevr('V','I','U',n,mat,n,vl,vu,il,iu,abstol,m,eigen_value,eigen_vect,n,isuppz,work,lwork,iwork,liwork,info)
  if(info.ne.0) then
    write(6,*) "problem in diagonalization",info
    stop
  endif

end subroutine diag
!----------------------------------------------------------
end program
