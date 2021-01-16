program alpha_fcc
  implicit none

  double precision tic, toc, pi, kb
  double precision rroot3, sigma,sigma6,sigma12, epsilon
  integer nsteps, nsteps_equilibration
  integer, parameter:: nc=4, n=4*nc**3, alpha=3 !integer n, nc !N=tot no. of molecules, NC tot no. of FCC units cells along each dirxn
  integer ix,iy,iz,iref, m
  double precision dt, tim, total_time, equilibration_time
  double precision, dimension(n+2):: rx,ry,rz,ex,ey,ez
  double precision, dimension(n)::  rx_cl, ry_cl, rz_cl
  double precision, dimension(2*n+1,3)::positions
  double precision, dimension(2*n+1)::mass,charge, A, C, neighbour
  double precision, dimension(2*n+1,3)::pos, vel,acc, pos_old, vel_old,pos_new
  double precision pot, pot_ss, pot_cs, pot_ab
  double precision cell, cell2, cell3
  double precision density, cube_len, box_len
  double precision bond_length_CH3Cl, bond_length_AH, bond_length_AB
  double precision A_ch3, A_cl, C_ch3, C_cl, mass_ch3, mass_cl
  double precision mass_A, mass_H, mass_B
  double precision ee(alpha,2)
  double precision r, r0, l
  double precision RAB,aa,bb,d_A,d_B,DA,n_A,n_B,cc
  double precision kin_energy, energy
  double precision rc, rcsq,rt, pot_rc_cs 
  double precision temperature, Average_temperature
  double precision tic1, toc1,tic2,toc2
  integer, parameter:: nbins=200
  double precision vdist(nbins), delv
  double precision boltz(2*n-2)
  !------------------------------------------------------------------!
  call cpu_time(tic)
  call setup_parameters
  call initial_conditions
  !call equilibrate
  call cpu_time(tic1)
  call evolve
  call cpu_time(toc1)
  call cpu_time(toc)
  !call boltzmann
  call distribution
  write(741,*) "time taken for evolve =",toc1-tic1
  write(741,*) "time taken for the program to run =",toc-tic
  !----------------------------------------------------------------!
  contains
  subroutine setup_parameters   ! all parameters are in atomic units
    implicit none
    integer z

    pi=dacos(-1.d0)
    kb=(8.617d0*10d-5)/27.211d0
    !dt = 2.d0
    dt = 5.d0  ! for obtaining initial configuration, time in atomic units
    total_time =8000 * dt
    !total_time =4000 * dt  !time for obtaining initial 
    equilibration_time = 2000.0 * dt
    write(6,*) dt, total_time 
    density=0.012*0.148184d0 
    temperature=247.d0    
    rroot3=1.d0/sqrt(3.d0)
    cube_len=28.d0*1.889
    box_len=(n/density)**(1/3.d0) ! calculating length of side of the unit cell/box_len
    cell=box_len/real(4)
    cell2=0.5d0*cell
    
    print *,"dt, cell, cell2, cube_len",dt/40.d0,cell, cell2, box_len, box_len/2.d0
    
    bond_length_CH3Cl =1.781d0*1.889d0
    bond_length_AH =1.d0*1.889d0
    bond_length_AB =2.7d0*1.889d0

    A_ch3=5136.1464d0
    A_cl=4173.7273d0

    C_ch3=14.1413d0
    C_cl=14.6465d0

    sigma=3.5*1.889d0
    sigma6=sigma**6
    sigma12=sigma6*sigma6
    epsilon=200.d0

    mass_ch3=15.035*1822.d0
    mass_cl=35.453*1822.d0
    mass_A=93*1822.d0
    mass_H=1*1822.d0
    mass_B=59*1822.d0
   
    r0=1.43*1.889d0
    l=0.125*1.889d0
    ee(1,1)=-0.5d0; ee(1,2)=-1.d0
    ee(2,1)=+0.5d0; ee(2,2)=+0.5d0
    ee(3,1)=+0.d0;  ee(3,2)=+0.5d0
    
    aa=11.2d0/1.889d0                       ! Hydrogen bonding potential parameters (au)
    bb=7.1e13/627.509d0
    d_A=0.95*1.889d0
    d_B=0.97*1.889d0
    DA=110.d0/627.509d0
    n_A=9.26d0/1.889d0
    n_B=11.42d0/1.889d0
    cc=0.776
    rc=13.8d0*1.889
    rcsq=rc**2
    rt=0.95*rc
    pot_rc_cs=4.d0*epsilon*((sigma12/rc**12)-(sigma6/rc**6))
    nsteps=nint(total_time/dt)  
    nsteps_equilibration=nint(equilibration_time/dt) 
    delv=(0.01d0-(-0.01d0))/real(nbins)    
 end subroutine setup_parameters

 !-----------------------------------------------------------------!
  subroutine initial_conditions
    implicit none
    !call VLJ
    !call RAB_r   !!plot V_AB vs r at R=2.7 A
    call fcc
    vel = 0.d0
    call compute_pot(pos,pot,acc)
    call compute_energy
    tim=0.0
  end subroutine initial_conditions
 !-----------------------------------------------------------------!
  subroutine equilibrate
    implicit none
    integer i, m  
    open(16,file='pot_kin.out')
    open(17,file="energy.out")
    open(18,file="vmd.xyz")

    do i=1,nsteps_equilibration
      call write_output
      m=mod(i,50) !! writing co-ordinates after 50 frames
      if(m==1) call write_vmd
      if(mod(i,500)==0) call thermal_velocities(temperature)
      call evolve_1step
      call restrict_atoms_inside_the_box
    end do
    close(16);close(17);close(18) 
    
  end subroutine equilibrate
 !-----------------------------------------------------------------!
  subroutine fcc
    implicit none

   open(10,file="X.xyz")
   open(11,file="Pos.xyz")
   open(12,file="positions.xyz")
   open(13,file="neighbour.out")
   open(14,file="charge_solvent.out")
   open(15,file="A_C.out")

   call build_unit_cell
   call build_lattice
   call write_r
   call generate_system_lattice  !will generate a lattice palcing all the molecules with appropriate bond lengths
  end subroutine fcc
 !--------------------------------------------------------------------!
  subroutine build_unit_cell
  !!http://www.ccl.net/cca/software/SOURCES/FORTRAN/allen-tildesley-book/f.23.shtml
    implicit none
    
    rx(1)=0.d0      ! building sublattice a
    ry(1)=0.d0
    rz(1)=0.d0
    ex(1)=rroot3
    ey(1)=rroot3
    ez(1)=rroot3

    rx(2)=cell2      ! building sublattice b
    ry(2)=cell2
    rz(2)=0.d0
    ex(2)=rroot3
    ey(2)=-rroot3
    ez(2)=-rroot3

    rx(3)=0.d0      ! building sublattice c
    ry(3)=cell2
    rz(3)=cell2
    ex(3)=-rroot3
    ey(3)=rroot3
    ez(3)=-rroot3
  
    rx(4)=cell2      ! building sublattice d
    ry(4)=0.d0
    rz(4)=cell2
    ex(4)=-rroot3
    ey(4)=-rroot3
    ez(4)=rroot3

  end subroutine 
 !----------------------------------------------------------------!
  subroutine build_lattice
    implicit none
    integer i 
    m=0

    do iz=1,nc
      do iy=1,nc
        do ix=1,nc
          do iref=1,4
            rx(iref+m)=rx(iref)+cell*real(ix-1)
            ry(iref+m)=ry(iref)+cell*real(iy-1)
            rz(iref+m)=rz(iref)+cell*real(iz-1)

            ex(iref+m)=ex(iref)
            ey(iref+m)=ey(iref)
            ez(iref+m)=ez(iref)
          enddo
          m=m+4
       enddo
     enddo
    enddo

    do i=1,n
      rx(i)=rx(i)-cube_len*0.5d0
      ry(i)=ry(i)-cube_len*0.5d0
      rz(i)=rz(i)-cube_len*0.5d0
    enddo

  end subroutine build_lattice
 !------------------------------------------------------------!
  subroutine write_r
    implicit none
    integer i,j

    write(10,*) n
    write(10,*)
    do i=1,n
      write(10,*) "X", rx(i), ry(i), rz(i)
    enddo
  end subroutine write_r
 !-----------------------------------------------------------!
  subroutine generate_system_lattice
    implicit none
    integer i,j,k,l
    double precision rij,rijsq,vec(3)
  
    vdist=0.d0  
    pos=0.d0;neighbour=0.d0
    A=0.d0;C=0.d0;mass=0.d0
    charge=0.d0
    
    !!open(200,file='initial-pos.xyz')     !! this was the initial file generated using the code from Allen tildesley
    !!read(200,*) pos
     open(580,file='pos19.xyz') !! this is the file obtained after resetting velocities=0, which still wasnt equilibrated,using this
                                !!and further resetting velocities =0,  
     read(580,*) pos
    j=1
    do i=1,2*(n-1),2
     !! pos(i,1)=rx(j)
     !! pos(i+1,1)=rx(j)+bond_length_CH3Cl*ex(j)

     !! pos(i,2)=ry(j)
     !! pos(i+1,2)=ry(j)+bond_length_CH3Cl*ey(j)

     !! pos(i,3)=rz(j)
     !! pos(i+1,3)=rz(j)+bond_length_CH3Cl*ez(j)

      charge(i)=+0.25
      charge(i+1)=-0.25

      A(i)=A_ch3
      A(i+1)=A_cl

      C(i)=C_ch3
      C(i+1)=C_cl

      mass(i)=mass_ch3
      mass(i+1)=mass_cl


      neighbour(i)=i+1
      neighbour(i+1)=i
      write(58,*) neighbour(i),neighbour(i+1),i,i+1
      j=j+1
    
    enddo
    
    mass(2*n-1)=mass_A
    mass(2*n)  =mass_H
    mass(2*n+1)=mass_B
    
  !-------- Chunk of code i needed initially to generate the lattice for the particles in the system--------!!!! 
   !! pos(2*n-1,1)=rx(n)     
   !! pos(2*n-1,2)=ry(n)
   !! pos(2*n-1,3)=rz(n)
   
   !! pos(2*n,1)=pos(2*n-1,1)+1.0*bond_length_AH*ex(n)
   !! pos(2*n,2)=pos(2*n-1,2)+1.0*bond_length_AH*ey(n)
   !! pos(2*n,3)=pos(2*n-1,3)+1.0*bond_length_AH*ez(n)
 
   !! pos(2*n+1,1)=pos(2*n-1,1)+1.0*bond_length_AB*ex(n)
   !! pos(2*n+1,2)=pos(2*n-1,2)+1.0*bond_length_AB*ey(n)
   !! pos(2*n+1,3)=pos(2*n-1,3)+1.0*bond_length_AB*ez(n)

    !!write(12,*)  2*(n)+1
    !!write(12,*)
    !!do i=1,2*(n-1),2
     !! write(12,*) "C", pos(i,:)
      !!write(12,*) "S", pos(i+1,:)
   !!enddo

   
   !!write(12,*) "O", pos(2*n-1,:)
  !! write(12,*) "H", pos(2*n,:)
  !! write(12,*) "N", pos(2*n+1,:)
  !------End of Chunk of code i needed initially to generate the lattice for the particles in the system--------!!!! 
    call distance(2*n-1,2*n,rijsq,vec)
    r=sqrt(rijsq)
    charge(2*n-1)=rdependent_charge(r,1)
    charge(2*n)  =rdependent_charge(r,2)
    charge(2*n+1)=rdependent_charge(r,3)
    write(123,*) r,charge(2*n-1), charge(2*n), charge(2*n+1)


    do i=1,2*n-2
      do j=2*n-1,2*n+1
        call distance(i,j,rijsq,vec)
        rij=sqrt(rijsq)
        if(rij<5.d0)write(800,*) tim,i,j,rij
      enddo
    enddo

  end subroutine generate_system_lattice
 !-----------------------------------------------------------!
   subroutine evolve
    implicit none
    integer i
   ! call check_acceleration
    open(16,file='pot_kin.out')
    open(17,file="energy.out")
    open(18,file="vmd.xyz")
    do i=1,nsteps
      call write_output
      if(mod(i,18000)==0)vel=0.d0   !used for obtaining intial configuration to perform MD on
      call write_vmd
      call cpu_time(tic2)
      call evolve_1step
      call cpu_time(toc2)
      write(741,*) 'time take for double do loop', toc2-tic1
      call restrict_atoms_inside_the_box
      call velocity_bin
      call boltzmann
    end do
    !!open(200,file='initial-pos.xyz')
    !write(200,*) pos
    !write(200,*)

    open(590,file='pos20.xyz')
    write(590,*) pos
  end subroutine evolve 
 !-----------------------------------------------------------!
  subroutine evolve_1step
    !! Evolution by one timestep dt using Velocity-Verlet method
    !! Velocity-Verlet - https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet
    implicit none
    double precision vec_rij(3)
    
    pos_old=pos !r(t) stored pos as they are required in Rattle_positions
    pos=pos+vel*dt+0.5*acc*dt*dt !r(t+dt)
    vel=vel+0.5*acc*dt
    vel_old=vel !needed in rattle subroutine
    pos_new=pos
    call rattle_position(pos_old,vel_old,pos_new,pos,vel)
   !call check_compute_pot(pos,pot,acc)         ! getting updated acceleration and potential
    call compute_pot(pos,pot,acc)         ! getting updated acceleration and potential
    vel=vel+0.5*acc*dt
    vel_old=vel
    call rattle_velocity(vel_old,pos,vel) 
    
    call compute_energy
    tim=tim+dt                       ! Updating current time
   
  end subroutine evolve_1step
 !-----------------------------------------------------------!
  subroutine rattle_position(pos_old,vel_old,pos_new,pos,vel)
    implicit none
    integer i,j,k
    double precision, intent(in):: pos_old(2*n+1,3),vel_old(2*n+1,3),pos_new(2*n+1,3)
    double precision, intent(out):: pos(2*n+1,3),vel(2*n+1,3)
    double precision lambda_r
    double precision vec_rij(3),vec_rij_prime(3), bl
    double precision a1,b2,c3, mu
    bl=bond_length_CH3Cl
    do i=1,2*n-2,2
      j=i+1

      do k=1,3 
        vec_rij(k)=pos_old(i,k)-pos_old(j,k)
        if(vec_rij(k)>box_len/2.d0)vec_rij(k)=vec_rij(k)-box_len
        if(vec_rij(k)<-box_len/2.d0)vec_rij(k)=vec_rij(k)+box_len
        vec_rij_prime(k)=pos_new(i,k)-pos_new(j,k) 
        if(vec_rij_prime(k)>box_len/2.d0)vec_rij_prime(k)=vec_rij_prime(k)-box_len
        if(vec_rij_prime(k)<-box_len/2.d0)vec_rij_prime(k)=vec_rij_prime(k)+box_len
      enddo
      mu=(1.d0/mass(i))+(1.d0/mass(j))
      a1=0.25*dt**4 *bl*bl*mu**2
      b2=sum(vec_rij * vec_rij_prime) * dt**2 * mu
      c3=sum(vec_rij_prime*vec_rij_prime)-bl**2
      lambda_r=(0.5/a1) * (-b2 + sqrt(b2*b2 - 4*a1*c3))
      pos(i,:)=pos_new(i,:)+0.5d0*lambda_r*vec_rij*dt*dt/mass(i)
      pos(j,:)=pos_new(j,:)-0.5d0*lambda_r*vec_rij*dt*dt/mass(j)
      vel(i,:)=vel_old(i,:)+0.5d0*dt*lambda_r*vec_rij/mass(i)
      vel(j,:)=vel_old(j,:)-0.5d0*dt*lambda_r*vec_rij/mass(j) 
   enddo
  ! stop
  end subroutine rattle_position
 !-----------------------------------------------------------!
  subroutine rattle_velocity(vel_old,pos,vel)
    implicit none
    integer i,j,k
    double precision, intent(in)::vel_old(2*n+1,3),pos(2*n+1,3)
    double precision, intent(out)::vel(2*n+1,3)
    double precision vec_rij(3), vel_ij_prime(3),lambda_v,bl,mu
    bl=bond_length_CH3Cl
    do i=1,2*n-2,2
      j=i+1
      !vel=vel_old+ 0.5*acc*dt    !v(t+dt)
      mu=(1.d0/mass(i))+(1.d0/mass(j))
      do k=1,3  
        vec_rij(k)=pos(i,k)-pos(j,k)
        if(vec_rij(k)>box_len/2.d0)vec_rij(k)=vec_rij(k)-box_len
        if(vec_rij(k)<-box_len/2.d0)vec_rij(k)=vec_rij(k)+box_len
        vel_ij_prime(k)=vel_old(i,k)-vel_old(j,k)
     enddo
     lambda_v=-sum(vec_rij*vel_ij_prime)/(0.5*dt*mu*bl*bl)
     vel(i,:)=vel_old(i,:)+0.5*dt*lambda_v*vec_rij/mass(i)
     vel(j,:)=vel_old(j,:)-0.5*dt*lambda_v*vec_rij/mass(j)
   enddo
  end subroutine rattle_velocity
 !-----------------------------------------------------------!
  subroutine write_output
   implicit none

   write(16,*) tim,kin_energy, pot, Average_temperature
   write(17,*) tim*0.02418d0, energy*2.19474d+5
  end subroutine write_output
 !-----------------------------------------------------------!
  subroutine write_vmd
   implicit none
   integer i
   write(18,*)  2*(n)+1
   write(18,*)
   do i=1,2*(n-1),2
     write(18,*) "C", pos(i,:)
     write(18,*) "S", pos(i+1,:)
   enddo
   
   write(18,*) "O", pos(2*n-1,:)
   write(18,*) "H", pos(2*n,:)
   write(18,*) "N", pos(2*n+1,:)

  end subroutine write_vmd
 !-----------------------------------------------------------! 
  subroutine compute_pot(pos,pot,acc)
    implicit none
    
    integer i,j,k,l,m
    double precision rijsq, vec(3), rij6,rij12, rij
    double precision, intent(in)::pos(2*n+1,3)
    double precision, intent(out)::acc(2*n+1,3), pot
    double precision pot_ss, pot_cs,pot_ab,dpot_drij_ss, dpot_drij_cs, dpot_drij
    double precision LJ, coul, dLJ_drij, dcoul_drij
    double precision LJ_sol, dLJ_drij_sol
    double precision vec_AB(3), vec_AH(3), dpot_dRAB, dpot_dr
    double precision pot_rc_ss
    double precision tic1,toc1
    pot=0.d0
    acc=0.d0
    pot_ss=0.d0
    pot_cs=0.d0
    pot_ab=0.d0
    dpot_drij=0.d0
    do i=1, 2*n-2   !loop over all solvent molecules
      do j=i+1, 2*n+1
        if(j.ne.neighbour(i)) then   !excludes interaction between atoms of same solvent molecule
          call distance(i,j,rijsq,vec)
          if(rijsq<rcsq) then
            rij6=rijsq**3
            rij12=rij6*rij6 
            rij=sqrt(rijsq)
            !if(i==1.and.j==2) write(42,*) rij
            write(150,*)i,j,rij,rc
            if(j.ne.2*n-1.and.j.ne.2*n.and.j.ne.2*n+1)then     ! Solvent-Solvent interactions!! if statment- excluding complex molecules
              pot_rc_ss= A(i)*A(j)/rc**12 - C(i)*C(j)/rc**6            !cut off for LJ pot solvent-solvent interactions
              LJ_sol= A(i)*A(j)/rij12 - C(i)*C(j)/rij6-pot_rc_ss
              dLJ_drij_sol=(1.d0/rijsq) * (-12.d0*A(i)*A(j)/rij12 + 6.d0*C(i)*C(j)/rij6)
            
              call evaluate_coulomb(i,j,rij,rijsq,rc,rt,coul,dcoul_drij)
            else
              LJ_sol=0.d0
              dLJ_drij_sol=0.d0
              coul=0.d0
              dcoul_drij=0.d0
            endif  
            pot_ss=LJ_sol+coul
            dpot_drij_ss=dLJ_drij_sol+dcoul_drij
            if(j==2*n-1.or.j==2*n+1)then             !! Complex-solvent Lennard-Jones interactions
              call evaluate_LJ_pot(rijsq,rij6,rij12,LJ,dLJ_drij)
              write(68,*) LJ, dLJ_drij
            else 
              LJ=0.d0
              dLJ_drij=0.d0         
            endif
            if(j==2*n.or.j==2*n-1.or.j==2*n+1)then
              call evaluate_coulomb(i,j,rij,rijsq,rc,rt,coul,dcoul_drij)        !Complex-solvent coulomb interactions
            else
              coul=0.d0
              dcoul_drij=0.d0
            endif
            pot_cs=LJ+coul
            dpot_drij_cs=dLJ_drij+dcoul_drij
        
            pot=pot+pot_ss+pot_cs
            dpot_drij=dpot_drij_ss+dpot_drij_cs
              !   if(i==2*n-1.or.j==2*n-1)write(40,*) i,j,dpot_drij
            acc(i,:)=acc(i,:)-1.d0/mass(i) * (dpot_drij*vec)
            acc(j,:)=acc(j,:)+1.d0/mass(j) * (dpot_drij*vec)
          endif
        endif
      enddo
    enddo
   
    call distance(3,4,rijsq,vec)
    rij=sqrt(rijsq)
    write(42,*) rij
           
    call distance(2*n-1,2*n+1,rijsq,vec)
    RAB=sqrt(rijsq)
    vec_AB=vec
    call distance(2*n-1,2*n,rijsq,vec)
    r=sqrt(rijsq)
    vec_AH=vec

    call  compute_potab(RAB,r,vec_AH,vec_AB,pot_ab, dpot_dRAB,dpot_dr)
    pot=pot+pot_ab
    write(1000,*)RAB,r, pot_ab, dpot_dRAB, dpot_dr
    write(1000,*) vec_AH, vec_AB
    write(1000,*)
    acc(2*n-1,:)=acc(2*n-1,:)-(1.d0/mass(2*n-1) * (dpot_dRAB * (vec_AB/RAB) + dpot_dr * (vec_AH/r))) 
    acc(2*n,:)=acc(2*n,:)- (1.d0/mass(2*n) *  (dpot_dRAB * 0.0 + dpot_dr * (-vec_AH/r))) 
    acc(2*n+1,:)=acc(2*n+1,:)-(1.d0/mass(2*n+1) * (dpot_dRAB * (-vec_AB/RAB) + dpot_dr * 0.0))
    write(502,*) tim,RAB, r

!    stop
  end subroutine compute_pot
 !-----------------------------------------------------------!
  subroutine check_compute_pot(pos,pot,acc)
    implicit none
    integer i,j,k,l,m
    double precision rijsq, vec(3), rij6,rij12, rij
    double precision, intent(in)::pos(2*n+1,3)
    double precision, intent(out)::acc(2*n+1,3), pot
    double precision pot_ss, pot_cs,pot_ab,dpot_drij_ss, dpot_drij_cs, dpot_drij
    double precision LJ, coul, dLJ_drij, dcoul_drij
    double precision LJ_sol, dLJ_drij_sol
    double precision vec_AB(3), vec_AH(3), dpot_dRAB, dpot_dr
    double precision pot_rc_ss
    pot=0.d0
    acc=0.d0
    pot_ss=0.d0
    pot_cs=0.d0
    pot_ab=0.d0
    dpot_drij=0.d0
    do i=1, 2*n-2   !loop over all solvent molecules
      do j=2*n-1, 2*n+1
        call distance(i,j,rijsq,vec)
        if(rijsq<rcsq) then
          rij6=rijsq**3
          rij12=rij6*rij6 
          rij=sqrt(rijsq)
          if(j==2*n-1.or.j==2*n+1)then             ! Complex-solvent Lennard-Jones interactions
            call evaluate_LJ_pot(rijsq,rij6,rij12,LJ,dLJ_drij)
            write(68,*) LJ, dLJ_drij
          else 
            LJ=0.d0
            dLJ_drij=0.d0         
          endif
          pot_cs=LJ
          dpot_drij_cs=dLJ_drij
          write(560,*)tim, i,j,pot_cs      
          write(560,*)  
          pot=pot+pot_cs
          dpot_drij=dpot_drij_cs
          !if(j==2*n) write(345,*) i,dpot_drij
          acc(i,:)=acc(i,:)-1.d0/mass(i) * (dpot_drij*vec)
          acc(j,:)=acc(j,:)+1.d0/mass(j) * (dpot_drij*vec)
          !if(j==2*n) write(345,*) i,dpot_drij, acc(j,:)
        endif
      enddo
    enddo
    call distance(3,4,rijsq,vec)
    rij=sqrt(rijsq)
    write(42,*) rij

  end subroutine check_compute_pot
 !-----------------------------------------------------------!
  subroutine VLJ
    implicit none
    integer i
    double precision rij,rijsq,rij6,rij12,LJ,dLJ_drij
    do i=1,100
      rij=5.0+i*10/100.d0
      rijsq=rij*rij
      rij6=rij**6
      rij12=rij**12
      call  evaluate_LJ_pot(rijsq,rij6,rij12,LJ,dLJ_drij)
      write(789,*) rij, LJ
   enddo
   !stop
  end subroutine VLJ
 !-----------------------------------------------------------!

  subroutine RAB_r       !plot of pot_ab versus r at R=2.7 Angstorms
    implicit none
    integer i,j
    double precision RAB,r
    double precision pot_ab, dpot_dRAB, dpot_dr
    double precision term_1, term_2, term_3
    double precision dterm1_RAB, dterm2_RAB, dterm3_RAB

    RAB=2.7d0*1.889    
    do i=1,100
      r=0.8+i*1.d0/100.d0
      r=r*1.889d0
      term_1=bb*exp(-aa*RAB)
      term_2=exp(-n_A*(r-d_A)**2/(2.d0*r))
      term_3=exp((-n_B*(RAB-r-d_B)**2)/(2.d0*(RAB-r)))
      pot_ab= term_1+DA*(1-term_2)+cc*DA*(1-term_3)
      write(600,*) tim,RAB,r,pot_ab
      dpot_dRAB=-aa*term_1+ cc*DA * (((-n_B*(RAB-r-d_B)*(RAB-r+d_B))/(2.d0*(RAB-r)**2)) * (-term_3))
      dpot_dr=DA*(((-n_A*(r-d_A)*(r+d_A))/(2.d0*r**2))*(-term_2))+cc*DA*(((n_B*(RAB-r-d_B)*(RAB-r+d_B))/(2.d0*(RAB-r)**2))*(-term_3))
      write(1200,*)RAB,r, pot_ab, dpot_dRAB, dpot_dr
    enddo
   stop
  end subroutine RAB_r
 !-----------------------------------------------------------!
  subroutine evaluate_LJ_pot(rijsq,rij6,rij12,LJ,dLJ_drij)
    implicit none
    double precision, intent(in) ::rijsq,rij6,rij12
    double precision, intent(out):: LJ, dLJ_drij
    pot_rc_cs=4.d0*epsilon*((sigma12/rc**12)-(sigma6/rc**6))
    LJ=4.d0*epsilon*(sigma12/rij12-sigma6/rij6)-pot_rc_cs
    dLJ_drij=(24.d0*epsilon/rijsq)*((sigma6/rij6)-(2.d0*sigma12/rij12))
    write(67,*) LJ, dLJ_drij
  end subroutine evaluate_LJ_pot
 !-----------------------------------------------------------!
  subroutine evaluate_coulomb(i,j,rij,rijsq,rc,rt,coul,dcoul_drij)
    implicit none
    integer,intent(in):: i,j
    double precision, intent(in):: rij, rijsq
    double precision, intent(out):: coul, dcoul_drij
    double precision T_rij, rt, rc, vec(3)
    double precision factor1,factor2,factor3,coul_prime

    charge(2*n-1)=rdependent_charge(r,1)
    charge(2*n)  =rdependent_charge(r,2)
    charge(2*n+1)=rdependent_charge(r,3)
    write(125,*) tim,r, charge(2*n-1), charge(2*n), charge(2*n+1)
    if(rij.le.rt)then
      coul=(charge(i)*charge(j))/rij
      T_rij=1.d0
      coul=coul*T_rij
      dcoul_drij=(-charge(i)*charge(j))/rijsq
      dcoul_drij=(1.d0/rij)*dcoul_drij
    endif

    if(rt.le.rij.and.rij.le.rc)then
      coul_prime=(charge(i)*charge(j))/rij
      T_rij=1.d0-((rc-rt)**(-3)*(rij-rt)**2*(3.d0*rc-rt-2.d0*rij))
      coul=coul_prime*T_rij
      dcoul_drij=(1.d0/rij)*(((-charge(i)*charge(j)/rijsq)*T_rij)-(coul_prime*(rc-rt)**(-3)*(-6*rijsq+6*rt*rij+6*rc*rij-6*rc*rt)))
    endif
   
    if(rc.le.rij)then
      coul=charge(i)*charge(j)/rij
      T_rij=0.d0
      coul=coul*T_rij
      dcoul_drij=0.d0
    endif
  end subroutine evaluate_coulomb
 !-----------------------------------------------------------!
  subroutine compute_potab(RAB,r,vec_AH,vec_AB,pot_ab, dpot_dRAB,dpot_dr)
    implicit none
    integer i,j
    double precision, intent(in)::RAB,r,vec_AH(3),vec_AB(3)
    double precision rijsq, vec(3)
    double precision, intent(out) :: pot_ab, dpot_dRAB, dpot_dr
    double precision term_1, term_2, term_3
    double precision dterm1_RAB, dterm2_RAB, dterm3_RAB
    
    term_1=bb*exp(-aa*RAB)
    term_2=exp(-n_A*(r-d_A)**2/(2.d0*r))
    term_3=exp((-n_B*(RAB-r-d_B)**2)/(2.d0*(RAB-r)))
    pot_ab=term_1+ DA*(1-term_2)+cc*DA*(1-term_3)
    write(500,*) tim,RAB,r
    dpot_dRAB=-aa*term_1+ cc*DA * (((-n_B*(RAB-r-d_B)*(RAB-r+d_B))/(2.d0*(RAB-r)**2)) * (-term_3))
    dpot_dr=DA*(((-n_A*(r-d_A)*(r+d_A))/(2.d0*r**2))*(-term_2))+cc*DA*(((n_B*(RAB-r-d_B)*(RAB-r+d_B))/(2.d0*(RAB-r)**2))*(-term_3))
  !! write(1000,*)RAB,r, pot_ab, dpot_dRAB, dpot_dr
  ! stop
  end subroutine compute_potab
 !-----------------------------------------------------------!
  subroutine distance(i,j,rijsq,vec)
    implicit none
    integer k
    integer, intent(in):: i,j
    double precision, intent(out):: rijsq, vec(3)
    rijsq=0.d0
    do k=1,3
      vec(k)=pos(i,k)-pos(j,k)
      if(vec(k)>box_len/2.d0)vec(k)=vec(k)-box_len
      if(vec(k)<-box_len/2.d0)vec(k)=vec(k)+box_len
      rijsq=rijsq+vec(k)**2
    enddo
  end subroutine distance
 !-----------------------------------------------------------!
  subroutine restrict_atoms_inside_the_box
    implicit none
    integer i,j
    do i=1,2*n+1
      do j=1,3
        if(pos(i,j)>box_len)pos(i,j)=pos(i,j)-box_len
        if(pos(i,j)<0.d0)pos(i,j)=pos(i,j)+box_len
      enddo
    enddo
  end subroutine restrict_atoms_inside_the_box
 !-----------------------------------------------------------!
  function rdependent_charge(r,i)
   implicit none
   integer, intent(in)::i
   double precision,intent(in) ::r
   double precision rdependent_charge
   double precision num,den, func_r
   num=r-r0
   den=sqrt((r-r0)**2+l*l)
   func_r=0.5d0*(1.d0+(num/den))
   rdependent_charge=(1.d0-func_r)*ee(i,1)+(func_r*ee(i,2))
   write(66,*) "func," ,func_r,1.d0-func_r, r, ee(i,1), ee(i,2) 
   write(66,*) "charge," ,rdependent_charge 
   write(66,*)
 end function
 !-----------------------------------------------------------!
  subroutine compute_energy
    implicit none
    integer i

    kin_energy=0.d0

    do i=1,2*n+1
      kin_energy=kin_energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
    enddo

    energy=kin_energy+pot

   ! Average_temperature= 2.d0/3.d0*kin_energy/(kb*real(n))

  end subroutine compute_energy
 !-----------------------------------------------------------!
  subroutine thermal_velocities(temperature)
    implicit none
    integer i,j
    double precision, intent(in):: temperature
    double precision  rnd

    do i=1,n
      do j=1,3
        call gaussian_random_number_generator(rnd)
        vel(i,j)=rnd*sqrt(kb*temperature/mass(i))
      enddo
   enddo

   do j=1,3
     vel(:,j)=vel(:,j)-sum(mass*vel(:,j))/sum(mass)
   enddo
  call compute_energy
  end subroutine thermal_velocities
 !-----------------------------------------------------------!
  subroutine gaussian_random_number_generator(rnd)
    implicit none

    double precision U1,U2
    double precision, intent(out):: rnd
    call random_number(U1)
    call random_number(U2)
    rnd= dsqrt(-2.d0*log(U1))*cos(2.d0*pi*U2)

   end subroutine gaussian_random_number_generator
 !----------------------------------------------------------------!
  subroutine check_acceleration
    implicit none
    integer i,j,k
    real*8 delx
    real*8 pot_sav,acc_sav_i(2*n+1,3),acc_sav_j(2*n-1,3),acc(2*n+1,3)
    double precision rijsq, rc, dLJ_drij, LJ,rij6,rij12, acc_1, acc_2,vec(3),rij, rij1,rij2    
    double precision pot_ab,dpot_dRAB,dpot_dr,vec_AB(3),vec_AH(3)
    do i=1,2*n-2
      do j=2*n-1,2*n+1
        do k=1,3
!    i=2
!    j=2*n-1
!          call initial_conditions    
          !call compute_pot(pos,pot,acc)
          call check_compute_pot(pos,pot,acc)
         ! call compute_potab(pot_ab, dpot_dRAB, dpot_dr,vec_AB, vec_AH)
         ! pot_sav=pot_ab
          pot_sav=pot
          acc_sav_i(i,k)=acc(i,k)
          !acc_sav_i(i,k)=dpot_dr
          !acc_sav_i(i,k)=dpot_dRAB
          !rij6=r
          ! rij6=RAB
          !write(6,*) "befor-sft",RAB,r
          delx=1.d-5
          pos(i,k)=pos(i,k)+delx
          !call compute_pot(pos,pot,acc)
          call check_compute_pot(pos,pot,acc)
          !call compute_potab(pot_ab, dpot_dRAB, dpot_dr,vec_AB, vec_AH)
          !rij12=r
          !rij12=RAB
          acc(i,k)=-1.d0/mass(i)*((pot-pot_sav)/delx)
          !acc(i,k)=((pot_ab-pot_sav)/(rij12-rij6))
          !write(6,*) "check_acc",RAB,r
          write(250,*) i,j,k,acc_sav_i(i,k),acc(i,k), acc(i,k)-acc_sav_i(i,k)
          pos(i,k)=pos(i,k)-delx   !changing position back in order to displace the next atom 
        enddo
      enddo
    enddo
 ! stop
  end subroutine check_acceleration
 !----------------------------------------------------------------!
  subroutine boltzmann
    implicit none
    integer i
    double precision boltz(2*n-2)
    !open(45, file='boltzmann_distribution.out')
    !do i=1,2*n+1
    do i=1,2*n-2
      if(mod(i,2)==0)then
        boltz=boltz+exp((-mass(i)*vel(i,1)*vel(i,1))/(2.d0*kb*temperature))
        !write(45,*) vel(i,1), boltz
      endif
      boltz=0.d0
   enddo
  end subroutine
 !----------------------------------------------------------------!
  subroutine write_boltz
   implicit none
   integer i
   double precision boltz(2*n-2)
   open(45, file='boltzmann_distribution.out')
   boltz=boltz/nsteps
   do i=1,2*n-2
     if(mod(i,2)==0)then
       write(45,*) vel(i,1), boltz
     endif
   enddo
  end subroutine
 !----------------------------------------------------------------!
  subroutine velocity_bin
    implicit none
    integer i,vi
    double precision v_diff
    !do i=1,2*n+1
    do i=1,2*n-2
      if(mod(i,2)==0)then
        v_diff=vel(i,1)-(-0.01d0)
        vi=nint(v_diff/delv)
        if(vi>0.and.vi<nbins)then
          vdist(vi)=vdist(vi)+1
        endif
      endif
    enddo
   end subroutine
  !-----------------------------------------------------------------
   subroutine distribution
     implicit none
     integer i
     double precision total,v
     vdist=vdist/nsteps
     vdist=vdist/sum(vdist)  !normalizing
     do i=1,nbins
       v=-0.01+delv*(i+0.5)
       write(423,*) v, vdist(i)
     enddo
   end subroutine distribution
  !------------------------------------------------------------------
end program alpha_fcc
