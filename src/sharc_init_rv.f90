program init_rv
implicit none
integer nat
real*8, allocatable :: geom(:,:),vel(:,:),noat(:),mass(:)
character*2, allocatable :: sat(:)
real*8 Esamp
integer traj0,trajf
integer nseed
integer, allocatable :: seed(:)

real*8, allocatable :: prob(:)
real*8 sumprob
real*8 Ek,mom(3),cm(3),masscm

integer i,j,k
integer itraj

write(6,*) "Number of atoms, and kinetic energy to be sampled (in eV)"
read(5,*) nat,Esamp
Esamp=Esamp/27.211
write(6,*) "Sampling kinetic energy (H,eV)",Esamp,Esamp*27.211
open(1,file="geom",status="old",iostat=i)
if (i.ne.0) stop "File geom does not exist"
allocate(sat(nat),noat(nat),mass(nat),geom(nat,3),vel(nat,3))
do i=1,nat
 read(1,*) sat(i),noat(i),(geom(i,j),j=1,3),mass(i)
enddo
close(1)
write(6,*) "Equilibrium geometry ready"

call random_seed(size=nseed)
allocate (seed(nseed))
write(6,*) "First trajectory, last and seed"
read(5,*) traj0,trajf,seed(1)

if (nseed.ne.1) then
 do i=2,nseed
  seed(i)=seed(i-1)+1
 enddo
endif

allocate(prob(3*nat))
open(1,file="geoms.xyz")
open(2,file="geoms.init")
vel=0.
sumprob=0.
i=0
write(1,*) nat
write(1,903) sumprob,sumprob,sumprob,Esamp," Traj ",i
do i=1,nat
 write(1,901) sat(i),(geom(i,j)*.5292,j=1,3),(vel(i,j),j=1,3)
enddo
write(2,*) nat
do i=1,nat
 write(2,901) sat(i),noat(i),mass(i)
 mass(i)=mass(i)*1822.888
enddo

!!! Center of masses
masscm=0.
do i=1,nat
 masscm=masscm+mass(i)
enddo
do j=1,3
 cm(j)=0.
 do i=1,nat
  cm(j)=cm(j)+mass(i)*geom(i,j)
 enddo
 cm(j)=cm(j)/masscm
enddo

 

call random_seed(put=seed)
do itraj=traj0,trajf
 call random_number(prob)
 sumprob=0.
 do i=1,3*nat
  prob(i)=prob(i)-.5
  sumprob=sumprob+abs(prob(i))
 enddo
 do i=1,3*nat
  prob(i)=prob(i)/sumprob
 enddo
 sumprob=0.
 do i=1,3*nat
  sumprob=sumprob+abs(prob(i))
 enddo
 k=0
 do i=1,nat
  do j=1,3
   k=k+1
   vel(i,j)=abs(prob(k))  !!! Energy for the coordinate with no units
   vel(i,j)=sqrt(2.*vel(i,j)/mass(i))
   if (prob(k).lt.0) vel(i,j)=-vel(i,j)
  enddo
 enddo

!!! Removing velocity of the center of masses using the total momentum
 do j=1,3
  mom(j)=0.
  do i=1,nat
   mom(j)=mom(j)+mass(i)*vel(i,j)
  enddo
 enddo
 do i=1,nat
  do j=1,3
   vel(i,j)=vel(i,j)-mom(j)/masscm
  enddo
 enddo

!!! Kinetic energy
 Ek=0.
 do i=1,nat
  do j=1,3
   Ek=Ek+mass(i)*vel(i,j)*vel(i,j)*.5
  enddo
 enddo
!!! Rescaling the velocity
 do i=1,nat
  do j=1,3
   vel(i,j)=sqrt(Esamp/Ek)*vel(i,j)
  enddo
 enddo
!!! Kinetic energy
 Ek=0.
 do i=1,nat
  do j=1,3
   Ek=Ek+mass(i)*vel(i,j)*vel(i,j)*.5
  enddo
 enddo

 sumprob=0.
 write(1,*) nat
 write(1,903) Esamp,sumprob,Esamp,Esamp," Traj ",itraj
 write(2,"(I6.6,A,3(x,E20.10e3))") itraj," Traj ",sumprob,Esamp,Ek
 do i=1,nat
  write(1,901) sat(i),(geom(i,j)*.5292,j=1,3),(vel(i,j)*1000.,j=1,3)
  write(2,900) (geom(i,j),j=1,3)
 enddo
 do i=1,nat
  write(2,900) (vel(i,j),j=1,3)
 enddo
 
enddo

900 format (1000(x,E20.10e3))
901 format (A2,1000(x,E20.10e3))
903 format (4(x,F20.10),A,I6.6)
end
   
