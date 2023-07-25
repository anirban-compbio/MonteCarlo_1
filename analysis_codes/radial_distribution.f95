!*** Basic Radial Distribution Code; !*** read data from folded trajectory
!*** Written By : Anirban Paul 
!*** user can further generalize the code; !*** Code can be modified further...
!*** gfortran radial_distribution.f95 
!*** ./a.out > rdf.dat

program rdf
implicit none

real, parameter ::pi = acos(-1.0)
integer, parameter :: n=125
integer :: i,j,ii,jj, nf, maxbin, bin
real :: x(n),y(n),z(n), xd,yd,zd,dr, delr, Rmin,Rmax, nideal, dens
real, parameter :: L = 5.63
real, allocatable :: hist (:)
character :: dummy
character(2):: at(n)


hist = 0.0
dens = n/L**3
maxbin = 50
delr = L/maxbin
allocate (hist(maxbin))

open (unit=22,file="../trajectory_folded.xyz")
nf = 4000
do i = 1,nf
		read (22,*) dummy
		read (22,*) dummy
			
			do j = 1, n
			read	(22,*) at(j), x(j),y(j),z(j)
			!print *, at(j), x(j),y(j),z(j)
			end do

if (i > 300) then
	do ii = 1,n-1
		do jj = ii+1, n
				xd = x(ii)-x(jj)
				yd = y(ii)-y(jj)
				zd = z(ii)-z(jj)

				xd = xd - L*nint (xd/L)
				yd = yd - L*nint (yd/L)
				zd = zd - L*nint (zd/L)
				
				dr = sqrt (xd**2 + yd**2 + zd**2)
				bin = int (dr/delr)+1
				
				hist (bin) = hist(bin) + 2
				
				if (bin ==1) print*,ii,jj,dr
				
		end do
	end do
end if				
end do

hist = hist / (real(nf-300)*n)

	do i = 1,maxbin
  Rmin = (i-1)*delr
	Rmax = Rmin + delr
	nideal = dens*4*pi*(Rmax**3 - Rmin**3)/3
  hist(i) = hist (i)/nideal   
  print *, (i-1)*delr, hist(i)
 	end do  


end program
