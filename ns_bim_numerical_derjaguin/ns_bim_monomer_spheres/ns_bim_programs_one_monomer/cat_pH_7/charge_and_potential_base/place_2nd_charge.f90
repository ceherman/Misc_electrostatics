module place_2nd_charge

  USE Pre_Constants
  USE Pre_csvformat
  implicit NONE
contains


subroutine print_second_charge_coordinates(charge1_coords, distance, charge2_coords)
! This is just an example illustrating the use of the subroutine place_point.
! place_point has two inputs:  spherical coordinates for the first charge
! given as an array of 3 elements (radius, zenith, and azimuth, respectively)
! and a separation distance between the two charges. It returns both the
! spherical and Cartesian coordinates of the second charge.
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: i
  real(kind = dp), intent(in) :: charge1_coords(3), distance
  real(kind = dp), intent(out) :: charge2_coords(3)
  real(kind = dp) :: main_sph_coords(3), second_sph_coords(3)

  call cart_to_sph(charge1_coords, main_sph_coords)

!print*,"Please enter the spherical coordinates of the first charge as"
!print*,"radius, zenith (theta), azimuth (phi)"
!read*, main_sph_coords(1), main_sph_coords(2), main_sph_coords(3)
!print*,"Please enter the separation distance in Angstroms"
!read*, distance

call place_point(main_sph_coords, distance, second_sph_coords, charge2_coords)

!do i = 1,3
!  print*,charge2_coords(i)
!end do

end subroutine

subroutine sph_to_cart(sph_coords, cart_coords)
! Computes Cartesian from spherical coordinates
! Input:  sph_coords (a 3-element array containing, radius in Angstroms, zenith
!            angle (theta) in radians, and the azimuthal angle (phi) in radians,
!            respectively)
! Output:  cart_coords (a 3-element array containing x, y, and z in Angstroms,
!             respectively)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(kind = dp) :: r, theta, phi
  real(kind = dp), dimension(3), intent(in)  :: sph_coords
  real(kind = dp), dimension(3), intent(out) :: cart_coords

  r     = sph_coords(1)
  theta = sph_coords(2)
  phi   = sph_coords(3)

  cart_coords(1) = r*sin(theta)*cos(phi) ! x
  cart_coords(2) = r*sin(theta)*sin(phi) ! y
  cart_coords(3) = r*cos(theta)          ! z
end subroutine


subroutine cart_to_sph(cart_coords, sph_coords)
! Computes spherical from Cartesian coordinates
! Input:  cart_coords (a 3-element array containing x, y, and z in Angstroms,
!            respectively)
! Output:  sph_coords (a 3-element array containing, radius in Angstroms, zenith
!             angle (theta) in radians, and the azimuthal angle (phi) in radians,
!             respectively)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(kind = dp) :: x, y, z, r, theta, phi
  real(kind = dp), dimension(3), intent(in) :: cart_coords
  real(kind = dp), dimension(3), intent(out)  :: sph_coords

  x = cart_coords(1)
  y = cart_coords(2)
  z = cart_coords(3)

  r = sqrt(x**2 + y**2 + z**2)
  theta = acos(z/r)

  if (x > 0) then
      phi = atan(y/x)
      if (phi < 0) then
          phi = phi + 2*pai
      end if
  else if (x == 0) then
      if (y == 0) then
          phi = 0
      else if (y > 0) then
          phi = pai/2
      else
          phi = 3*pai/2
      end if
  else
      phi = atan(y/x) + pai
  end if

  sph_coords(1) = r
  sph_coords(2) = theta
  sph_coords(3) = phi
end subroutine


subroutine coplanar_vectors(sph_coords, b, c)
! Computes orthonormal vectors in the plane defined by the point given
! in sph_coords and the normal vector passing through this point and the origin
! Inputs: sph_coords (a 3-element array containing, radius in Angstroms, zenith
!            angle (theta) in radians, and the azimuthal angle (phi) in radians,
!            respectively)
! Outputs:  b and c (3-element arrays containing x, y, and z components (in
!             Angstroms) of orthonormal vectors in the plane described above)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  real(kind = dp) :: a1, a2, a3, b1, b2, b3, c1, c2, c3, denom_12, denom_13, &
                     denom_23, max_denom, magnitude_c
  real(kind = dp), dimension(3), intent(in)  :: sph_coords
  real(kind = dp), dimension(3)              :: cart_coords
  real(kind = dp), dimension(3), intent(out) :: b, c

  call sph_to_cart(sph_coords, cart_coords)
  a1 = cart_coords(1)
  a2 = cart_coords(2)
  a3 = cart_coords(3)

  denom_12 = sqrt(a1**2 + a2**2)
  denom_13 = sqrt(a1**2 + a3**2)
  denom_23 = sqrt(a2**2 + a3**2)
  max_denom = max(denom_12, denom_13, denom_23)

  if (denom_23 == max_denom) then
      b1 = 0.0
      b2 = -a3/denom_23
      b3 = a2/denom_23
  else if (denom_13 == max_denom) then
      b1 = a3/denom_13
      b2 = 0.0
      b3 = -a1/denom_13
  else
      b1 = -a2/denom_12
      b2 = a1/denom_12
      b3 = 0.0
  end if

  c1 =  a2*b3 - a3*b2
  c2 = -a1*b3 + a3*b1
  c3 =  a1*b2 - a2*b1

  magnitude_c = sqrt(c1**2 + c2**2 + c3**2)
  c1 = c1/magnitude_c
  c2 = c2/magnitude_c
  c3 = c3/magnitude_c

  b = [b1, b2, b3]
  c = [c1, c2, c3]
end subroutine


subroutine place_point(main_sph_coords, distance, second_sph_coords, second_cart_coords)
! Randomly places a charge on the circle of radius 'distance' in the plane
! defined by the point given in sph_coords and the normal vector passing through
! this point and the origin. The circle is centerd at the point given in sph_coords.
! Inputs: main_sph_coords (a 3-element array containing, radius in Angstroms, zenith
!            angle (theta) in radians, and the azimuthal angle (phi) in radians,
!            respectively, of the first charge)
!         distance (separation distance between the first and second point Charges
!            in Angstroms)
! Outputs:  second_sph_coords (a 3-element array containing, radius in Angstroms, zenith
!              angle (theta) in radians, and the azimuthal angle (phi) in radians,
!              respectively, of the second charge)
!           second_cart_coords (corresponding Cartesian coordinates of the second
!              charge, in Angstroms)
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer :: i
  real(kind = dp) :: t
  real(kind = dp), dimension(3) :: a, b, c
  real(kind = dp), dimension(3), intent(in)  :: main_sph_coords
  real(kind = dp), intent(in)                :: distance
  real(kind = dp), dimension(3), intent(out) :: second_sph_coords, second_cart_coords

  call sph_to_cart(main_sph_coords, a)
  call coplanar_vectors(main_sph_coords, b, c)

  call init_random_seed()
  call random_number(t)
  t = 2*pai*t

  second_cart_coords = [0, 0, 0]
  do i = 1, 3
      second_cart_coords(i) = a(i) + distance*cos(t)*b(i) + distance*sin(t)*c(i)
  end do

  call cart_to_sph(second_cart_coords, second_sph_coords)
end subroutine


subroutine init_random_seed()
  ! Used to make random number generation portable
  ! This code was obtained on 2019-11-24 from:
  ! https://stackoverflow.com/questions/37304793/random-number-generator-produces-same-sequence-even-though-its-seeded
   use iso_fortran_env, only: int64
   implicit none
   integer, allocatable :: seed(:)
   integer :: i, n, un, istat, dt(8), pid
   integer(int64) :: t

   call random_seed(size = n)
   allocate(seed(n))
   ! First try if the OS provides a random number generator
   open(newunit=un, file="/dev/urandom", access="stream", &
        form="unformatted", action="read", status="old", iostat=istat)
   if (istat == 0) then
      read(un) seed
      close(un)
   else
      ! Fallback to XOR:ing the current time and pid. The PID is
      ! useful in case one launches multiple instances of the same
      ! program in parallel.
      call system_clock(t)
      if (t == 0) then
         call date_and_time(values=dt)
         t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
              + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
              + dt(3) * 24_int64 * 60 * 60 * 1000 &
              + dt(5) * 60 * 60 * 1000 &
              + dt(6) * 60 * 1000 + dt(7) * 1000 &
              + dt(8)
      end if
      pid = getpid()
      t = ieor(t, int(pid, kind(t)))
      do i = 1, n
         seed(i) = lcg(t)
      end do
   end if
   call random_seed(put=seed)
 contains
   ! This simple PRNG might not be good enough for real work, but is
   ! sufficient for seeding a better PRNG.
   function lcg(s)
     integer :: lcg
     integer(int64) :: s
     if (s == 0) then
        s = 104729
     else
        s = mod(s, 4294967296_int64)
     end if
     s = mod(s * 279470273_int64, 4294967291_int64)
     lcg = int(mod(s, int(huge(0), int64)), kind(0))
   end function lcg
 end subroutine init_random_seed

end module
