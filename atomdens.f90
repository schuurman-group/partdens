!
! Build the atomic density using sums of hydrogenic wavefunctions with
! and integrated screening parameter
!
module atomdens
  use accuracy

contains

!
! Build the electron density from the sum of the squares of
! the atomic wavefunction
!
subroutine psisq(r, nr, Z, chg, norb, niter, thrsh, psi2)
  integer(ik) :: nr, Z, chg, norb, niter
  real(rk)    :: thrsh, r(nr)
  real(ark)   :: psi(norb, nr), psi2(nr)

  call psi_integ(r, nr, Z, chg, norb, niter, thrsh, psi)
  psi2 = sum(psi**2, dim=1)
end subroutine psisq

!
! Build the components of an atomic wavefunction using the integrated
! Slater's rules for screening S_n(r)
!
subroutine psi_integ(r, nr, Z, chg, norb, niter, thrsh, psi)
  integer(ik) :: i, j, nr, Z, chg, norb, ne, niter
  integer(ik) :: nnum(norb), lnum(norb), nocc(norb)
  real(rk)    :: dr, thrsh, selfsc, norm, cj
  real(rk)    :: r(nr), zeff(nr), last(nr), resid(nr), fr(nr)
  real(rk)    :: scr(norb, nr)
  real(ark)   :: psi(norb,nr)

  ne = Z - chg
  call get_occ(ne, norb, nnum, lnum, nocc)
  dr = r(2) - r(1)
  call psi_slater(r, nr, Z, chg, norb, psi)
  last = r**2*sum(psi**2, dim=1)
  do i = 1,niter
    do j = 1,norb
      scr(j,:) = cumsum(r**2*psi(j,:)**2, nr)*dr
    end do
    do j = 1,norb
      selfsc = 0.5*(nocc(j) - 1) / nocc(j)
      zeff = Z - sum(scr(1:j-1,:), dim=1) - selfsc*scr(j,:)
      fr = r_nl(r, nr, zeff, nnum(j), lnum(j))
      norm = sum(r**2*fr**2)*dr
      cj = sqrt(nocc(j) / norm)
      psi(j,:) = cj*fr
    end do
    resid = abs(r**2*sum(psi**2, dim=1) - last)
    if (maxval(resid) < thrsh) then
      return
    end if
    last = r**2*sum(psi**2, dim=1)
  end do
end subroutine psi_integ

!
! Build the components of an atomic wavefunction using Slater's
! rules for screen S_n
!
subroutine psi_slater(r, nr, Z, chg, norb, psi)
  integer(ik) :: i, nr, Z, chg, norb, ne, nnum(norb), lnum(norb), nocc(norb)
  real(rk)    :: ci, scr(norb), r(nr), zeff(nr)
  real(ark)   :: psi(norb,nr)

  ne = Z - chg
  call get_occ(ne, norb, nnum, lnum, nocc)
  scr = scrni(nocc, norb)
  do i = 1,norb
    zeff(:) = Z - scr(i)
    ci = sqrt(real(nocc(i)))
    psi(i,:) = ci*r_nl(r, nr, zeff, nnum(i), lnum(i))
  end do
end subroutine psi_slater

!
! Get the radial hydrogenic wavefunction for an n, l pair
!
function r_nl(r, nr, Z, n, l)
  integer(ik) :: i, j, nr, n, l
  real(rk)    :: ak(n-l), ak2(2*(n-l)-1), pwr(2*(n-l)-1)
  real(rk)    :: r(nr), Z(nr), rho(nr), aterm(nr), norm(nr)
  real(ark)   :: r_nl(nr)

  rho = 2*Z*r/n
  ! get the expansion coefficients
  ak(1) = 1
  aterm = 1
  do i = 1,n-l-1
    ak(i+1) = ak(i)*(i + l - n)/(i*(i + 2*l + 1))
    aterm = aterm + ak(i+1)*rho**i
  end do
  ! get the square of the expansion coefficients
  ak2(:) = 0
  do i = 1,n-l
    do j = 1,n-l
      ak2(i+j-1) = ak2(i+j-1) + ak(i)*ak(j)
    end do
  end do
  ! get the factorials of the powers of r
  pwr = (/(i, i=0,2*(n-l)-2)/) + 2*(1 + l)
  do i = 2,2*l+1
    pwr(1) = pwr(1)*i
  end do
  do i = 2,2*(n-l)-1
    pwr(i) = pwr(i)*pwr(i-1)
  end do
  ! get the normalization constant
  norm = sqrt((2*Z/n)**3/sum(pwr*ak2))
  ! assemble parts
  r_nl = norm*rho**l*aterm*exp(-rho/2)
end function r_nl

!
! Get a list of integer screening constants using Slater's rules
!
function scrni(nocc, norb)
  integer(ik) :: i, j, norb, nocc(norb)
  real(rk)    :: scrni(norb)

  do i = 1,norb
    scrni(i) = 0.5*(nocc(i) - 1)
    do j = 1,i-1
      scrni(i) = scrni(i) + nocc(j)
    end do
  end do
end function scrni

!
! Get the number of occupied or partially occupied orbitals by the
! number of electrons
!
function get_norb(ne)
  integer(ik) :: l, ne, maxocc, nl, get_norb

  get_norb = 0
  maxocc = 0
  nl = 0
  do
    nl = nl + 1
    do l = (nl - 1)/2,0,-1
      get_norb = get_norb + 1
      maxocc = maxocc + 4*l + 2
      if (maxocc >= ne) then
          return
      end if
    end do
  end do
end function get_norb

!
! Get the orbital occupancy and quantum numbers by the number
! electrons
!
subroutine get_occ(ne, norb, nnum, lnum, nocc)
  integer(ik) :: i, l, nin, maxocc ,ne, norb, nl
  integer(ik) :: nnum(norb), lnum(norb), nocc(norb)

  nl = 0
  i = 0
  maxocc = 0
  do
    nl = nl + 1
    do l = (nl - 1)/2,0,-1
      i = i + 1
      nin = 4*l + 2
      maxocc = maxocc + nin
      nnum(i) = nl - l
      lnum(i) = l
      if (maxocc >= ne) then
        nocc(i) = nin + ne - maxocc
        return
      else
        nocc(i) = nin
      end if
    end do
  end do
end subroutine get_occ

!
! Generate an equally spaced list of n elements between specified numbers
!
function linspace(first, last, n)
  integer(ik) :: i, n
  real(rk)    :: first, last, step, list(n), linspace(n)

  list = (/(i, i=0,n-1)/)
  step = (last - first) / (n - 1)
  linspace = list*step + first
end function linspace

!
! Return the cumulative sum of a 1D array
!
function cumsum(arr, nx)
  integer(ik) :: i, nx
  real(rk)    :: arr(nx), cumsum(nx)

  cumsum(1) = arr(1)
  do i = 2,nx
    cumsum(i) = arr(i) + cumsum(i-1)
  end do
end function cumsum

end module atomdens
