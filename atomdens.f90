!
! Build the atomic density using sums of hydrogenic wavefunctions with
! an integrated form of Slater's rules
!
module atomdens
  use accuracy
  use atoms
  use math

contains

!
! Build the electron density from the sum of the squares of
! the atomic wavefunction
!
subroutine psisq(r, wr, nr, Z, chg, norb, niter, thrsh, psi2)
  integer(ik) :: nr, Z, chg, norb, niter
  real(rk)    :: thrsh, r(nr), wr(nr)
  real(ark)   :: psi2(nr)
  !
  real(ark)   :: psi(norb, nr)

  call psi_integ(r, wr, nr, Z, chg, norb, niter, thrsh, psi)
  psi2 = sum(psi**2, dim=1)
end subroutine psisq

!
! Build the components of an atomic wavefunction using the integrated
! Slater's rules for screening S_n(r)
! In a trivial even-spaced grid, weights are wr = dr*4*pi*r^2
!
subroutine psi_integ(r, wr, nr, Z, chg, norb, niter, thrsh, psi)
  integer(ik) :: nr, Z, chg, norb, niter
  real(rk)    :: thrsh, r(nr), wr(nr)
  real(ark)   :: psi(norb,nr), wpsi2j(nr)
  !
  integer(ik) :: i, j, ne, nnum(norb), lnum(norb), nocc(norb)
  real(rk)    :: norm, cj, zeff(nr)
  real(rk)    :: last(nr), resid(nr), fr(nr), scr(norb, nr), sfac(norb, nr)

  ne = Z - chg
  if (ne <= 0) then
    psi(:,:) = 0
    return
  end if
  call get_occ(ne, norb, nnum, lnum, nocc)
  ! initial guess
  call psi_slater(r, nr, Z, chg, norb, psi)
  last = wr*sum(psi**2, dim=1)
  do i = 1,niter
    ! set up integrated screening factor
    do j = 1,norb
      wpsi2j = wr*psi(j,:)**2
      scr(j,:) = cumsum(wpsi2j, nr) / nocc(j)
    end do
    sfac = scrni(nocc, nnum, lnum, norb, scr, nr)
    do j = 1,norb
      ! calculate Z_eff
      zeff = Z - sfac(j,:)
      ! find the single electron contribution and renormalize
      fr = r_nl(r, nr, zeff, nnum(j), lnum(j))
      norm = sum(wr*fr**2)
      cj = sqrt(nocc(j) / norm)
      ! add renormalized component to psi
      psi(j,:) = cj*fr
    end do
    ! check for convergence
    resid = abs(wr*sum(psi**2, dim=1) - last)
    if (maxval(resid) < thrsh) then
      return
    end if
    last = wr*sum(psi**2, dim=1)
  end do
end subroutine psi_integ

!
! Build the components of an atomic wavefunction using Slater's
! rules for screening S_n
!
subroutine psi_slater(r, nr, Z, chg, norb, psi)
  integer(ik) :: nr, Z, chg, norb
  real(rk)    :: r(nr)
  real(ark)   :: psi(norb,nr)
  !
  integer(ik) :: i, ne, nnum(norb), lnum(norb), nocc(norb)
  real(rk)    :: ci, scr(norb, nr), sfac(norb, nr), zeff(nr)

  ne = Z - chg
  if (ne <= 0) then
    psi(:,:) = 0
    return
  end if
  call get_occ(ne, norb, nnum, lnum, nocc)
  scr(:,:) = 1
  sfac = scrni(nocc, nnum, lnum, norb, scr, nr)
  do i = 1,norb
    zeff = Z - sfac(i,:)
    ci = sqrt(real(nocc(i)))
    psi(i,:) = ci*r_nl(r, nr, zeff, nnum(i), lnum(i))
  end do
end subroutine psi_slater

!
! Get the radial hydrogenic wavefunction for an n, l pair
!
function r_nl(r, nr, Z, n, l)
  integer(ik) :: nr, n, l
  real(rk)    :: r(nr), Z(nr)
  real(ark)   :: r_nl(nr)
  !
  integer(ik) :: i, j
  real(rk)    :: rho(nr), aterm(nr), norm(nr)
  real(rk)    :: ak(n-l), ak2(2*(n-l)-1), pwr(2*(n-l)-1)

  rho = 2*Z*r/n
  ! get the expansion coefficients
  ak(1) = 1
  aterm(:) = 1
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
  norm = sqrt((2*Z/n)**3/(4*pi*sum(pwr*ak2))) ! added in 4pi
  ! assemble parts
  r_nl = norm*rho**l*aterm*exp(-rho/2)
end function r_nl

!
! Get r-dependent screening parameters using Slater's rules
! See Slater, J. C. Phys. Rev. 1930, 36, 57-64.
!
function scrni(nocc, nnum, lnum, norb, scr, nr)
  integer(ik) :: norb, nr, nnum(norb), lnum(norb), nocc(norb)
  real(rk)    :: scrni(norb, nr), scr(norb, nr)
  !
  integer(ik) :: i, j

  do i = 1,norb
    if (nnum(i) == 1) then
      scrni(i,:) = 0.30*(nocc(i) - 1)*scr(i,:)
    else if (i < norb .and. lnum(i) == 0) then
      scrni(i,:) = 0.35*((nocc(i) - 1)*scr(i,:) + nocc(i+1)*scr(i+1,:))
    else
      scrni(i,:) = 0.35*(nocc(i) - 1)*scr(i,:)
    end if
    do j = 1,i-1
      if (lnum(i) < 2 .and. nnum(j) == nnum(i)) then
        scrni(i,:) = scrni(i,:) + 0.35*nocc(j)*scr(j,:)
      else if (lnum(i) < 2 .and. nnum(j) == nnum(i) - 1) then
        scrni(i,:) = scrni(i,:) + 0.85*nocc(j)*scr(j,:)
      else
        scrni(i,:) = scrni(i,:) + nocc(j)*scr(j,:)
      end if
    end do
  end do
end function scrni

!
! Get the number of occupied or partially occupied orbitals by the
! number of electrons
!
function get_norb(ne)
  integer(ik) :: ne, get_norb
  !
  integer(ik) :: l, maxocc, nl

  get_norb = 0
  maxocc = 0
  nl = 0
  step_nl: do
    nl = nl + 1
    do l = (nl - 1)/2,0,-1
      get_norb = get_norb + 1
      maxocc = maxocc + 4*l + 2
      if (maxocc >= ne) then
        return
      end if
    end do
  end do step_nl
end function get_norb

!
! Get the orbital occupancy and quantum numbers by the number of
! electrons
!
subroutine get_occ(ne, norb, nnum, lnum, nocc)
  integer(ik) :: ne, norb, nnum(norb), lnum(norb), nocc(norb)
  !
  integer(ik) :: i, j, l, nin, maxocc, nl, tmpn, tmpl, tmpo

  nl = 0
  i = 0
  maxocc = 0
  ! scan through increasing n+l starting from the highest l value
  do while (maxocc < ne)
    nl = nl + 1
    do l = (nl - 1)/2,0,-1
      i = i + 1
      nin = 4*l + 2
      maxocc = maxocc + nin
      nnum(i) = nl - l
      lnum(i) = l
      if (maxocc >= ne) then
        nocc(i) = nin + ne - maxocc
        exit
      else
        nocc(i) = nin
      end if
    end do
  end do
  ! sort the orbitals by n
  do i = 2,norb
    if (nnum(i) < nnum(i-1)) then
      tmpn = nnum(i)
      tmpl = lnum(i)
      tmpo = nocc(i)
      do j = i-1,2,-1
        if (nnum(j) == tmpn) then
          nnum(j+1) = tmpn
          lnum(j+1) = tmpl
          nocc(j+1) = tmpo
          exit
        else
          nnum(j+1) = nnum(j)
          lnum(j+1) = lnum(j)
          nocc(j+1) = nocc(j)
        end if
      end do
    end if
  end do
end subroutine get_occ

!
! Return a radial grid with weights
!
subroutine rlegendre(nrad, Z, rad, wgt)
    integer(ik)  :: nrad, Z
    real(rk)     :: rad(nrad), wgt(nrad)
    !
    character(2) :: atype
    real(rk)     :: r(nrad), w(nrad), rat

    call MathGetQuadrature('Legendre', nrad, r, w)
    atype = AtomElementSymbol(Z*1.0_rk)
    rat = 0.5_rk * AtomCovalentR(atype) / abohr
    ! map from (-1, 1) to (0, ...)
    rad = rat * (1 + r) / (1 - r)
    ! scale weights appropriately
    wgt = 8._rk * pi * w * rat * (rad / (1 - r))**2
    !wgt = pi * w * rat * (rad / (1 - r))**2
end subroutine rlegendre

!
!  Interpolate a point on a function f(x)
!
function interp(xpt, x, fx, nx, ityp, ordr)
    character(3) :: ityp
    integer(ik)  :: nx, ordr
    real(rk)     :: xpt, x(nx)
    real(ark)    :: fx(nx), interp
    !
    integer(ik)  :: i, ix, ind, n(ordr+1)

    ! find the nearest index
    find_ind: do ix = 2,nx
        if (xpt < x(ix) .or. ix == nx) then
            ind = ix - 1
            exit
        end if
    end do find_ind
    ! get the range of indices for interpolation
    if (ind < ordr/2 + 1) then
        n = (/(i, i = 1,ordr+1)/)
    else if (ind+ordr/2+mod(ordr,2) > nx) then
        n = (/(i, i = nx-ordr,nx)/)
    else
        n = (/(i, i = ind-ordr/2,ind+ordr/2+mod(ordr,2))/)
    end if
    select case (ityp)
        case ("pol")
            ! polynomial interpolation
            interp = 0.
            do i = 1,ordr+1
                ! don't bother adding zeros
                if (fx(n(i)) /= 0) then
                    interp = interp + fx(n(i))*x_li(xpt, x, nx, i, n, ordr)
                end if
            end do
        case ("exp")
            ! exponential interpolation
            interp = 1.
            do i = 1,ordr+1
                ! a f(x) value of 0 kills the interpolation
                if (fx(n(i)) == 0) then
                    interp = 0
                    return
                else
                    interp = interp * fx(n(i))**x_li(xpt, x, nx, i, n, ordr)
                end if
            end do
        case default
            write (out,"('interp: Error unrecognized type ',s)") ityp
            stop "interp - unrecognized type"
    end select
end function interp

!
!  Find the x weighting factor for interpolation
!
function x_li(xpt, x, nx, i, n, ordr)
    integer(ik) :: i, ordr, nx, n(ordr+1)
    real(rk)    :: xpt, x(nx), x_li
    !
    integer(ik) :: j
    x_li = 1.
    do j = 1,ordr+1
        if (j /= i) then
            x_li = x_li * (xpt - x(n(j))) / (x(n(i)) - x(n(j)))
        end if
    end do
end function x_li

!
! Return the cumulative sum of a 1D array
!
function cumsum(arr, nx)
  integer(ik) :: nx
  real(rk)    :: arr(nx), cumsum(nx)
  !
  integer(ik) :: i

  cumsum(1) = arr(1)
  do i = 2,nx
    cumsum(i) = arr(i) + cumsum(i-1)
  end do
end function cumsum

end module atomdens
