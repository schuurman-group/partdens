!
!  A module to write gamess data to a molden file (courtesy of Simon Neville)
!
module molden

  implicit none

contains

!######################################################################
! write_molden: writes norb orbitals (MOs, NTOs, etc.) held in the
!               array orb to the molden file named filename
!######################################################################

  subroutine write_molden(gam,filename,nao,norb,orb,ener,occ)

    use accuracy
    use gamess_internal

    implicit none

    integer                          :: norb,nao
    integer                          :: imolden
    integer                          :: i,j,k,pk,p1,p2,np,iang,iorb,&
                                        count
    real(rk), dimension(nao,norb)    :: orb,orb1
    real(rk), dimension(norb)        :: ener,occ
    real(rk)                         :: alpha,coeff
    real(rk), dimension(10,10)       :: ftransmat
    character(len=*)                 :: filename
    character(len=1), dimension(0:4) :: shlbl
    type(gam_structure)              :: gam

!-----------------------------------------------------------------------
! Set shell labels
!-----------------------------------------------------------------------
    shlbl(0:4)=(/ 's','p','d','f','g' /)

!-----------------------------------------------------------------------
! Transformation matrix for the reordering of the f-functions to
! comply with the molden ordering of:
!
!  1    2    3    4    5    6    7    8    9   10
! xxx, yyy, zzz, xyy, xxy, xxz, xzz, yzz, yyz, xyz
!-----------------------------------------------------------------------
    ftransmat=0.0d0
    ftransmat(1,1)=1.0d0
    ftransmat(2,2)=1.0d0
    ftransmat(3,3)=1.0d0
    ftransmat(6,4)=1.0d0
    ftransmat(4,5)=1.0d0
    ftransmat(5,6)=1.0d0
    ftransmat(8,7)=1.0d0
    ftransmat(9,8)=1.0d0
    ftransmat(7,9)=1.0d0
    ftransmat(10,10)=1.0d0

!-----------------------------------------------------------------------
! Make a copy of the orbitals and rearrange any f-functions to
! correspond to the molden ordering
!-----------------------------------------------------------------------
    orb1=orb

    count=0
    do i=1,gam%natoms

       do j=1,gam%atoms(i)%nshell

          ! Angular momentum quantum number and primitive indices
          ! for the current shell
          iang=gam%atoms(i)%sh_l(j)
          p1=count+1
          p2=count+gam_orbcnt(iang)

          ! If we are at a shell of f-functions, then rearrange
          ! to correspond to the molden ordering
          if (iang.eq.3) then
             orb1(p1:p2,:)=matmul(ftransmat,orb1(p1:p2,:))
          endif

          do k=ang_loc(iang),ang_loc(iang)+gam_orbcnt(iang)-1
             count=count+1
          enddo

       enddo

    enddo

!----------------------------------------------------------------------
! Open the molden file
!----------------------------------------------------------------------
    imolden = 11
    open(imolden,file=filename,form='formatted',status='unknown')

!----------------------------------------------------------------------
! Preamble
!----------------------------------------------------------------------
    write(imolden,'(a)') '[Molden Format]'
    write(imolden,'(a)') '[Title]'

!-----------------------------------------------------------------------
! Cartesian coordinates
!-----------------------------------------------------------------------
    write(imolden,'(a)') '[Atoms] Angs'
    do i=1,gam%natoms
       write(imolden,'(1x,a2,2x,i3,2x,i3,3(2x,F10.7))') &
            gam%atoms(i)%name,i,int(gam%atoms(i)%znuc),&
            (gam%atoms(i)%xyz(j),j=1,3)
    enddo

!-----------------------------------------------------------------------
! AO basis set section
!-----------------------------------------------------------------------
    write(imolden,'(a)') '[GTO]'

    ! Loop over atoms
    do i=1,gam%natoms

       ! Primitive counter for the current atom: pk
       pk=0

       ! Atom number
       if (i.gt.1) write(imolden,*)
       write(imolden,'(2x,i3,2x,a1)') i,'0'

       ! Loop over shells for the current atom
       do j=1,gam%atoms(i)%nshell

          ! Angular momentum quantum number and primitive indices
          ! for the current shell
          iang=gam%atoms(i)%sh_l(j)
          p1=gam%atoms(i)%sh_p(j)
          p2=gam%atoms(i)%sh_p(j+1)
          np=p2-p1

          ! Label and no. primitives for the current shell
          write(imolden,'(2x,a1,2x,i3,2x,a4)') shlbl(iang),&
               np,'1.00'

          ! Primitive exponents and coefficients for the current
          ! shell.
          ! Note that we output the original, unscaled,
          ! non-normalised primitive coefficients
          do k=1,np
             pk=pk+1
             alpha=gam%atoms(i)%p_zet(p1-1+k)
             coeff=gam%atoms(i)%p_c_orig(p1-1+k)
             write(imolden,*) alpha,coeff
          enddo

       enddo

    enddo

!----------------------------------------------------------------------
! Orbitals
!----------------------------------------------------------------------
    write(imolden,'(/,a)') '[MO]'

    do iorb=1,norb

       write(imolden,'(a)') 'Sym= 1'
       write(imolden,'(a,x,F10.6)') 'Ene=',ener(iorb)
       write(imolden,'(a)') 'Spin= Alpha'
       write(imolden,'(a,x,F10.6)') 'Occup= ',occ(iorb)

       count=0
       do i=1,gam%natoms

          do j=1,gam%atoms(i)%nshell

             ! Angular momentum quantum number and primitive indices
             ! for the current shell
             iang=gam%atoms(i)%sh_l(j)
             p1=count+1
             p2=count+gam_orbcnt(iang)

             do k=ang_loc(iang),ang_loc(iang)+gam_orbcnt(iang)-1
                count=count+1
                ! Note that here k tells us what type of function we
                ! are at (s, px, py,...) and count indexes the current
                ! AO basis function
                !
                ! Also, the coefficient that we output has to
                ! correspond to a normalised AO - hence the
                ! multiplication by ang_c(k)
                write(imolden,*) count,orb1(count,iorb)*ang_c(k)
             enddo

           enddo
        enddo

    enddo

!----------------------------------------------------------------------
! Close the molden file
!----------------------------------------------------------------------
    close(imolden)

    return

  end subroutine write_molden

!######################################################################

end module molden
