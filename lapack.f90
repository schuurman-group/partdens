module lapack
!
!  Simplistic type-agnostic LAPACK and LINPACK interface
!
  use accuracy
  implicit none

  interface lapack_gelss
    module procedure lapack_cgelss
    module procedure lapack_zgelss
    module procedure lapack_sgelss
    module procedure lapack_dgelss
  end interface ! lapack_gelss

  interface lapack_stev
    module procedure lapack_sstev
    module procedure lapack_dstev
  end interface ! lapack_stev

  interface lapack_sterf
    module procedure lapack_dsterf
    module procedure lapack_ssterf
  end interface ! lapack_sterf

  interface lapack_geev
    module procedure lapack_cgeev
    module procedure lapack_zgeev
  end interface ! lapack_geev

  interface lapack_heev
    module procedure lapack_cheev
    module procedure lapack_zheev
  end interface ! lapack_heev

  interface lapack_syev
    module procedure lapack_dsyev
    module procedure lapack_ssyev
  end interface ! lapack_syev

  interface lapack_ginverse
    module procedure lapack_ginverse_real
    module procedure lapack_ginverse_double
    module procedure lapack_ginverse_complex
    module procedure lapack_ginverse_doublecomplex
  end interface ! lapack_ginverse

  interface lapack_svd
    module procedure lapack_dgesvd
    module procedure lapack_zgesvd
  end interface ! lapack_svd  

  interface linpack_determinant
    module procedure linpack_determinant_double
  end interface ! linpack_determinant

  contains

  subroutine lapack_cgelss(a,b)
    complex, intent(inout) :: a(:,:)
    complex, intent(inout) :: b(:,:)

    external cgelss
    real                   :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex                :: work (50*max(size(a,dim=1),size(a,dim=2)))
    real                   :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call cgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_cgelss

  subroutine lapack_zgelss(a,b)
    complex(kind=drk), intent(inout) :: a(:,:)
    complex(kind=drk), intent(inout) :: b(:,:)

    external zgelss
    double precision       :: s    (   min(size(a,dim=1),size(a,dim=2)))
    complex(kind=drk)         :: work (50*max(size(a,dim=1),size(a,dim=2)))
    double precision       :: rwork( 5*min(size(a,dim=1),size(a,dim=2)))
    integer                :: rank, info
    integer                :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call zgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), rwork, info)

    if (info/=0) then
      write (out,"(' cgelss returned ',i8)") info
      stop 'lapack_cgelss - cgelss failed'
    end if
  end subroutine lapack_zgelss

  subroutine lapack_sgelss(a,b)
    real, intent(inout) :: a(:,:)
    real, intent(inout) :: b(:,:)

    external sgelss
    real                :: s    (   min(size(a,dim=1),size(a,dim=2)))
    real                :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call sgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0*spacing(1.0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' sgelss returned ',i8)") info
      stop 'lapack_sgelss - sgelss failed'
    end if
  end subroutine lapack_sgelss

  subroutine lapack_dgelss(a,b)
    double precision, intent(inout) :: a(:,:)
    double precision, intent(inout) :: b(:,:)

    external dgelss
    double precision    :: s    (   min(size(a,dim=1),size(a,dim=2)))
    double precision    :: work (50*max(size(a,dim=1),size(a,dim=2),size(b,dim=2)))
    integer             :: rank, info
    integer             :: na1, na2, nb1, nb2
    
    na1 = size(a,dim=1) ; na2 = size(a,dim=2)
    nb1 = size(b,dim=1) ; nb2 = size(b,dim=2)
    call dgelss(na1,na2,nb2,a(1:na1,1:na2),na1,b(1:nb1,1:nb2),nb1, &
                s, 100.0d0*spacing(1.0d0), rank, work, size(work), info)

    if (info/=0) then
      write (out,"(' dgelss returned ',i8)") info
      stop 'lapack_dgelss - dgelss failed'
    end if
  end subroutine lapack_dgelss

  subroutine lapack_sstev(d,e,z)
    real, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                  ! Out: Eigenvalues, ascending order
    real, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                  ! Out: Destroyed
    real, intent(out)   :: z(:,:) ! Out: Eigenvectors

    real    :: work(max(1,2*size(d)-2))
    integer :: info
    integer :: nz1, nz2
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call sstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' sstev returned ',i8)") info
      stop 'lapack_sstev - sstev failed'
    end if
  end subroutine lapack_sstev

  subroutine lapack_dstev(d,e,z)
    double precision, intent(inout) :: d(:)   ! In:  Diagonal elements of the matrix 
                                              ! Out: Eigenvalues, ascending order
    double precision, intent(inout) :: e(:)   ! In:  Sub-/super-diagonal elements of the matrix
                                              ! Out: Destroyed
    double precision, intent(out)   :: z(:,:) ! Out: Eigenvectors

    double precision :: work(max(1,2*size(d)-2))
    integer          :: info
    integer          :: nz1, nz2
    
    nz1 = size(z,dim=1) ; nz2 = size(z,dim=2)
    call dstev('V',size(d),d,e,z(1:nz1,1:nz2),nz1,work,info)

    if (info/=0) then
      write (out,"(' dstev returned ',i8)") info
      stop 'lapack_dstev - dstev failed'
    end if
  end subroutine lapack_dstev

  subroutine lapack_cgeev(h,e)
    complex, intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                      ! Out: Eigenvectors
    complex, intent(out)   :: e(:)    ! Out: Eigenvalues

    complex :: work(50*size(h,dim=2))
    real    :: rwork(3*size(h,dim=2))
    complex :: vl(1,1)
    complex :: vr(size(h,dim=2),size(h,dim=2))
    integer :: info
    integer :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cgeev('N','V',nh2,h(1:nh1,1:nh2),nh1,e(:),vl,1,vr,nh2,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cgeev returned ',i8)") info
      stop 'lapack_cgeev - cgeev failed'
    end if
    h(1:nh2,1:nh2) = vr
  end subroutine lapack_cgeev

  subroutine lapack_zgeev(h,e)
    complex(kind=drk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    complex(kind=drk), intent(out)   :: e(:)    ! Out: Eigenvalues

    complex(kind=drk)   :: work(50*size(h,dim=2))
    double precision :: rwork(3*size(h,dim=2))
    complex(kind=drk)   :: vl(1,1)
    complex(kind=drk)   :: vr(size(h,dim=2),size(h,dim=2))
    integer :: info
    integer :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zgeev('N','V',nh2,h(1:nh1,1:nh2),nh1,e(:),vl,1,vr,nh2,work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zgeev returned ',i8)") info
      stop 'lapack_zgeev - zgeev failed'
    end if
    h(1:nh2,1:nh2) = vr
  end subroutine lapack_zgeev

  subroutine lapack_cheev(h,e)
    complex, intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                      ! Out: Eigenvectors
    real, intent(out)   :: e(:)       ! Out: Eigenvalues

    complex :: work(50*size(h,dim=1))
    real    :: rwork(3*size(h,dim=1))
    integer :: info
    integer :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call cheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' cheev returned ',i8)") info
      stop 'lapack_cheev - cheev failed'
    end if
  end subroutine lapack_cheev

  subroutine lapack_zheev(h,e)
    complex(kind=drk), intent(inout) :: h(:,:)  ! In:  Hermitian matrix to be diagonalized
                                             ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)  ! Out: Eigenvalues

    complex(kind=drk)   :: work(50*size(h,dim=1))
    double precision :: rwork(3*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call zheev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),rwork,info)
    if (info/=0) then
      write (out,"(' zheev returned ',i8)") info
      stop 'lapack_zheev - zheev failed'
    end if
  end subroutine lapack_zheev

  subroutine lapack_dsyev(h,e)
    double precision, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                               ! Out: Eigenvectors
    double precision, intent(out)   :: e(:)    ! Out: Eigenvalues

    double precision :: work(50*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call dsyev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),info)
    if (info/=0) then
      write (out,"(' dsyev returned ',i8)") info
      stop 'lapack_dsyev - dsyev failed'
    end if
  end subroutine lapack_dsyev

  subroutine lapack_ssyev(h,e)
    real, intent(inout) :: h(:,:)  ! In:  symmetric matrix to be diagonalized
                                   ! Out: Eigenvectors
    real, intent(out)   :: e(:)    ! Out: Eigenvalues

    real             :: work(50*size(h,dim=1))
    integer          :: info
    integer          :: nh1, nh2
    
    nh1 = size(h,dim=1) ; nh2 = size(h,dim=2)
    call ssyev('V','U',nh1,h(1:nh1,1:nh2),nh1,e(:),work,size(work),info)
    if (info/=0) then
      write (out,"(' ssyev returned ',i8)") info
      stop 'lapack_ssyev - ssyev failed'
    end if
  end subroutine lapack_ssyev

  subroutine lapack_dgesvd(a,s,u,vth)
    double precision, intent(inout) :: a  (:,:)  ! In:  Matrix to be decomposed
                                                 ! Out: Content destroyed
    double precision, intent(out)   :: s  (:)    ! Out: Singular values
    double precision, intent(out)   :: u  (:,:)  ! Out: Left singular vectors
    double precision, intent(out)   :: vth(:,:)  ! Out: Right singular vectors, transposed & conjugated
                                                 ! The overall result is A = U S VTH

    character(len=1)              :: jobu        ! Either 'A' (all) or 'S' (singular), depending on the
    character(len=1)              :: jobvth      ! sizes of u and vth arrays
    double precision              :: lwq(1)
    double precision, allocatable :: work(:)
    integer                       :: info, lwork
    integer                       :: m, n, lda, ldu, ldvth, nsing
   
    m     = size(a,dim=1) 
    n     = size(a,dim=2) 
    lda   = m
    nsing = min(m,n)
    ldu   = size(u,dim=1)
    ldvth = size(vth,dim=1)
    !
    if (size(s)<nsing) stop 'lapack%lapack_dgesvd - array s is too small'
    !
    if (size(u,dim=1)<m    ) stop 'lapack%lapack_dgesvd - array u is too small (1)'
    if (size(u,dim=2)<nsing) stop 'lapack%lapack_dgesvd - array u is too small (2)'
    jobu = 'S'
    if (size(u,dim=2)>=m) jobu = 'A'
    !
    if (size(vth,dim=2)<n    ) stop 'lapack%lapack_dgesvd - array vth is too small (2)'
    if (size(vth,dim=1)<nsing) stop 'lapack%lapack_dgesvd - array vth is too small (1)'
    jobvth = 'S'
    if (size(vth,dim=1)>=n) jobvth = 'A'
    !
    call dgesvd(jobu,jobvth,m,n,a,lda,s,u,ldu,vth,ldvth,lwq, -1,   info)
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8,' for workspace query')") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    lwork = 1+nint(lwq(1))
    allocate (work(lwork),stat=info)
    if (info/=0) then
      write (out,"(' Error ',i8,' allocating ',i10,'-element array for dgesvd')") info, lwork
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    !
    call dgesvd(jobu,jobvth,m,n,a,lda,s,u,ldu,vth,ldvth,work,lwork,info)
    if (info/=0) then
      write (out,"(' dgesvd returned ',i8)") info
      stop 'lapack_dgesvd - dgesvd failed'
    end if
    deallocate (work)
  end subroutine lapack_dgesvd

  subroutine lapack_zgesvd(a,s,u,vth)
    complex(kind=drk), intent(inout) :: a  (:,:)  ! In:  Matrix to be decomposed
                                                  ! Out: Content destroyed
    double precision, intent(out)    :: s  (:)    ! Out: Singular values
    complex(kind=drk), intent(out)   :: u  (:,:)  ! Out: Left singular vectors
    complex(kind=drk), intent(out)   :: vth(:,:)  ! Out: Right singular vectors, transposed & conjugated
                                                  ! The overall result is A = U S VTH

    complex(kind=drk)   :: work(50*max(size(a,dim=1),size(a,dim=2)))
    double precision :: rwork(5*size(a,dim=1))
    integer          :: info, lwork
    integer          :: m, n, lda, ldu, ldvth
   
    m     = size(a,dim=1) ; lda = m
    n     = size(a,dim=2) 
    ldu   = size(u,dim=1)
    ldvth = size(vth,dim=1)
    lwork = size(work)
    if (size(s)<min(m,n))               stop 'lapack%lapack_zgesvd - array s is too small'
    if (ldu<m .or. size(u,dim=2)<m)     stop 'lapack%lapack_zgesvd - array u is too small'
    if (ldvth<n .or. size(vth,dim=2)<n) stop 'lapack%lapack_zgesvd - array vth is too small'
    !
    call zgesvd('A','A',m,n,a,lda,s,u,ldu,vth,ldvth,work,lwork,rwork,info)
    if (info/=0) then
      write (out,"(' zgesvd returned ',i8)") info
      stop 'lapack_zgesvd - zgesvd failed'
    end if
  end subroutine lapack_zgesvd

  subroutine lapack_dsterf(a,b)
    double precision, intent(inout) :: a(:) ! In: Diagonal elements of the tri-diagonal matrix
                                            ! Out: Eigenvalues, in the ascending order
    double precision, intent(inout) :: b(:) ! In: Sub-diagonal elements of the tri-diagonal matirx
                                            ! Out: Destroyed
    !
    integer :: na, nb
    integer :: info
    !
    na = size(a)
    nb = size(b)
    if (na/=nb+1) then
      write (out,"('lapack_dsterf: inconsistent array sizes: diagonal ',i6,' subdiagonal ',i6)") na, nb
      stop 'lapack_dsterf - bad input'
    end if
    call dsterf(na,a,b,info)
    if (info/=0) then
      write (out,"(' dsterf returned ',i8)") info
      stop 'lapack_dsterf - dsterf failed'
    end if
  end subroutine lapack_dsterf

  subroutine lapack_ssterf(a,b)
    real, intent(inout) :: a(:) ! In: Diagonal elements of the tri-diagonal matrix
                                ! Out: Eigenvalues, in the ascending order
    real, intent(inout) :: b(:) ! In: Sub-diagonal elements of the tri-diagonal matirx
                                ! Out: Destroyed
    !
    integer :: na, nb
    integer :: info
    !
    na = size(a)
    nb = size(b)
    if (na/=nb+1) then
      write (out,"('lapack_ssterf: inconsistent array sizes: diagonal ',i6,' subdiagonal ',i6)") na, nb
      stop 'lapack_ssterf - bad input'
    end if
    call ssterf(na,a,b,info)
    if (info/=0) then
      write (out,"(' ssterf returned ',i8)") info
      stop 'lapack_ssterf - ssterf failed'
    end if
  end subroutine lapack_ssterf

  subroutine lapack_ginverse_real(amat,power_)
    real, intent(inout)        :: amat(:,:) ! In: matrix to invert
                                            ! Out: generalized inverse of the matrix
    real, intent(in), optional :: power_    ! In: power of the inverse; 1/2 if omitted
                                            !     (corresponding to A^(-1))
    !
    real :: eval(size(amat,dim=1))
    real :: evec(size(amat,dim=1),size(amat,dim=1))
    real :: eps
    real :: power
    !
    power = 0.5
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))
    !
  end subroutine lapack_ginverse_real

  subroutine lapack_ginverse_complex(amat,power_)
    complex, intent(inout)     :: amat(:,:) ! In: matrix to invert
                                            ! Out: generalized inverse of the matrix
    real, intent(in), optional :: power_    ! In: power of the inverse; 1/2 if omitted
                                            !     (corresponding to A^(-1))
    !
    real    :: eval(size(amat,dim=1))
    complex :: evec(size(amat,dim=1),size(amat,dim=1))
    real    :: eps
    real    :: power
    !
    power = 0.5
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_heev(evec(:,:),eval(:))
    eps = 100.0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(conjg(evec)))
    !
  end subroutine lapack_ginverse_complex

  subroutine lapack_ginverse_double(amat,power_)
    double precision, intent(inout) :: amat(:,:)     ! In: matrix to invert
                                                     ! Out: generalized inverse of the matrix
    double precision, intent(in), optional :: power_ ! In: power of the inverse; 1/2 if omitted
                                                     !     (corresponding to A^(-1))

    double precision :: eval(size(amat,dim=1))
    double precision :: evec(size(amat,dim=1),size(amat,dim=1))
    double precision :: eps
    double precision :: power
    !
    power = 0.5d0
    if (present(power_)) power = power_
    !
    evec = amat
    call lapack_syev(evec(:,:),eval(:))
    eps = 100.0d0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0d0 / eval**power
    elsewhere  
      eval = 0.0d0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(evec))
    !
  end subroutine lapack_ginverse_double

  subroutine lapack_ginverse_doublecomplex(amat,power_)
    complex(kind=drk), intent(inout)     :: amat(:,:)   ! In: matrix to invert
                                                     ! Out: generalized inverse of the matrix
    double precision, intent(in), optional :: power_ ! In: power of the inverse; 1/2 if omitted
                                                     !     (corresponding to A^(-1))
    !
    double precision :: eval(size(amat,dim=1))
    complex(kind=drk)   :: evec(size(amat,dim=1),size(amat,dim=1))
    double precision :: eps
    double precision :: power
    !
    power = 0.5d0
    if (present(power_)) power = power_
    !
    evec  = amat
    call lapack_heev(evec(:,:),eval(:))
    eps = 100.0d0*spacing(maxval(eval))
    where (abs(eval)>eps) 
      eval = 1.0d0 / eval**power
    elsewhere  
      eval = 0.0
    end where
    evec = evec * spread(eval,dim=1,ncopies=size(evec,dim=1))
    amat = matmul(evec,transpose(conjg(evec)))
    !
  end subroutine lapack_ginverse_doublecomplex

  function linpack_determinant_double(mat) result(det)
    double precision, intent(in)  :: mat(:,:) ! Matrix to compute determinant for
    double precision              :: det
    !
    double precision, allocatable :: tm(:,:)
    integer                       :: order, info
    integer                       :: ipvt(size(mat,dim=1))
    double precision              :: work(size(mat,dim=1))
    double precision              :: detx(2)
    external                      :: dgefa, dgedi
    !
    order = size(mat,dim=1)
    if (size(mat,dim=2)/=order) then
      write (out,"('Determinant requested for a non-square matrix: ',2i5,'. Bummer.')") &
             order, size(mat,dim=2)
      stop 'lapack%linpack_determinant_double - bad input'
    end if
    !
    allocate (tm(order,order),stat=info)
    if (info/=0) then
      write (out,"('Error ',i5,' allocating order-',i5,' matrix.')") info, order
      stop 'lapack%linpack_determinant_double - no memory'
    end if
    tm = mat
    !
    call dgefa(tm,order,order,ipvt,info)
    !
    call dgedi(tm,order,order,ipvt,detx,work,10)
    !
    ! tm = mat
    ! call lapack_dsyev(tm,work)
    ! write (out,"('Diagonalization gives: ',g40.20/' linpack gives ',g40.20/' Diff = ',g40.20)") &
    !        product(work), detx(1) * 10.0d0**detx(2), product(work) - detx(1) * 10.0d0**detx(2)
    !
    deallocate(tm,stat=info)
    if (info/=0) then
      write (out,"('Error ',i5,' deallocating order-',i5,' matrix.')") info, order
      stop 'lapack%linpack_determinant_double - memory deallocation failed'
    end if
    !
    det = detx(1) * 10.0d0**detx(2)
  end function linpack_determinant_double

end module lapack
