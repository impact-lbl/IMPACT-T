!----------------------------------------------------------------
! fftbench: benchmark and agreement check for the two IMPACT-T
! FFT backends (bundled FFTPACK 5.1 vs FFTW3).
!
! Exercises the two primitives, four1 (complex) and realft (real),
! that every FFT in the space-charge and wakefield solvers funnels
! through (see src/Appl/Field.f90 and src/Func/FFT.f90), using the
! same batched per-pencil call pattern as the solver: a batch of
! `size` 1D transforms per timing unit, i.e. one plane of a
! doubled space-charge grid.
!
! Timing is a forward+inverse round trip (the solver always does
! both), reported as microseconds per single 1D transform, after a
! warm-up call so that FFTW plan creation is excluded -- plans are
! cached and amortized over the whole run in production, too.
!
! Agreement is the maximum relative difference (L-inf, normalized
! by the L-inf norm of the FFTPACK result) between the two backends
! on identical inputs, tracked per function and direction and
! summarized at the end.
!
! Build:  make -C utilities/fftbench      (needs FFTW3; override
!         FC/FFTW_HOME as needed, see the Makefile)
! Run:    ./fftbench [size...]            sizes must be even; the
!         defaults 64 128 256 512 are the doubled grids for
!         32^3..256^3 space-charge meshes.
!
! The two backend modules are normally interchangeable (same name,
! the build system picks one); the Makefile here renames them via
! sed so both can live in this one comparison binary.
!----------------------------------------------------------------
program fftbench
  use fftpackbackendclass, only: four1_fftpack => four1, &
                                 realft_fftpack => realft
  use fftwbackendclass, only: four1_fftw => four1, &
                              realft_fftw => realft
  implicit none

  integer, parameter :: FOUR1_FWD = 1, FOUR1_INV = 2, &
                        REALFT_FWD = 3, REALFT_INV = 4
  character(len=14), parameter :: diffname(4) = &
      (/ 'four1  forward', 'four1  inverse', &
         'realft forward', 'realft inverse' /)
  ! Each size uses two cached plans per backend function; keep the
  ! total under the FFTW backend's MAXPLANS.
  integer, parameter :: MAXSIZES = 16
  integer :: sizes(MAXSIZES), nsizes
  real*8 :: maxdiff(4)
  integer*8 :: seed
  integer :: i, n, ierr
  character(len=32) :: arg

  nsizes = 0
  do i = 1, command_argument_count()
    call get_command_argument(i, arg)
    read(arg, *, iostat=ierr) n
    if (ierr .ne. 0 .or. n .lt. 4 .or. mod(n, 2) .ne. 0) then
      print*, 'fftbench: sizes must be even integers >= 4, got ', &
              trim(arg)
      stop
    endif
    if (nsizes .ge. MAXSIZES) then
      print*, 'fftbench: at most ', MAXSIZES, ' sizes'
      stop
    endif
    nsizes = nsizes + 1
    sizes(nsizes) = n
  enddo
  if (nsizes .eq. 0) then
    nsizes = 4
    sizes(1:4) = (/ 64, 128, 256, 512 /)
  endif

  seed = 20260728
  maxdiff = 0.0d0

  write(*,'(a)') 'fftbench: time per 1D transform '// &
      '(forward+inverse round trip, batched as in the solver)'
  write(*,'(a)') '  size batch func     fftpack(us)    fftw(us)'// &
      '  speedup  maxdiff(fwd)  maxdiff(inv)'
  do i = 1, nsizes
    call bench_four1(sizes(i))
    call bench_realft(sizes(i))
  enddo

  write(*,'(a)')
  write(*,'(a)') 'maximum relative difference per function '// &
      '(across all sizes):'
  do i = 1, 4
    write(*,'(2x,a,a,es10.2)') diffname(i), ' : ', maxdiff(i)
  enddo

contains

  subroutine bench_four1(n)
    integer, intent(in) :: n
    real*8, allocatable :: tmpl(:,:), a(:,:), b(:,:), spec(:,:)
    real*8 :: dfwd, dinv, tp, tw
    integer :: j, m

    m = n
    allocate(tmpl(2*n,m), a(2*n,m), b(2*n,m), spec(2*n,m))
    call fill(tmpl)

    a = tmpl
    b = tmpl
    do j = 1, m
      call four1_fftpack(a(:,j), n, -1)
      call four1_fftw(b(:,j), n, -1)
    enddo
    dfwd = reldiff(a, b)
    maxdiff(FOUR1_FWD) = max(maxdiff(FOUR1_FWD), dfwd)

    ! inverse agreement, on a common forward spectrum
    spec = a
    a = spec
    b = spec
    do j = 1, m
      call four1_fftpack(a(:,j), n, 1)
      call four1_fftw(b(:,j), n, 1)
    enddo
    dinv = reldiff(a, b)
    maxdiff(FOUR1_INV) = max(maxdiff(FOUR1_INV), dinv)

    tp = time_four1(1, n, m, tmpl)
    tw = time_four1(2, n, m, tmpl)
    write(*,'(i6,i6,2x,a,2f12.3,f9.2,2es14.2)') &
        n, m, 'four1 ', tp, tw, tp/tw, dfwd, dinv
  end subroutine bench_four1

  subroutine bench_realft(n)
    integer, intent(in) :: n
    real*8, allocatable :: tmpl(:,:), a(:,:), b(:,:), spec(:,:)
    real*8 :: dfwd, dinv, tp, tw
    integer :: j, m

    m = n
    allocate(tmpl(n,m), a(n,m), b(n,m), spec(n,m))
    call fill(tmpl)

    a = tmpl
    b = tmpl
    do j = 1, m
      call realft_fftpack(a(:,j), n, 1)
      call realft_fftw(b(:,j), n, 1)
    enddo
    dfwd = reldiff(a, b)
    maxdiff(REALFT_FWD) = max(maxdiff(REALFT_FWD), dfwd)

    spec = a
    a = spec
    b = spec
    do j = 1, m
      call realft_fftpack(a(:,j), n, -1)
      call realft_fftw(b(:,j), n, -1)
    enddo
    dinv = reldiff(a, b)
    maxdiff(REALFT_INV) = max(maxdiff(REALFT_INV), dinv)

    tp = time_realft(1, n, m, tmpl)
    tw = time_realft(2, n, m, tmpl)
    write(*,'(i6,i6,2x,a,2f12.3,f9.2,2es14.2)') &
        n, m, 'realft', tp, tw, tp/tw, dfwd, dinv
  end subroutine bench_realft

  ! Seconds -> microseconds per single 1D transform for backend
  ! `which` (1 = FFTPACK, 2 = FFTW). Repeats forward+inverse round
  ! trips (renormalizing so values stay bounded) and grows the
  ! repeat count until the measurement is long enough to trust.
  function time_four1(which, n, m, tmpl) result(us)
    integer, intent(in) :: which, n, m
    real*8, intent(in) :: tmpl(:,:)
    real*8 :: us, secs, s
    real*8, allocatable :: x(:,:)
    integer*8 :: iters, it, t0, t1, rate
    integer :: j

    allocate(x(2*n,m))
    x = tmpl
    call four1_which(which, x(:,1), n, -1)
    call four1_which(which, x(:,1), n, 1)
    s = 1.0d0/dble(n)
    iters = 1
    do
      x = tmpl
      call system_clock(t0)
      do it = 1, iters
        do j = 1, m
          call four1_which(which, x(:,j), n, -1)
        enddo
        do j = 1, m
          call four1_which(which, x(:,j), n, 1)
        enddo
        x = x*s
      enddo
      call system_clock(t1, rate)
      secs = dble(t1 - t0)/dble(rate)
      if (secs .gt. 0.2d0 .or. iters .ge. 16777216_8) exit
      iters = iters*2
    enddo
    us = secs/(2.0d0*dble(iters)*dble(m))*1.0d6
  end function time_four1

  function time_realft(which, n, m, tmpl) result(us)
    integer, intent(in) :: which, n, m
    real*8, intent(in) :: tmpl(:,:)
    real*8 :: us, secs, s
    real*8, allocatable :: x(:,:)
    integer*8 :: iters, it, t0, t1, rate
    integer :: j

    allocate(x(n,m))
    x = tmpl
    call realft_which(which, x(:,1), n, 1)
    call realft_which(which, x(:,1), n, -1)
    ! realft's inverse returns n/2 times the true inverse
    s = 2.0d0/dble(n)
    iters = 1
    do
      x = tmpl
      call system_clock(t0)
      do it = 1, iters
        do j = 1, m
          call realft_which(which, x(:,j), n, 1)
        enddo
        do j = 1, m
          call realft_which(which, x(:,j), n, -1)
        enddo
        x = x*s
      enddo
      call system_clock(t1, rate)
      secs = dble(t1 - t0)/dble(rate)
      if (secs .gt. 0.2d0 .or. iters .ge. 16777216_8) exit
      iters = iters*2
    enddo
    us = secs/(2.0d0*dble(iters)*dble(m))*1.0d6
  end function time_realft

  subroutine four1_which(which, data, n, isign)
    integer, intent(in) :: which, n, isign
    real*8 :: data(2*n)
    if (which .eq. 1) then
      call four1_fftpack(data, n, isign)
    else
      call four1_fftw(data, n, isign)
    endif
  end subroutine four1_which

  subroutine realft_which(which, data, n, isign)
    integer, intent(in) :: which, n, isign
    real*8 :: data(n)
    if (which .eq. 1) then
      call realft_fftpack(data, n, isign)
    else
      call realft_fftw(data, n, isign)
    endif
  end subroutine realft_which

  function reldiff(a, b) result(d)
    real*8, intent(in) :: a(:,:), b(:,:)
    real*8 :: d
    d = maxval(abs(a - b))/max(maxval(abs(a)), 1.0d-300)
  end function reldiff

  subroutine fill(x)
    real*8, intent(out) :: x(:,:)
    integer :: j, k
    do j = 1, size(x, 2)
      do k = 1, size(x, 1)
        x(k,j) = rand01()
      enddo
    enddo
  end subroutine fill

  ! Park-Miller LCG in [-1,1): deterministic, compiler-independent
  function rand01() result(r)
    real*8 :: r
    seed = mod(seed*48271_8, 2147483647_8)
    r = 2.0d0*dble(seed)/2147483647.0d0 - 1.0d0
  end function rand01

end program fftbench
