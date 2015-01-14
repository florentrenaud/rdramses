program rdramses

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Florent Renaud    !!!!!    since 5 Aug 2011     !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Based on Romain Teyssier's amr2map and part2map !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Documentation, help:                                  !!
!!   http://irfu.cea.fr/Pisp/florent.renaud/rdramses.php !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
#ifdef FFTW
  include "fftw3.f"
#endif

#ifdef MPIMODE
  include 'mpif.h'

  integer mpi_rank, mpi_size, mpi_ierror, mpi_status(MPI_STATUS_SIZE)
  real(kind=8) mpi_dx, mpi_dy, mpi_tmpxmin, mpi_tmpymin, mpi_tmpxmax, mpi_tmpymax
  character(len=128) mpi_outname
#endif

  integer mpi_ncol, mpi_nrow ! this must be outside of the ifdef MPIMODE because of the reading of the parameters



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,parameter::NPARTMAX=60000000 ! increase this value if necessary !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! and don't forget to re-compile!  !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  integer::i, j, k, n, impi, icpu, ilevel, iidim, ivar, ind
  character(len=5)::nchar,ncharcpu,nlev
  character(len=128)::filename, repository, out='', outval, outvalunit
  character(len=14)::void_str
  real(KIND=8)::xmin=0, xmax=1, ymin=0, ymax=1, zmin=0, zmax=1, cube=0, cubex=0, cubey=0, cubez=0, xcen=0, ycen=0, zcen=0, rhomin=1e-20, rhomax=1e20
  real(kind=8)::xxmin, xxmax, yymin, yymax, zzmin, zzmax
  real(kind=8)::offsetx=0, offsety=0, offsetz=0
  character(len=1)::dir
  integer::lmax=0, typ=1
  logical::kpc=.false., hydro=.false., makemap=.true., maxrho=.false., maxval=.false., sum=.false., psd=.false., psdm=.false., grav=.false., gravinp=.false., q=.false., p=.false., part=.false., id=.false., ns=.false., cen=.false.

  integer::ncpu, ndim, nx, ny, nz, nlevelmax, ngridmax, nboundary, ngrid_current
  integer::twotondim, levelmin, bit_length, maxdom, ndom, ncpu_read, nvarh
  integer::ngrida
  integer::idim, jdim, kdim, imin, imax, jmin, jmax, kmin, kmax
  integer::ix, iy, iz
  integer::nx_full, ny_full, nz_full
  integer::nx_sample=0, ny_sample=0
  real(kind=8)::dkey, order_min, dx, dmax, dxline, boxlen, t, ljeans, mass, pi
  real::weight
  character(len=80)::ordering
  logical::ok
  
  integer,dimension(1:8)::idom, jdom, kdom, cpu_min, cpu_max
  real(kind=8),dimension(1:3)::xbound
  real(kind=8),dimension(1:8,1:3)::xc
  integer,dimension(:,:),allocatable::ngridfile, ngridlevel, ngridbound
  integer,dimension(:),allocatable::cpu_list
  real(kind=8),dimension(:),allocatable::bound_key
  real(kind=8),dimension(1:8)::bounding_min, bounding_max
  logical,dimension(:),allocatable::cpu_read

  real(kind=8),dimension(:,:),allocatable::x, xg
  real(kind=8),dimension(:,:,:),allocatable::var, varg
  real(kind=8),dimension(:),allocatable::rho, map
  logical,dimension(:),allocatable::ref
  integer,dimension(:,:),allocatable::son

  type level
    integer::ilevel
    integer::ngrid
    real(KIND=8),dimension(:,:),pointer::map
    real(KIND=8),dimension(:,:),pointer::rho
    integer::imin
    integer::imax
    integer::jmin
    integer::jmax
    integer::kmin
    integer::kmax
  end type level

  type(level),dimension(1:100)::grid

  real(KIND=4),dimension(:,:),allocatable::tmpmap
  integer::nxmap, nymap
  real,dimension(:,:),allocatable::map2


  ! random sub-cells particles
  integer,parameter::seed=9876543
  real(KIND=8)::rx,ry,rz
  real(kind=8)::rand=0.0
  integer::npartcell

  ! pdf
  integer::ipdf, pdf=0
  real(kind=8)::lrhoampli, lrhomin, lrhomax
  real(kind=8),dimension(:),allocatable::pdfhist

  ! psd
  integer*8 plan_forward, plan_backward
  integer::ipsd, jpsd, binpsd, psdmap
  character(len=2)::charlmax
  real(kind=8)::fpsd, psdnorm, freq
  integer,dimension(:),allocatable::npsd
  real(kind=8),dimension(:),allocatable::ppsd
  complex(kind=8),dimension(:,:),allocatable::psdout, psdin

  ! gravity, tides
  real(kind=8)::fgrav
  real(kind=8),dimension(:,:,:,:),allocatable::tt, sh
  real(kind=8)::lambda1,lambda2,lambda3

  ! velocity dispersion
  real(kind=8),dimension(:),allocatable::vdisp2, vdisp1, vdispm
  real(kind=8)::vdisptmp

  ! turbulence
  real(kind=8),dimension(:,:),allocatable::nablav,nablav1,nablav2,nablav3


  ! particles
  real(kind=8),dimension(1:NPARTMAX,1:3)::xselect
  real(kind=8),dimension(1:NPARTMAX)::mselect
  integer::npart, nselect
  integer::pstep=1, pix=512
  real(kind=8)::pamin=0.0D0, pamax=1.0D20, pmmin=0.0D0, pmmax=1.0D20
  real(kind=8),dimension(:,:),allocatable::v ! x already defined
  real(kind=8),dimension(:),allocatable::age, m
  integer,dimension(:),allocatable::part_id
  integer::bx, by
  real(kind=8)::dy, agetmp, atmp ! dx already defined

  ! sfr
  integer::sfrn, isfr
  real(kind=8)::eta_sn=0.0, f_w, sfr=0.0, agesn = 10.0 ! agesn in Myr
  real(kind=8),dimension(:),allocatable::sfrhist

#ifdef FITS
  integer::status, unit, blocksize, bitpix, naxis
  integer,dimension(2)::naxes
  integer::group,fpixel,nelements
  logical::simple,extend
#endif
  
  integer::iargc
  character(len=8)::opt
  character(len=128)::arg
  namelist /rdr/ out, lmax, xmin, xmax, ymin, ymax, zmin, zmax, xcen, ycen, zcen, kpc, cube, cubex, cubey, cubez, cen, hydro, q, rhomin, rhomax, rand, grav, pdf, dir, maxval, maxrho, sum, typ, psd, psdm
! particle nml

  real(kind=8)::scale_l,scale_t,scale_d, scale_lkpc, scale_lpc, scale_dhcc, scale_tmyr, scale_tyr, scale_msun, scale_vkms, scale_temk, scale_forc, scale_surf, scale_map


  call srand(seed)
  pi = acos(-1.0D0)

! order of cells : left-bottom-rear, right-bottom-rear, left-top-rear, rigth-top-rear, left-bottom-front, right-bottom-front, left-top-front, right-top-front
!=======================================================================
!!! Read parameters and options

    n = iargc()
    i = 1
    do while(i.le.n)
      call getarg(i,opt)
      select case (opt)
        case ('-inp')
          call getarg(i+1,arg)
          repository = trim(arg)
        case ('-nml')
          call getarg(i+1,arg)
          open(1,file=trim(arg))
          read(1,rdr)
          close(1)
        case ('-out')
          call getarg(i+1,arg)
          out = trim(arg)
        case ('-lmax')
          call getarg(i+1,arg)
          read (arg,*) lmax
        case ('-xmin')
          call getarg(i+1,arg)
          read (arg,*) xmin
        case ('-xmax')
          call getarg(i+1,arg)
          read (arg,*) xmax
        case ('-ymin')
          call getarg(i+1,arg)
          read (arg,*) ymin
        case ('-ymax')
          call getarg(i+1,arg)
          read (arg,*) ymax
        case ('-zmin')
          call getarg(i+1,arg)
          read (arg,*) zmin
        case ('-zmax')
          call getarg(i+1,arg)
          read (arg,*) zmax
        case ('-xcen')
          call getarg(i+1,arg)
          read (arg,*) xcen
        case ('-ycen')
          call getarg(i+1,arg)
          read (arg,*) ycen
        case ('-zcen')
          call getarg(i+1,arg)
          read (arg,*) zcen
        case ('-kpc')
          kpc = .true.
          i = i-1
        case ('-cube')
          call getarg(i+1,arg)
          read (arg,*) cube
        case ('-cen')
          cen = .true.
          i = i-1
        case ('-d')
          call getarg(i+1,arg)
          read (arg,*) cubex
          call getarg(i+2,arg)
          read (arg,*) cubey
          call getarg(i+3,arg)
          read (arg,*) cubez
          i = i+2
        case ('-hydro')
          hydro = .true.
          i = i-1
        case ('-q')
          q = .true.
          i = i-1
        case ('-rhomin')
          call getarg(i+1,arg)
          read (arg,*) rhomin
        case ('-rhomax')
          call getarg(i+1,arg)
          read (arg,*) rhomax
        case ('-rand')
          call getarg(i+1,arg)
          read (arg,*) rand
        case ('-grav')
          grav=.true.
          i = i-1
        case ('-pdf')
          call getarg(i+1,arg)
          read (arg,*) pdf
        case ('-dir')
          call getarg(i+1,arg)
          dir = trim(arg) 
        case ('-maxval')
          maxval = .true.
          i = i-1
        case ('-maxrho')
          maxrho = .true.
          i = i-1
        case ('-sum')
          sum = .true.
          i = i-1
        case ('-typ')
          call getarg(i+1,arg)
          read (arg,*) typ
        case ('-psd')
          psd = .true.
          i = i-1
        case ('-psdm')
          psd = .true.
          psdm = .true.
          i = i-1
        case ('-p')
          p = .true.
          i = i-1
        case ('-pix')
          call getarg(i+1,arg)
          read (arg,*) pix
        case ('-part')
          part = .true.
          i = i-1
        case ('-ns')
          ns = .true.
          i = i-1
        case ('-id')
          id = .true.
          i = i-1
        case ('-sfr')
          call getarg(i+1,arg)
          read (arg,*) sfr
        case ('-eta_sn')
          call getarg(i+1,arg)
          read (arg,*) eta_sn
        case ('-f_w')
          call getarg(i+1,arg)
          read (arg,*) f_w
        case ('-pstep')
          call getarg(i+1,arg)
          read (arg,*) pstep
       case ('-pmmin')
          call getarg(i+1,arg)
          read (arg,*) pmmin
       case ('-pmmax')
          call getarg(i+1,arg)
          read (arg,*) pmmax
       case ('-pamin')
          call getarg(i+1,arg)
          read (arg,*) pamin
       case ('-pamax')
          call getarg(i+1,arg)
          read (arg,*) pamax


        case ('-ncol')
          call getarg(i+1,arg)
          read (arg,*) mpi_ncol
        case ('-nrow')
          call getarg(i+1,arg)
          read (arg,*) mpi_nrow


        case default
          print '("unknown option ",a8," ignored")', opt
          i = i-1
      end select
      i = i+2
    end do


    ! set variables (must be done here in case of nml reading)
    if(grav) gravinp = .true.
    if(psdm) psd = .true.
    if(ns .OR. id) part = .true.
    if(part) p = .true.
    if(sfr > 0) then
      ns = .true.
      p = .true.
    endif
    if(hydro .OR. q .OR. (rand > 0) .OR. grav .OR. (pdf > 0) .OR. ns .OR. id .OR. (sfr > 0)) makemap = .false.
    if(cube > 0 .or. cubex > 0) kpc = .true.
    if(len(TRIM(out)).ne.0) out = '_'//TRIM(out)

!=======================================================================  

!!! read info, compute units
  i=INDEX(repository,'output_')
  nchar=repository(i+7:i+13)

  filename=TRIM(repository)//'/info_'//TRIM(nchar)//'.txt'
  inquire(file=filename, exist=ok)
  if(.not. ok)then
    write(*,*) "Error: ", trim(filename), " not found"
    stop
  endif

  open(unit=1, file=filename, form='formatted', status='old')
  read(1,'(A14,I11)') void_str, ncpu
  read(1,'(A14,I11)') void_str, ndim
  read(1,'(A14,I11)'), void_str, levelmin
  read(1,'(A14,I11)'), void_str, nlevelmax
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,'(A14,E23.15)'), void_str, boxlen
  read(1,'(A14,E23.15)'), void_str, t
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,'(A14,E23.15)'), void_str, scale_l
  read(1,'(A14,E23.15)'), void_str, scale_d
  read(1,'(A14,E23.15)'), void_str, scale_t
  read(1,*)
  read(1,'(A14,A80)'), void_str, ordering
  read(1,*)

  if(ndim .NE. 3) then
    close(1)
    write(*,*) 'Unfortunately, rdramses is not designed to read non-3D output. Sorry!'
    stop
  endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Unit conversion

! At this point, scale_l, scale_d and scale_t contains the conversion factors from RAMSES output to CGS
! ramses * scale_l => cm
! ramses * scale_d => g/cm^3
! ramses * scale_t => s
!!! rdramses computes other conversion factors toward useful units, like kpc and Myr
!!! E.g. to convert the snapshot timestamp from RAMSES units to Myr: t * scale_tmyr
!!! E.g. to convert the user-defined minimum age in Myr into RAMES units: pamin / scale_tmyr


!from units.f90
!  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3
!  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)
!  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)


!! length in kpc
scale_lkpc = scale_l / 3.085677581282D21 ! (length [cm]) / (3e21 [cm/kpc])

!! length in pc
scale_lpc = scale_l / 3.085677581282D18 ! (length [cm]) / (3e18 [cm/pc])

!! density in hydrogen atoms per cm^3
scale_dhcc = scale_d / 1.66D-24 * 0.76 ! (density [g/cm^3]) / (mass_of_hydrogen [g]) * (hydrogen_aboundance [%])

!! time in Myr
scale_tmyr = scale_t / 3.15576D13 ! (time [s]) / (3e13 [s/Myr])

!! time in yr
scale_tyr = scale_t / 3.15576D7 ! (time [s]) / (3e7 [s/yr])

!! mass in Msun
scale_msun = scale_d * scale_l**3 / 1.9891D33 ! (density [g/cm^3]) * (length [cm])^3 / (solar_mass [g/Msun])

!! velocity in km/s
scale_vkms = scale_l / scale_t / 1D5 ! (length [cm]) / (time [s]) / (1e5 [cm/km])

!! temperature in Kelvin
scale_temk = (scale_l / scale_t)**2 * 1.66D-24 / 1.3806200D-16 ! (velocity [cm/s])^2 * (mass_of_hydrogen [g]) / boltzmann

!! acceleration (=specific force) in kpc/Myr^2
scale_forc = scale_l / (scale_t)**2 / 3.085677581282D21 * (3.15576D13)**2 ! (length [cm]) / (time [s])^2 / (3e21 [cm/kpc]) * (3e13 [s/Myr])^2

!! surface_density in Msun/pc^2
scale_surf = scale_d * scale_l  / 1.9891D33 * (3.085677581282D18)**2 ! (density [g/cm^3]) * (length [cm]) / (solar_mass [g/Msun]) * (3e18 [cm/pc])^2


!! conversion of user-defined values into RAMSES units
  rand = rand / scale_msun
  rhomin = rhomin / scale_dhcc
  rhomax = rhomax / scale_dhcc
  pmmin = pmmin / scale_msun
  pmmax = pmmax / scale_msun
  pamin = pamin / scale_tmyr
  pamax = pamax / scale_tmyr
  agesn = agesn / scale_tmyr
  sfr = sfr / scale_tmyr
  xcen = xcen / scale_lkpc
  ycen = ycen / scale_lkpc
  zcen = zcen / scale_lkpc
  cubex = cubex / scale_lkpc
  cubey = cubey / scale_lkpc
  cubez = cubez / scale_lkpc
  cube = cube / scale_lkpc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if(lmax == 0) lmax = nlevelmax

  if(.not. p) then ! hydro / grav only
    allocate(cpu_list(1:ncpu))
    if(TRIM(ordering).eq.'hilbert')then
      allocate(bound_key(0:ncpu))
      allocate(cpu_read(1:ncpu))
      cpu_read=.false.
      do impi=1,ncpu
        read(1,'(I8,1X,E23.15,1X,E23.15)') i, bound_key(impi-1), bound_key(impi)
      end do
    endif
  endif ! end hydro/grav only

  close(1) ! close info file


!!! define extraction volume
  if(kpc) then
    if(cubex == 0.0) cubex = cube
    if(cubey == 0.0) cubey = cube
    if(cubez == 0.0) cubez = cube

    if(cubex > 0.0) then
      cubex = cubex/2.
      cubey = cubey/2.
      cubez = cubez/2.
      xmax = xcen + cubex
      xmin = xcen - cubex
      ymax = ycen + cubey
      ymin = ycen - cubey
      zmax = zcen + cubez
      zmin = zcen - cubez
    endif

    xmin = xmin/boxlen + 0.5
    xmax = xmax/boxlen + 0.5
    ymin = ymin/boxlen + 0.5
    ymax = ymax/boxlen + 0.5
    zmin = zmin/boxlen + 0.5
    zmax = zmax/boxlen + 0.5
  endif ! xmin and co are now in fraction of boxlen (between 0 and 1)

  if(cen) then
    offsetx = (xmin+xmax)/2.
    offsety = (ymin+ymax)/2.
    offsetz = (zmin+zmax)/2.
  else
    offsetx = 0.5
    offsety = 0.5
    offsetz = 0.5
  endif



!!! extraction volume splitting
#ifdef MPIMODE
  call MPI_INIT(mpi_ierror)

  call MPI_COMM_SIZE(MPI_COMM_WORLD, mpi_size, mpi_ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, mpi_rank, mpi_ierror)

  if(mpi_rank == 0)then
    if(mpi_size /= mpi_ncol*mpi_nrow) then
      write(*,*) 'Error: ncol * nrow must = Ncpu'
    endif
    
    mpi_dx = (xmax-xmin)/mpi_ncol ! size of sub-region
    mpi_dy = (ymax-ymin)/mpi_nrow
    
    do i=1, mpi_size-1  
      mpi_tmpxmin = xmin+mod(i,mpi_ncol)*mpi_dx
      mpi_tmpymin = ymin+floor(i/float(mpi_ncol))*mpi_dy
      mpi_tmpxmax = mpi_tmpxmin+mpi_dx
      mpi_tmpymax = mpi_tmpymin+mpi_dy

      call MPI_SEND(mpi_tmpxmin,1,MPI_DOUBLE_PRECISION,i,101,MPI_COMM_WORLD,mpi_ierror)
      call MPI_SEND(mpi_tmpymin,1,MPI_DOUBLE_PRECISION,i,101,MPI_COMM_WORLD,mpi_ierror)
      call MPI_SEND(mpi_tmpxmax,1,MPI_DOUBLE_PRECISION,i,101,MPI_COMM_WORLD,mpi_ierror)
      call MPI_SEND(mpi_tmpymax,1,MPI_DOUBLE_PRECISION,i,101,MPI_COMM_WORLD,mpi_ierror)
    enddo

    xmax=xmin+mpi_dx
    ymax=ymin+mpi_dy
  else
    call MPI_RECV(xmin,1,MPI_DOUBLE_PRECISION,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierror)
    call MPI_RECV(ymin,1,MPI_DOUBLE_PRECISION,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierror)
    call MPI_RECV(xmax,1,MPI_DOUBLE_PRECISION,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierror)
    call MPI_RECV(ymax,1,MPI_DOUBLE_PRECISION,0,101,MPI_COMM_WORLD,MPI_STATUS_IGNORE,mpi_ierror)
  endif

  write(mpi_outname, '(I5.5)')  mpi_rank
  out = out//'_cpu'//mpi_outname
#endif


!! projections and suffix
  if(makemap) then 
    if(.not. p) then
      if(maxval) out = '_mv'//out
      if(maxrho) out = '_mr'//out
      if(sum) out = '_sum'//out
    endif
  endif

  ! projections. (xmin = intrinsic, xxmin = projected)
  ! needed even if not makemap
  select case (dir)
    case ('x')
      out = '_x'//out
      idim=3
      jdim=2
      kdim=1
      xxmin=zmin
      xxmax=zmax
      yymin=ymin
      yymax=ymax
      zzmin=xmin
      zzmax=xmax
    case ('y')
      out = '_y'//out
      idim=1
      jdim=3
      kdim=2
      xxmin=xmin
      xxmax=xmax
      yymin=zmin
      yymax=zmax
      zzmin=ymin
      zzmax=ymax
    case default
      idim=1
      jdim=2
      kdim=3
      xxmin=xmin
      xxmax=xmax
      yymin=ymin
      yymax=ymax
      zzmin=zmin
      zzmax=zmax
  end select


  if(makemap .OR. q) then 
    if(p) then
      outval = 'part'
      outvalunit = 'Msun/pc3'
    else
      select case (typ)
        case (-1)
          outval = 'cpu'
          outvalunit = ' '
          scale_map = 1.0
        case (0)
          outval = 'level'
          outvalunit = ''
          scale_map = 1.0
        case (2) ! x-velocity
          outval = 'vx'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (3) ! y-velocity
          outval = 'vy'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (4) ! z-velocity
          outval = 'vz'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (5) ! Pressure
          outval = 'P'
          outvalunit = '?'
          scale_map = 1.0
        case (6) ! Passive scalar
          outval = 'passive'
          outvalunit = '?'
          scale_map = 1.0
        case (7) ! Temperature
          outval = 'T'
          outvalunit = 'K'
          scale_map = scale_temk
        case (8) ! speed of sound
          outval = 'cs'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (9) ! Jean's length
          outval = 'Lj'
          outvalunit = 'pc'
          scale_map = scale_lpc
        case (10) ! velocity dispersion
          outval = 'vdisp'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (11) ! planar velocity
          outval = 'v2D'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (12) ! 3D velocity
          outval = 'v3D'
          outvalunit = 'km/s'
          scale_map = scale_vkms
        case (13) ! 2D gravitational force
          outval = 'f2D'
          outvalunit = 'kpc/Myr^2'
          scale_map = scale_forc
          gravinp=.true.
        case (14) ! 3D-gravitational force
          outval = 'f3D'
          outvalunit = 'kpc/Myr^2'
          scale_map = scale_forc
          gravinp=.true.
        case (15) ! x-gravitational force
          outval = 'fx'
          outvalunit = 'kpc/Myr^2'
          scale_map = scale_forc
          gravinp=.true.
        case (16) ! y-gravitational force
          outval = 'fy'
          outvalunit = 'kpc/Myr^2'
          scale_map = scale_forc
          gravinp=.true.
        case (17) ! z-gravitational force
          outval = 'fz'
          outvalunit = 'kpc/Myr^2'
          scale_map = scale_forc
          gravinp=.true.
        case (18) ! tidal force (max eigenvalue)
          outval = 'tide'
          outvalunit = 'Myr^-2'
          scale_map = 1./(scale_tmyr**2)
          gravinp=.true.
        case (19) ! tidal force (max eigenvalue) normalized to gravity
          outval = 'tideg'
          outvalunit = ' '
          scale_map = 1.0
          gravinp=.true.
        case (20) ! shear
          outval = 'shear'
          outvalunit = 'Myr^-1'
          scale_map = 1./scale_tmyr
        case (21) ! shear normalized to self-gravity
          outval = 'shearg'
          outvalunit = ' '
          scale_map = 1.0
        case (22) ! compressive mode
          outval = 'compturb'
          outvalunit = 'km/s/pc'
          scale_map = scale_vkms / scale_lpc
        case (23) ! solenoidal mode
          outval = 'solturb'
          outvalunit = 'km/s/pc'
          scale_map = scale_vkms / scale_lpc
        case default ! density
          outval = 'rho'
          outvalunit = 'cm^-3'
          scale_map = scale_dhcc
      end select
    endif  

    out = '_'//TRIM(outval)//out
  endif

  if(cen) out = '_centered'//out

  if (lmax < 10) then
    write(nlev,'(I1)') lmax
  else
    write(nlev,'(I2)') lmax
  endif
  if(.not. p) out = '_l'//TRIM(nlev)//out


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! BEGIN PART ONLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(p) then  

    filename=TRIM(repository)//'/part_'//TRIM(nchar)//'.out00001'
    inquire(file=filename, exist=ok)
    if(.not. ok)then
      write(*,*) "Error: ", trim(filename), " not found"
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif
  
    if(sfr > 0.0) then
      sfrn = 1+floor(t / sfr)
      allocate(sfrhist(1:sfrn))
      sfrhist = 0.0
    endif

    if(part) then
      if(ns) then
        open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.ns')
      else
        open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.part')
      endif

      if(id) open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.partid')
    endif

    nselect = 0

    do k=1,ncpu ! Loop over CPU files
      write(ncharcpu,fmt='(I5.5)') k
      filename=TRIM(repository)//'/part_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
      open(unit=1, file=trim(filename), status='old', form='unformatted')
      read(1)
      read(1) ndim
      read(1) npart
      read(1)
      read(1)
      read(1)
      read(1)
      read(1)
        
      allocate(m(1:npart), x(1:npart,1:ndim), v(1:npart,1:ndim), part_id(1:npart))
      if(ns) allocate(age(1:npart))
                        
      do i=1,ndim
        read(1) m
        x(1:npart,i) = m/boxlen !m-0.5*boxlen
      end do 

      do i=1,ndim
        read(1) m
        v(1:npart,i) = m
      end do

      read(1) m
      read(1) part_id

      if(ns) then
        read(1)
        read(1) age
      end if


      ! Loop over particles
      do i=1, npart, pstep
        if(ns) then 
          agetmp = t-age(i) ! compute age from date of birth
          atmp = age(i) ! only used in the selection below no need to set for non-ns
        endif

        ! particle selection (id AND mass AND position AND (age IF ns) )
        if( part_id(i) >= 0 .AND. m(i) < pmmax .AND. m(i) > pmmin .AND. x(i,1) > xmin .AND. x(i,1) < xmax .AND. x(i,2) > ymin .AND. x(i,2) < ymax .AND. x(i,3) > zmin .AND. x(i,3) < zmax .AND. ( (ns .EQV. .false.) .OR. ( ns .EQV. .true. .AND. atmp > 0 .AND. agetmp > pamin .AND. agetmp < pamax ) ) ) then

          if(part) then
            if(ns) then
              write(3,fmt='(8(1x,E15.8))'), (x(i,1)-offsetx)*boxlen*scale_lkpc, (x(i,2)-offsety)*boxlen*scale_lkpc, (x(i,3)-offsetz)*boxlen*scale_lkpc, v(i,:)*scale_vkms, m(i)*scale_msun, agetmp*scale_tmyr
            else
              write(3,fmt='(7(1x,E15.8))'), (x(i,1)-offsetx)*boxlen*scale_lkpc, (x(i,2)-offsety)*boxlen*scale_lkpc, (x(i,3)-offsetz)*boxlen*scale_lkpc, v(i,:)*scale_vkms, m(i)*scale_msun
            endif
          endif

          if(id) then
            write(3,fmt='(I8.8)'), part_id(i)
          endif

          if(sfr > 0.0) then
            isfr = floor( age(i) / sfr ) +1

            if(eta_sn > 0.0 .AND. ( (f_w == 0.0 .AND. agetmp > agesn) .OR. (f_w > 0.0) ) ) then ! SN feedback
              ! kinetic feedback: "GMC" particles have been created at the time of birth and carry eta_sn of the mass. Correction must be made for all the stars.
              ! thermal feedback: the mass of the SNe is reduced at the time of the blast. Correction must be made for the stars older than 10 Myr only.
              sfrhist(isfr) = sfrhist(isfr) + m(i) / (1.-eta_sn)  ! mass_at_birth = mass_post_feedback / (1.-eta_sn)
            else
		      sfrhist(isfr) = sfrhist(isfr) + m(i)
            endif
          endif

          if(makemap) then ! store data for fits
            nselect = nselect+1
            xselect(nselect,1) = x(i,1)
            xselect(nselect,2) = x(i,2)
            xselect(nselect,3) = x(i,3)
            mselect(nselect) = m(i)
          endif

        end if ! end particle selection
      end do ! end loop over particles
    
      deallocate(m,x,v,part_id)
      if(ns) deallocate(age)

      close(1)  
    end do ! end loop over CPU files
 
!! SFR output  
    if(sfr > 0.0) then
      open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.sfr')
      do i=1, sfrn
        write(3,*) (i-1)*sfr*scale_tmyr, sfrhist(i)*scale_msun/(sfr*scale_tyr) ! left side of the bin [Myr], bin value [Msun/yr]
      end do
      deallocate(sfrhist)
    endif

    if(part .OR. ns .OR. id .OR. (sfr > 0.0) ) then
      close(3)
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif


!! particle map
    npart = nselect

    nxmap = pix
    nymap = INT(nxmap*(yymax-yymin)/(xxmax-xxmin))
    if(nymap > floor(nxmap*(yymax-yymin)/(xxmax-xxmin))) nymap = nymap+1

    allocate(map2(nxmap,nymap))
    map2 = 0.
  
    ! pixel size
    dx = (xxmax-xxmin)/nxmap
    dy = (yymax-yymin)/nymap

    do i=1, npart
      bx = floor((xselect(i,idim)-xxmin) / dx)+1
      by = floor((xselect(i,jdim)-yymin) / dy)+1

      map2(bx,by) = map2(bx,by) + mselect(i)*pstep ! simple non-smoothing version
    enddo

    map2 = map2 / (dx*dy * boxlen * boxlen) * scale_surf ! [Msun/pc^2]

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! END PART ONLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! BEGIN HYDRO/GRAV ONLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  else

!! open output files
    if(hydro) open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.hydro')
    if(rand > 0.0) open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.rand')
    if(grav) open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.grav')
    if(q) open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.q')

    if(pdf > 0) then
      lrhomin = log10(rhomin)
      lrhomax = log10(rhomax)
      lrhoampli = (lrhomax-lrhomin) / pdf
      allocate(pdfhist(1:pdf))
      pdfhist = 0.0
      typ = 1
    endif
  
!!! check for input
    if(gravinp) then
      filename=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out00001'
      inquire(file=filename, exist=ok)
      if(.not. ok)then
        write(*,*) "Error: ", trim(filename), " not found"
#ifdef MPIMODE
        call MPI_FINALIZE(mpi_ierror)
#endif
        stop
      endif
    endif
    
    filename=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out00001'
    inquire(file=filename, exist=ok)
    if(.not. ok)then
      write(*,*) "Error: ", trim(filename), " not found"
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif
  
    filename=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out00001'
    inquire(file=filename, exist=ok)
    if(.not. ok)then
      write(*,*) "Error: ", trim(filename), " not found"
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif
  
    open(unit=1, file=filename, status='old', form='unformatted')
    read(1) ncpu
    read(1) ndim
    read(1) nx, ny, nz
    read(1) nlevelmax
    read(1) ngridmax
    read(1) nboundary
    read(1) ngrid_current
    read(1) boxlen
    close(1)
  
    twotondim=2**ndim
    xbound=(/dble(nx/2),dble(ny/2),dble(nz/2)/)
  
    allocate(ngridfile(1:ncpu+nboundary,1:nlevelmax))
    allocate(ngridlevel(1:ncpu,1:nlevelmax))
    if(nboundary > 0) allocate(ngridbound(1:nboundary,1:nlevelmax))
  
    if(TRIM(ordering).eq.'hilbert') then ! if Hilbert
      dmax=max(xmax-xmin,ymax-ymin,zmax-zmin)
      do ilevel=1,lmax
        dx=0.5d0**ilevel
        if(dx.lt.dmax) exit
      end do
    
      bit_length=ilevel-1
      maxdom=2**bit_length
      imin=0; imax=0; jmin=0; jmax=0; kmin=0; kmax=0

      if(bit_length>0)then
        imin=int(xmin*dble(maxdom))
        imax=imin+1
        jmin=int(ymin*dble(maxdom))
        jmax=jmin+1
        kmin=int(zmin*dble(maxdom))
        kmax=kmin+1
      end if
  
      dkey=(dble(2**(nlevelmax+1)/dble(maxdom)))**ndim
      
      ndom=1
      if(bit_length>0)ndom=8
      idom(1)=imin
      idom(2)=imax
      idom(3)=imin
      idom(4)=imax
      idom(5)=imin
      idom(6)=imax
      idom(7)=imin
      idom(8)=imax
      
      jdom(1)=jmin
      jdom(2)=jmin
      jdom(3)=jmax
      jdom(4)=jmax
      jdom(5)=jmin
      jdom(6)=jmin
      jdom(7)=jmax
      jdom(8)=jmax
      
      kdom(1)=kmin
      kdom(2)=kmin
      kdom(3)=kmin
      kdom(4)=kmin
      kdom(5)=kmax
      kdom(6)=kmax
      kdom(7)=kmax
      kdom(8)=kmax
      
      do i=1,ndom
        if(bit_length>0)then
          call hilbert3d(idom(i),jdom(i),kdom(i),order_min,bit_length,1)
        else
          order_min=0.0d0
        endif
        bounding_min(i)=(order_min)*dkey
        bounding_max(i)=(order_min+1.0D0)*dkey
      end do
  
      cpu_min=0
      cpu_max=0
      do impi=1,ncpu
        do i=1,ndom
          if (bound_key(impi-1).le.bounding_min(i).and.bound_key(impi).gt.bounding_min(i)) cpu_min(i)=impi
          if (bound_key(impi-1).lt.bounding_max(i).and.bound_key(impi).ge.bounding_max(i)) cpu_max(i)=impi
        end do
      end do
       
      ncpu_read=0
      do i=1,ndom
        do j=cpu_min(i), cpu_max(i)
          if(.not. cpu_read(j))then
            ncpu_read=ncpu_read+1
            cpu_list(ncpu_read)=j
            cpu_read(j)=.true.
          endif
        enddo
      enddo
    else
      ncpu_read=ncpu
      do j=1,ncpu
        cpu_list(j)=j
      end do
    end  if ! end if Hilbert
  
    ! Compute hierarchy
    do ilevel=1,lmax
      nx_full=2**ilevel
      ny_full=nx_full
      nz_full=nx_full
      imin=int(xxmin*dble(nx_full))+1
      imax=int(xxmax*dble(nx_full))+1
      jmin=int(yymin*dble(ny_full))+1
      jmax=int(yymax*dble(ny_full))+1
      allocate(grid(ilevel)%map(imin:imax,jmin:jmax))
      allocate(grid(ilevel)%rho(imin:imax,jmin:jmax))
      grid(ilevel)%map(:,:)=0.0
      grid(ilevel)%rho(:,:)=0.0
      grid(ilevel)%imin=imin
      grid(ilevel)%imax=imax
      grid(ilevel)%jmin=jmin
      grid(ilevel)%jmax=jmax    
      grid(ilevel)%kmin=int(zzmin*dble(nz_full))+1
      grid(ilevel)%kmax=int(zzmax*dble(nz_full))+1
    end do
  
  
!! Read data
    ! Loop over cpu files
    do k=1,ncpu_read
      icpu=cpu_list(k)
      write(ncharcpu,'(I5.5)') icpu
      
      ! Open AMR file and skip header
      filename=TRIM(repository)//'/amr_'//TRIM(nchar)//'.out'//TRIM(ncharcpu)
      open(unit=1, file=filename, status='old', form='unformatted')
      do i=1,21
        read(1)
      end do
      ! Read grid numbers
      read(1) ngridlevel
      ngridfile(1:ncpu,1:nlevelmax)=ngridlevel
      read(1)
      if(nboundary>0)then
        read(1)
        read(1)
        read(1) ngridbound
        ngridfile(ncpu+1:ncpu+nboundary,1:nlevelmax)=ngridbound
      endif
      read(1)     
      read(1) ! comment this line for old stuff
      if(TRIM(ordering).eq.'bisection')then
        do i=1,5
          read(1)
        end do
      else
        read(1)
      endif
      read(1)
      read(1)
      read(1)
  
      ! Open HYDRO file and skip header
      open(unit=2, file=TRIM(repository)//'/hydro_'//TRIM(nchar)//'.out'//TRIM(ncharcpu), status='old', form='unformatted')
      read(2)
      read(2) nvarh
      read(2)
      read(2)
      read(2)
      read(2)
  
      if(gravinp) then
        ! Open GRAV file and skip header
        open(unit=7, file=TRIM(repository)//'/grav_'//TRIM(nchar)//'.out'//TRIM(ncharcpu), status='old', form='unformatted')
        read(7)
        read(7)
        read(7)
        read(7)
      endif
      
      ! Loop over levels
      do ilevel=1, lmax
        ! Geometry
        dx=0.5**ilevel
        dxline=1
        if(ndim==3) dxline=dx
        nx_full=2**ilevel
        ny_full=2**ilevel
        nz_full=2**ilevel
  
        do ind=1,twotondim
          iz=(ind-1)/4
          iy=(ind-1-4*iz)/2
          ix=(ind-1-2*iy-4*iz)
          xc(ind,1)=(dble(ix)-0.5D0)*dx
          xc(ind,2)=(dble(iy)-0.5D0)*dx
          xc(ind,3)=(dble(iz)-0.5D0)*dx
        end do
  
        ! Allocate work arrays
        ngrida=ngridfile(icpu,ilevel)
        grid(ilevel)%ngrid=ngrida
        if(ngrida>0)then
          allocate(xg(1:ngrida,1:ndim))
          allocate(son(1:ngrida,1:twotondim))
          allocate(var(1:ngrida,1:twotondim,1:nvarh))
          allocate(x  (1:ngrida,1:ndim))
          allocate(rho(1:ngrida))
          allocate(map(1:ngrida))
          allocate(ref(1:ngrida))
  
          if(typ == 10) then
            allocate(vdisp2(1:ngrida))
            allocate(vdisp1(1:ngrida))
            allocate(vdispm(1:ngrida))
          endif
  
          if(gravinp) then
            allocate(varg(1:ngrida,1:twotondim,1:ndim))
            if(typ == 18 .or. typ == 19) allocate(tt(1:ngrida,1:twotondim,1:3,1:3))
          endif
  
          if(typ == 20 .or. typ == 21) allocate(sh(1:ngrida,1:twotondim,1:3,1:3))
  
          if(typ == 22) allocate(nablav(1:ngrida,1:twotondim))
  
          if(typ == 23) then
            allocate(nablav1(1:ngrida,1:twotondim))
  	        allocate(nablav2(1:ngrida,1:twotondim))
            allocate(nablav3(1:ngrida,1:twotondim))
          endif
  
        endif
  
        ! Loop over domains
        do j=1,nboundary+ncpu
          ! Read AMR data
          if(ngridfile(j,ilevel)>0)then
            read(1) ! Skip grid index
            read(1) ! Skip next index
            read(1) ! Skip prev index
            ! Read grid center
            do iidim=1,ndim
              if(j.eq.icpu)then
                read(1) xg(:,iidim)
              else
                read(1)
              endif
            end do
            read(1) ! Skip father index
            do ind=1,2*ndim
              read(1) ! Skip nbor index
            end do
            ! Read son index
            do ind=1,twotondim
              if(j.eq.icpu)then
                read(1) son(:,ind)
              else
                read(1)
              end if
            end do
            ! Skip cpu map
            do ind=1,twotondim
              read(1)
            end do
            ! Skip refinement map
            do ind=1,twotondim
              read(1)
            end do
          endif
  
          if(gravinp) then
            ! Read grav data
            read(7)
            read(7)
            if(ngridfile(j,ilevel)>0)then
              ! Read grav variables
              do ind=1,twotondim
                do ivar=1,ndim
                  if(j.eq.icpu)then
                    read(7) varg(:,ind,ivar)
                  else
                    read(7)
                  end if
                end do
              end do
            end if
          endif
          
          
          ! Read hydro data
          read(2)
          read(2)
          if(ngridfile(j,ilevel)>0)then
            ! Read hydro variables
            do ind=1,twotondim
              do ivar=1,nvarh
                if(j.eq.icpu)then
                  read(2) var(:,ind,ivar)
                else
                  read(2)
                end if
              end do
            end do
          end if
        end do
        ! end loop over domains
        
        ! Compute map
        if(ngrida>0)then
  
          ! reset dispersion variables
          if(typ == 10) then
            vdisp2 = 0.0
            vdisp1 = 0.0
            vdispm = 0.0
          endif
  
          ! Loop over cells
          do ind=1,twotondim
            do i=1,ngrida
              ! Compute cell center
              x(i,1)=(xg(i,1)+xc(ind,1)-xbound(1))
              x(i,2)=(xg(i,2)+xc(ind,2)-xbound(2))
              if(ndim>2)x(i,3)=(xg(i,3)+xc(ind,3)-xbound(3))
              ref(i)=son(i,ind)>0.and.ilevel<lmax
  
              if(.not.ref(i)) then
                ix=int(x(i,idim)*dble(nx_full))+1
                iy=int(x(i,jdim)*dble(ny_full))+1
                iz=int(x(i,kdim)*dble(nz_full))+1
  
                ! volume selection
                if(ix>=grid(ilevel)%imin.and.iy>=grid(ilevel)%jmin.and.ix<=grid(ilevel)%imax.and.iy<=grid(ilevel)%jmax.and.iz>=grid(ilevel)%kmin.and.iz<=grid(ilevel)%kmax) then
  
                  if(makemap .or. q) then
                    ! Extract variable
                    select case (typ)
                      case (-1)
                        map(i) = icpu
                      case (0)
                        map(i) = ilevel
                      case (1) ! Density
                        map(i) = var(i,ind,1)
                      case (2) ! x-velocity
                        map(i) = var(i,ind,2)
                      case (3) ! y-velocity
                        map(i) = var(i,ind,3)
                      case (4) ! z-velocity
                        map(i) = var(i,ind,4)
                      case (5) ! Pressure
                        map(i) = var(i,ind,5)
                      case (6) ! Passive scalar
                        map(i) = var(i,ind,6)
                      case (7) ! Temperature
                        map(i) = var(i,ind,5)/var(i,ind,1)
                      case (8) ! Speed of sound = SQRT(gamma * P/rho)
                        map(i) = sqrt(5./3.* var(i,ind,5) / var(i,ind,1) )
                      case (9) ! Jeans length = sqrt(gamma * P / rho) * sqrt(pi / G / rho))   (* density to get the density-weighted value)
                        map(i) = sqrt(var(i,ind,5)) / var(i,ind,1) * 2.28823 ! * sqrt(5/3 * 3.14) = 2.28823
                      case (10) ! Velocity dispersion
                        if(ind==1) then ! compute the dispersion in the oct only once.
                          do j=1, twotondim
                            vdisp2(i) = vdisp2(i) + (var(i,j,2)**2+var(i,j,3)**2+var(i,j,4)**2)*var(i,j,1)
                            vdisp1(i) = vdisp1(i) + sqrt(var(i,j,2)**2+var(i,j,3)**2+var(i,j,4)**2)*var(i,j,1)
                            vdispm(i) = vdispm(i) + var(i,j,1)
                          end do
                        endif
                        map(i) = sqrt(abs(vdisp2(i) / vdispm(i) - (vdisp1(i) / vdispm(i))**2)) ! apply a different density weight for all cells of the oct
                      case (11) ! 2D velocity
                        map(i) = sqrt(var(i,ind,idim+1)**2+var(i,ind,jdim+1)**2)
                      case (12) ! 3D velocity
                        map(i) = sqrt(var(i,ind,2)**2+var(i,ind,3)**2+var(i,ind,4)**2)
                      case (13) ! 2D gravitational force
                        map(i) = sqrt(varg(i,ind,1)**2+varg(i,ind,2)**2)
                      case (14) ! 3D gravitational force
                        map(i) = sqrt(varg(i,ind,1)**2+varg(i,ind,2)**2+varg(i,ind,3)**2)
                      case (15) ! x-gravitational force
                        map(i) = varg(i,ind,1)
                      case (16) ! y-gravitational force
                        map(i) = varg(i,ind,2)
                      case (17) ! z-gravitational force
                        map(i) = varg(i,ind,3)
                      case (18) ! Tidal force
                        ! will be computed below
                      case (19) ! Tidal force, normalized to gravity
                        ! will be computed below
                      case (20) ! Shear
                        ! will be computed below
                      case (21) ! Shear, normalized to self-gravity
                        ! will be computed below
                      case (22) ! Compressive turbulence
                        nablav(i,ind) = 0
                        ! dvx/dx
                        nablav(i,ind) = nablav(i,ind) + var(i,2,2)-var(i,1,2)
                        nablav(i,ind) = nablav(i,ind) + var(i,4,2)-var(i,3,2)
                        nablav(i,ind) = nablav(i,ind) + var(i,6,2)-var(i,5,2)
                        nablav(i,ind) = nablav(i,ind) + var(i,8,2)-var(i,7,2)
                        ! dvy/dy
                        nablav(i,ind) = nablav(i,ind) + var(i,3,3)-var(i,1,3)
                        nablav(i,ind) = nablav(i,ind) + var(i,4,3)-var(i,2,3)
                        nablav(i,ind) = nablav(i,ind) + var(i,7,3)-var(i,5,3)
                        nablav(i,ind) = nablav(i,ind) + var(i,8,3)-var(i,6,3)
                        ! dvz/dz
                        nablav(i,ind) = nablav(i,ind) + var(i,5,4)-var(i,1,4)
                        nablav(i,ind) = nablav(i,ind) + var(i,6,4)-var(i,2,4)
                        nablav(i,ind) = nablav(i,ind) + var(i,7,4)-var(i,3,4)
                        nablav(i,ind) = nablav(i,ind) + var(i,8,4)-var(i,4,4)
      
                        nablav(i,ind) = nablav(i,ind) / 4. ! each axis was counted 4 times, average over cells in grid
                        map(i) = nablav(i,ind) / (2.*boxlen*dx)
                      case (23) ! Solenoidal turbulence
                        nablav1(i,ind) = 0
                        nablav1(i,ind) = nablav1(i,ind) + (var(i,3,4)-var(i,1,4)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) + (var(i,4,4)-var(i,2,4)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) + (var(i,7,4)-var(i,5,4)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) + (var(i,8,4)-var(i,6,4)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) - (var(i,5,3)+var(i,1,3)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) - (var(i,6,3)+var(i,2,3)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) - (var(i,7,3)+var(i,3,3)) / 4.
                        nablav1(i,ind) = nablav1(i,ind) - (var(i,8,3)+var(i,4,3)) / 4.
      
                        nablav2(i,ind) = 0
                        nablav2(i,ind) = nablav2(i,ind) + (var(i,5,2)-var(i,1,2)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) + (var(i,6,2)-var(i,2,2)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) + (var(i,7,2)-var(i,3,2)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) + (var(i,8,2)-var(i,4,2)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) - (var(i,2,4)+var(i,1,4)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) - (var(i,4,4)+var(i,3,4)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) - (var(i,6,4)+var(i,5,4)) / 4.
                        nablav2(i,ind) = nablav2(i,ind) - (var(i,8,4)+var(i,7,4)) / 4.
      
                        nablav3(i,ind) = 0
                        nablav3(i,ind) = nablav3(i,ind) + (var(i,5,3)-var(i,1,3)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) + (var(i,6,3)-var(i,2,3)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) + (var(i,7,3)-var(i,3,3)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) + (var(i,8,3)-var(i,4,3)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) - (var(i,2,2)+var(i,1,2)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) - (var(i,4,2)+var(i,3,2)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) - (var(i,6,2)+var(i,5,2)) / 4.
                        nablav3(i,ind) = nablav3(i,ind) - (var(i,8,2)+var(i,7,2)) / 4.
      
                        nablav1(i,ind) = sqrt(nablav1(i,ind)**2+nablav2(i,ind)**2+nablav3(i,ind)**2)
                        map(i) = nablav1(i,ind) /(2.*boxlen*dx)
                    end select
                  endif ! end makemap
        
  
                  if(typ == 18 .or. typ == 19) then ! tidal field
                    do j=1,3
                      if(mod(ind,2) == 0) then
                        tt(i,ind,1,j) = varg(i,ind,j)-varg(i,ind-1,j)
                      else
                        tt(i,ind,1,j) = varg(i,ind+1,j)-varg(i,ind,j)
                      endif
                      if(mod(ind,4) == 0 .or. mod(ind,4)==3 ) then
                        tt(i,ind,2,j) = varg(i,ind,j)-varg(i,ind-2,j)
                      else
                        tt(i,ind,2,j) = varg(i,ind+2,j)-varg(i,ind,j)
                      endif
                      if(ind/8.0 > 0.5) then
                        tt(i,ind,3,j) = varg(i,ind,j)-varg(i,ind-4,j)
                      else
                        tt(i,ind,3,j) = varg(i,ind+4,j)-varg(i,ind,j)
                      endif
                    enddo
                    tt(i,ind,:,:) = tt(i,ind,:,:) / (boxlen*dx)

                    ! put tidal tensor in diagonal form (see Renaud et al. 2008)
                    call eigenval(tt(i,ind,1,1),tt(i,ind,1,2),tt(i,ind,1,3),tt(i,ind,2,2),tt(i,ind,2,3),tt(i,ind,3,3),lambda1,lambda2,lambda3) ! warning: lambda's are only dF and NOT dF/dx.
                    if(typ == 18) then
                      map(i) = lambda1 / (boxlen*dx)
                    else
                      fgrav = sqrt(varg(i,ind,1)**2+varg(i,ind,2)**2+varg(i,ind,3)**2)
!                     map(i) = sqrt(lambda1**2+lambda2**2+lambda3**2) / fgrav ! ratio of tidal to total grav force
                      map(i) = lambda1 / fgrav ! ratio of tidal to total grav force
                    endif
                  endif ! end tides
  
  
                  if(typ == 20 .or. typ == 21) then ! shear
                    if(mod(ind,2) == 0) then ! right
                      sh(i,ind,1,2) = 0.5*(var(i,ind,3)-var(i,ind-1,3)) ! d(vy)/dx
                      sh(i,ind,1,3) = 0.5*(var(i,ind,4)-var(i,ind-1,4)) ! d(vz)/dx
                    else ! left
                      sh(i,ind,1,2) = 0.5*(var(i,ind+1,3)-var(i,ind,3)) ! d(vy)/dx
                      sh(i,ind,1,3) = 0.5*(var(i,ind+1,4)-var(i,ind,4)) ! d(vz)/dx
                    endif
                    if(mod(ind,4) == 0 .or. mod(ind,4)==3) then ! top
                      sh(i,ind,1,2) = sh(i,ind,1,2) + 0.5*(var(i,ind,2)-var(i,ind-2,2)) ! d(vx)/dy
                      sh(i,ind,2,3) = 0.5*(var(i,ind,4)-var(i,ind-2,4)) ! d(vz)/dy
                    else ! bottom
                      sh(i,ind,1,2) = sh(i,ind,1,2) + 0.5*(var(i,ind+2,2)-var(i,ind,2)) ! d(vx)/dy
                      sh(i,ind,2,3) = 0.5*(var(i,ind+2,4)-var(i,ind,4)) ! d(vz)/dy
                    endif
                    if(ind > 4) then ! front
                      sh(i,ind,1,3) = sh(i,ind,1,3) + 0.5*(var(i,ind,2)-var(i,ind-4,2)) ! d(vx)/dz
                      sh(i,ind,2,3) = sh(i,ind,2,3) + 0.5*(var(i,ind,3)-var(i,ind-4,3)) ! d(vy)/dz
                    else ! rear
                      sh(i,ind,1,3) = sh(i,ind,1,3) + 0.5*(var(i,ind+4,2)-var(i,ind,2)) ! d(vx)/dz
                      sh(i,ind,2,3) = sh(i,ind,2,3) + 0.5*(var(i,ind+4,3)-var(i,ind,3)) ! d(vy)/dz
                    endif
                    ! shear seen a centrifugal force = Omega^2 r
                    ! Omega = dv /dr = sh_from_above / dr 
                    ! centrifugal / selfgrav = Omega^2 / 4/3 pi G rho
                    ! centrifugal / fgrav = Omega^2 r / fgrav = (sh_from_above^2 / r) / fgrav 
                    sh(i,ind,:,:) = sh(i,ind,:,:)/(boxlen*dx)  ! dv / dr = Omega

                    if(typ == 20) then
                      map(i) = sqrt(sh(i,ind,1,2)**2+sh(i,ind,1,3)**2+sh(i,ind,2,3)**2) ! Omega
                    else
                      map(i) = (sh(i,ind,1,2)**2+sh(i,ind,1,3)**2+sh(i,ind,2,3)**2)/4.18879 / var(i,ind,1) ! omega^2 / (4/3 pi G rho) (all in Ramses units because ratio is dimensionless)
                    endif
                  endif ! end shear
  
  
                  if(var(i,ind,1) > rhomin .AND. var(i,ind,1) < rhomax) then ! density selection (ascii modes only)
  
                    mass=var(i,ind,1)*(boxlen*dx)**3

                    if(q) then
                      ! q output: x, y, z, map
                      write(3, '(4e20.6)') (x(i,1)-offsetx)*boxlen*scale_lkpc, (x(i,2)-offsety)*boxlen*scale_lkpc, (x(i,3)-offsetz)*boxlen*scale_lkpc, map(i)*scale_map
                    endif

  
                    if(hydro) then
                      ! normal output: x, y, z, vx, vy, vz, rho, lev, mass, T
                      write(3, '(7e20.6,1I5,2e20.6)') (x(i,1)-offsetx)*boxlen*scale_lkpc, (x(i,2)-offsety)*boxlen*scale_lkpc, (x(i,3)-offsetz)*boxlen*scale_lkpc, var(i,ind,2)*scale_vkms, var(i,ind,3)*scale_vkms, var(i,ind,4)*scale_vkms, var(i,ind,1)*scale_dhcc, ilevel, mass*scale_msun, var(i,ind,5)/var(i,ind,1)*scale_temk
                    endif
  
                    if(rand > 0.0) then
                      npartcell = floor(mass/rand)
                      do j=1, npartcell
                        call random_number(rx)
                        call random_number(ry)
                        call random_number(rz)
                        rx = (rx-0.5)*dx
                        ry = (ry-0.5)*dx
                        rz = (rz-0.5)*dx
                        ! random mode output: x, y, z
                        write(3, '(3e20.6)') (x(i,1)-offsetx+rx)*boxlen*scale_lkpc,(x(i,2)-offsety+ry)*boxlen*scale_lkpc,(x(i,3)-offsetz+rz)*boxlen*scale_lkpc
                      enddo
                    endif
  
  
                    if(pdf > 0) then
                      ipdf = int( (log10(var(i,ind,1))-lrhomin) / lrhoampli )+1
                      pdfhist(ipdf) = pdfhist(ipdf) + mass
                    endif
  
                    if(grav) then
                      ! hydro output + gravitational force: x, y, z, vx, vy, vz, rho, lev, mass, T, fx, fy, fz
                      write(3, '(7e20.6,1I5,5e20.6)') (x(i,1)-offsetx)*boxlen*scale_lkpc, (x(i,2)-offsetx)*boxlen*scale_lkpc, (x(i,3)-offsetx)*boxlen*scale_lkpc, var(i,ind,2)*scale_vkms, var(i,ind,3)*scale_vkms, var(i,ind,4)*scale_vkms, var(i,ind,1)*scale_dhcc, ilevel, mass*scale_msun, var(i,ind,5)/var(i,ind,1)*scale_temk, varg(i,ind,1)*scale_forc, varg(i,ind,2)*scale_forc, varg(i,ind,3)*scale_forc
                    endif
  
                  endif ! end density selection
  
  
                  if(makemap)then
                    !!! compute the projection
                    if(ndim==3)then
                      weight=(min(x(i,kdim)+dx/2.,zzmax)-max(x(i,kdim)-dx/2.,zzmin))/dx
                      weight=min(1.0d0,max(weight,0.0d0))
                    else
                      weight=1.0
                    endif
  
                    if(maxval)then
                      if(grid(ilevel)%map(ix,iy)<map(i))then
                        grid(ilevel)%map(ix,iy)=map(i) ! update the variable map
                        grid(ilevel)%rho(ix,iy)=var(i,ind,1)  ! not important ! update the weight map with the density at *this* position
                      endif
                    else
                      if(maxrho)then
                        if(grid(ilevel)%rho(ix,iy)<var(i,ind,1))then
                          grid(ilevel)%map(ix,iy)=map(i) ! update the variable map
                          grid(ilevel)%rho(ix,iy)=var(i,ind,1) ! update the weight map with the density at *this* position
                        endif
                      else
                        if(sum) then
                          grid(ilevel)%map(ix,iy)=grid(ilevel)%map(ix,iy)+map(i)
                          grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+var(i,ind,1)
                        else ! density weighted average
                          grid(ilevel)%map(ix,iy)=grid(ilevel)%map(ix,iy)+map(i)*var(i,ind,1)*dxline*weight/(zzmax-zzmin)
                          grid(ilevel)%rho(ix,iy)=grid(ilevel)%rho(ix,iy)+var(i,ind,1)*dxline*weight/(zzmax-zzmin)
                        endif
                      endif
                    endif
                  endif ! endif makemap
  
      
                endif ! end of volume selection   
              end if ! endif refined
            end do
          end do ! End loop over cell
  
  
          deallocate(xg, son, var, ref, rho, map, x)
  
  
          if(typ == 10) deallocate(vdisp2, vdisp1, vdispm)
  
          if(gravinp) then
            deallocate(varg)
            if(typ == 18 .or. typ == 19) deallocate(tt)
          endif
  
          if(typ == 20 .or. typ == 21) deallocate(sh)
  
          if(typ == 22) deallocate(nablav)
  
          if(typ == 23) then
            deallocate(nablav1)
  	      deallocate(nablav2)
  	      deallocate(nablav3)
          endif
  
        end if
      end do ! End loop over levels
  
      close(1)
      close(2)
      if(grav) close(7)
    end do ! End loop over cpu
  
  
  
!!! close ascii file
    if(hydro .or. grav) then
      close(3)
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif
  
  
!!! PDF output  
    if(pdf > 0) then
      open(3, file='rdr_'//TRIM(nchar)//TRIM(out)//'.pdf')
      do ipdf=1, pdf
        write(3,*) scale_dhcc*10**( (ipdf-1) * lrhoampli + lrhomin ), scale_msun*pdfhist(ipdf) ! left side of the bin, bin value
      end do
      deallocate(pdfhist)
      close(3)
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif
  
  
    if(makemap)then
      !!! compute the average for projections
      nx_full=2**lmax
      ny_full=nx_full
      imin=int(xxmin*dble(nx_full))+1
      imax=int(xxmax*dble(nx_full))
      jmin=int(yymin*dble(ny_full))+1
      jmax=int(yymax*dble(ny_full))
  
      do ix=imin,imax
        xmin=((ix-0.5)/2**lmax)
        do iy=jmin,jmax
          ymin=((iy-0.5)/2**lmax)
          do ilevel=1,lmax-1
            ndom=2**ilevel
            i=int(xmin*ndom)+1
            j=int(ymin*ndom)+1
            if(maxval) then
              if(grid(lmax)%map(ix,iy)<grid(ilevel)%map(i,j))then
                grid(lmax)%map(ix,iy)=grid(ilevel)%map(i,j) ! update the variable map
                grid(lmax)%rho(ix,iy)=grid(ilevel)%rho(i,j) ! update the weight map with the density at *this* position
              endif
            else
              if(maxrho)then
                if(grid(lmax)%rho(ix,iy)<grid(ilevel)%rho(i,j))then
                  grid(lmax)%map(ix,iy)=grid(ilevel)%map(i,j) ! update the variable map
                  grid(lmax)%rho(ix,iy)=grid(ilevel)%rho(i,j) ! update the weight map with the density at *this* position
                endif
              else ! average or sum
                grid(lmax)%map(ix,iy)=grid(lmax)%map(ix,iy) + grid(ilevel)%map(i,j)
                grid(lmax)%rho(ix,iy)=grid(lmax)%rho(ix,iy) + grid(ilevel)%rho(i,j)
              endif
            endif
          end do
        end do
      end do
  
      if(nx_sample==0)then
        nxmap=imax-imin+1 ! projected average density
        nymap=jmax-jmin+1
        allocate(tmpmap(nxmap,nymap))
        if(maxval .or. maxrho .or. sum) then
          tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)
        else ! average
          tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)/grid(lmax)%rho(imin:imax,jmin:jmax)
        endif    
      else ! this is never done ?  (nx_sample = 0 and is never changed)
        if(ny_sample==0) ny_sample = (jmax-jmin+1)/(imax-imin+1) !nx_sample
        allocate(tmpmap(0:nx_sample,0:ny_sample))
        nxmap=nx_sample+1
        nymap=ny_sample+1
        do i=0,nx_sample
          ix=int(dble(i)/dble(nx_sample)*dble(imax-imin+1))+imin
          ix=min(ix,imax)
          do j=0,ny_sample
            iy=int(dble(j)/dble(ny_sample)*dble(jmax-jmin+1))+jmin
            iy=min(iy,jmax)
            if(maxval .or. maxrho .or. sum) then
              tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)
            else ! average
              tmpmap=grid(lmax)%map(imin:imax,jmin:jmax)/grid(lmax)%rho(imin:imax,jmin:jmax)
            endif    
          end do
        end do
      endif
      allocate(map2(nxmap,nymap))
      map2 = tmpmap*scale_map
    endif ! end makemap
  
  
!!! PSD
#ifdef FFTW
    if(psd) then
      allocate(psdin(nxmap,nymap))
      allocate(psdout(nxmap,nymap))
      psdnorm = 1.0D0/dble(nxmap*nymap)**2
      do j=1, nymap
        do i=1, nxmap
          psdin(i,j) = map2(i,j)
          !psdin(i,j) = map2(i,j)*(1-cos(2.*pi/nymap*j))*(1-cos(2.*pi/nxmap*i)) ! apodization fucntion (because the signal is not periodic)
        end do
      end do
  
      psdmap = 1+ int(sqrt(real(nxmap/2+1)**2+real(nymap/2+1)**2))
      allocate(ppsd(psdmap))
      allocate(npsd(psdmap))
      ppsd = 0.0
      npsd = 0  
      call dfftw_plan_dft_2d(plan_forward, nxmap, nymap, psdin, psdout, FFTW_FORWARD, FFTW_ESTIMATE)
      call dfftw_execute(plan_forward)
      call dfftw_destroy_plan(plan_forward)
  
      fpsd = 1.0D0/((xxmax-xxmin)*boxlen)  ! step in frequency = sampling frequency / N   (same for both axes)
  
      write(charlmax,'(I2.2)') lmax
      if(psdm) then
        open(2, file='psdm_'//TRIM(nchar)//'_l'//TRIM(charlmax)//TRIM(out)//'.psdm')
#ifdef MPIMODE
        if(mpi_rank == 0) open(11, file='psdm_l'//TRIM(charlmax)//TRIM(out)//'.wavenumber')
#endif
      endif
  
      do j=0, nymap/2
        do i=0, nxmap-1
          binpsd = 1 + int(sqrt(real(min(i,nxmap-i))**2 + real(j)**2))
          ppsd(binpsd) = ppsd(binpsd) + 4.*abs(psdout(i+1,j+1))**2 ! sum power
          npsd(binpsd) = npsd(binpsd) + 2
          if(psdm) then
#ifdef MPIMODE
            write(2,*) abs(psdout(i+1,j+1))**2*psdnorm
            if(mpi_rank == 0) write(11,*) binpsd*fpsd/scale_lpc, scale_lpc/(binpsd*fpsd)
#else
            write(2,*) binpsd*fpsd/scale_lpc, abs(psdout(i+1,j+1))**2*psdnorm, scale_lpc/(binpsd*fpsd)
#endif
  		endif
        end do
      end do
  
      if(psdm) then
        close(2)
#ifdef MPIMODE
        if(mpi_rank == 0) close(11)
#endif
      endif
  
      if(.not. psdm) then
        open(1, file='psd_'//TRIM(nchar)//'_l'//TRIM(charlmax)//TRIM(out)//'.psd')
        do i=0, psdmap-1
          write(1,*) i*fpsd/scale_lpc, ppsd(i+1)*psdnorm/npsd(i+1), scale_lpc/(i*fpsd)
        end do
        close(1)
      endif
      
#ifdef MPIMODE
      call MPI_FINALIZE(mpi_ierror)
#endif
      stop
    endif ! end psd
#endif
! end FFTW

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! END HYDRO ONLY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! write fits (grid or particles)
if(makemap) then
#ifdef FITS
    status=0

    filename = 'rdr_'//TRIM(nchar)//TRIM(out)//'.fits'
    call deletefile(filename,status)
    call ftgiou(unit,status)
    blocksize=1
    call ftinit(unit,filename,blocksize,status)
    simple=.true.
    bitpix=-32
    naxis=2
    naxes(1)=nxmap
    naxes(2)=nymap

    extend=.true.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  
    call ftpkyd(unit,'time',t*scale_tmyr,6,'time',status)
    call ftpkyd(unit,'boxlen',boxlen,6,'boxlen',status)
    call ftpkyj(unit,'idim',idim,'idim',status)
    call ftpkyj(unit,'jdim',jdim,'jdim',status)
    call ftpkyj(unit,'kdim',kdim,'kdim',status)
    call ftpkyd(unit,'xmin',xxmin*scale_lkpc,6,'xmin',status)
    call ftpkyd(unit,'xmax',xxmax*scale_lkpc,6,'xmax',status)
    call ftpkyd(unit,'ymin',yymin*scale_lkpc,6,'ymin',status)
    call ftpkyd(unit,'ymax',yymax*scale_lkpc,6,'ymax',status)
    call ftpkyd(unit,'zmin',zzmin*scale_lkpc,6,'zmin',status)
    call ftpkyd(unit,'zmax',zzmax*scale_lkpc,6,'zmax',status)
    call ftpkyj(unit,'lmax',lmax,'lmax',status)
    call ftpkys(unit,'outval',outval,'value',status)
    call ftpkys(unit,'outvalunit',outvalunit,'value unit',status)
    call ftpkyl(unit,'maxval',maxval,'maxval',status)
    call ftpkyl(unit,'maxrho',maxrho,'maxrho',status)
  
    group=1
    fpixel=1
    nelements=naxes(1)*naxes(2)
    call ftppre(unit,group,fpixel,nelements,MAP2,status)
    call ftclos(unit, status)
    call ftfiou(unit, status)
#endif
  write(*,*) TRIM(filename)//' has been created.'
  endif ! end makemap

!!! terminate mpi
#ifdef MPIMODE
  call MPI_FINALIZE(mpi_ierror)
#endif

end program rdramses

!=======================================================================
!=======================================================================
!=======================================================================

#ifdef FITS
subroutine deletefile(filename,status) !  Delete a FITS file

  integer::status,unit,blocksize
  character(*)::filename
  
  if (status .gt. 0) return

  call ftgiou(unit,status) ! Get an unused Logical Unit Number
  call ftopen(unit,filename,1,blocksize,status) ! Try to open the file
  
  if (status .eq. 0)then ! file is opened: delete it 
    call ftdelt(unit,status)
  else if (status .eq. 103)then ! file doesn't exist: reset status and clear errors
    status=0
    call ftcmsg
  else ! there was some other error opening the file: delete the file anyway
    status=0
    call ftcmsg
    call ftdelt(unit,status)
  end if
  
  call ftfiou(unit, status) ! Free the unit number

end
#endif


!=======================================================================
!=======================================================================
!=======================================================================

subroutine hilbert3d(x,y,z,order,bit_length,npoint)
  implicit none

  integer,intent(in)::bit_length, npoint
  integer,intent(in),dimension(1:npoint)::x, y, z
  real(kind=8),intent(out),dimension(1:npoint)::order

  logical,dimension(0:3*bit_length-1)::i_bit_mask
  logical,dimension(0:1*bit_length-1)::x_bit_mask, y_bit_mask, z_bit_mask
  integer,dimension(0:7,0:1,0:11)::state_diagram
  integer::i, ip, cstate, nstate, b0, b1, b2, sdigit, hdigit

  if(bit_length>bit_size(bit_length))then
    write(*,*)'Maximum bit length=',bit_size(bit_length)
    write(*,*)'stop in hilbert3d'
    stop
  endif

  state_diagram = RESHAPE( (/   1, 2, 3, 2, 4, 5, 3, 5,&
                            &   0, 1, 3, 2, 7, 6, 4, 5,&
                            &   2, 6, 0, 7, 8, 8, 0, 7,&
                            &   0, 7, 1, 6, 3, 4, 2, 5,&
                            &   0, 9,10, 9, 1, 1,11,11,&
                            &   0, 3, 7, 4, 1, 2, 6, 5,&
                            &   6, 0, 6,11, 9, 0, 9, 8,&
                            &   2, 3, 1, 0, 5, 4, 6, 7,&
                            &  11,11, 0, 7, 5, 9, 0, 7,&
                            &   4, 3, 5, 2, 7, 0, 6, 1,&
                            &   4, 4, 8, 8, 0, 6,10, 6,&
                            &   6, 5, 1, 2, 7, 4, 0, 3,&
                            &   5, 7, 5, 3, 1, 1,11,11,&
                            &   4, 7, 3, 0, 5, 6, 2, 1,&
                            &   6, 1, 6,10, 9, 4, 9,10,&
                            &   6, 7, 5, 4, 1, 0, 2, 3,&
                            &  10, 3, 1, 1,10, 3, 5, 9,&
                            &   2, 5, 3, 4, 1, 6, 0, 7,&
                            &   4, 4, 8, 8, 2, 7, 2, 3,&
                            &   2, 1, 5, 6, 3, 0, 4, 7,&
                            &   7, 2,11, 2, 7, 5, 8, 5,&
                            &   4, 5, 7, 6, 3, 2, 0, 1,&
                            &  10, 3, 2, 6,10, 3, 4, 4,&
                            &   6, 1, 7, 0, 5, 2, 4, 3 /), &
                            & (/8 ,2, 12 /) )

  do ip=1,npoint
    ! convert to binary
    do i=0,bit_length-1
      x_bit_mask(i)=btest(x(ip),i)
      y_bit_mask(i)=btest(y(ip),i)
      z_bit_mask(i)=btest(z(ip),i)
    enddo

    ! interleave bits
    do i=0,bit_length-1
      i_bit_mask(3*i+2)=x_bit_mask(i)
      i_bit_mask(3*i+1)=y_bit_mask(i)
      i_bit_mask(3*i  )=z_bit_mask(i)
    end do

    ! build Hilbert ordering using state diagram
    cstate=0
    do i=bit_length-1,0,-1
      b2=0
      if(i_bit_mask(3*i+2))b2=1
      b1=0
      if(i_bit_mask(3*i+1))b1=1
      b0=0
      if(i_bit_mask(3*i  ))b0=1
      
      sdigit=b2*4+b1*2+b0
      nstate=state_diagram(sdigit,0,cstate)
      hdigit=state_diagram(sdigit,1,cstate)
      i_bit_mask(3*i+2)=btest(hdigit,2)
      i_bit_mask(3*i+1)=btest(hdigit,1)
      i_bit_mask(3*i  )=btest(hdigit,0)
      cstate=nstate
    enddo

    ! save Hilbert key as double precision real
    order(ip)=0.
    do i=0,3*bit_length-1
      b0=0 ; if(i_bit_mask(i))b0=1
      order(ip)=order(ip)+dble(b0)*dble(2)**i
    end do
  end do

end subroutine hilbert3d

!=======================================================================
!=======================================================================
!=======================================================================

subroutine eigenval(t0,t1,t2,t3,t4,t5,lambda1,lambda2,lambda3)
  implicit none
  real(kind=8)::c0,c1,c2,p,q,r,s,d,phi,cphi,sphi,tmp
  real(kind=8),intent(in)::t0,t1,t2,t3,t4,t5
  real(kind=8),intent(out)::lambda1,lambda2,lambda3

  c2 = -t0-t3-t5
  c1 = t0*t3 + t0*t5 + t3*t5 - t1**2 - t2**2 - t4**2
  c0 = t0*t4**2 + t3*t2**2 + t5*t1**2 - t0*t3*t5 - 2*t2*t1*t4
	
  p = c2**2 - 3.*c1
  q = -27./2.0*c0 - c2**3 + 9./2.0*c2*c1
	
  if(p<0) then
    r=0
  else
    r = sqrt(p)/3.0
  endif
	
  s = -c2/3.0
  d = p**3 - q**2

  if(q == 0) then
    phi = 0.523598 ! PI /6.0
  else
    if (d<0) then
      phi = 0.0
    else	
      phi = 1./3. * atan(sqrt(d)/q)
    endif
  endif
	
  cphi = cos(phi)
  sphi = sin(phi)
  lambda1 = r * 2.*cphi + s
  lambda2 = r * (-cphi - sqrt(3.)*sphi) + s
  lambda3 = r * (-cphi + sqrt(3.)*sphi) + s

  if(lambda3 > lambda2) then
    tmp = lambda3
    lambda3 = lambda2
    lambda2 = tmp
  endif
  if(lambda2 > lambda1) then
    tmp = lambda2
    lambda2 = lambda1
    lambda1 = tmp
  endif
  if(lambda3 > lambda2) then
    tmp = lambda3
    lambda3 = lambda2
    lambda2 = tmp
  endif

end subroutine eigenval

