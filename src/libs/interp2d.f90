MODULE mo_interp2d


#ifdef ECRADOFF_DP
  use, intrinsic :: iso_fortran_env, realk => real64,     intk => int32
  use, intrinsic :: iso_c_binding, c_realk => c_double, c_intk => c_int32_t
#else
  use, intrinsic :: iso_fortran_env, realk => real32,    intk => int32
  use, intrinsic :: iso_c_binding, c_realk => c_float, c_intk => c_int32_t
#endif

  implicit none
  private
  public :: interp_2d

contains

  SUBROUTINE throw_error(error_string)
    character(*), intent(in) :: error_string
    write(*,*) error_string
  END SUBROUTINE

  FUNCTION get_dist(x1, x2) RESULT(dist)
  ! periodic-safe difference

  real(realk), intent(in) :: x1, x2
  real(realk) :: xmax, xmin, dist

  xmax = max(x1,x2)
  xmin = min(x1,x2)
  
  dist = xmax - xmin
  ! Takes the shortest path around th globe.
  if (dist > 180) then
    dist = xmin+360 - xmax
  endif  

  END FUNCTION get_dist

  FUNCTION get_geodist(dphi, phim, dlam) RESULT(geodist)

  real(realk), intent(in) :: dphi, phim, dlam
  real(realk), parameter :: pi=acos(-1.0)
  real(realk) :: geodist

  geodist = sqrt((sin(dphi*pi/180)+cos(dlam*pi/180))**2 + (cos(phim*2*pi/180)*sin(dlam*pi/180)**2))
  
  END FUNCTION get_geodist

  SUBROUTINE get_interp_coeffs_reduced(ysrc, xsrc, ny_src, redpts_src, tgtx, tgty, &
                                     & xysrc, xyidxsrc, atol_in)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

    integer(intk), intent(in), value :: ny_src

      real(realk), intent(in), dimension(ny_src) :: ysrc
      real(realk), intent(in), dimension(ny_src) :: xsrc  
    integer(intk), intent(in), dimension(ny_src) :: redpts_src

    real(realk), intent(in) :: tgtx, tgty
    
    real(realk), intent(out), dimension(8) :: xysrc
    integer(intk), intent(out), dimension(4) :: xyidxsrc


    real(realk), intent(in), optional :: atol_in
    real(realk) :: atol = 1.e-6

    integer(intk) :: jlow, jmid, jhigh, ilow1, ihigh1, ilow2, ihigh2
    integer(intk) :: jj
    logical :: foundy
    logical :: yisascending

    if (present(atol_in)) then
      atol = atol_in
    endif

    if (ysrc(ny_src) > ysrc(1)) then
      yisascending = .true.
    else
      yisascending = .false.
    endif

    foundy = .true.

    if (yisascending .and. (tgty < ysrc(1))) then
      jlow = 1
      jhigh = 1
    else if (yisascending .and. (tgty > ysrc(ny_src))) then
      jlow = ny_src
      jhigh = ny_src
    else if ((.not. yisascending) .and. (tgty < ysrc(ny_src))) then
      jlow = ny_src
      jhigh = ny_src
    else if ((.not. yisascending) .and. (tgty > ysrc(1))) then
      jlow = 1
      jhigh = 1
    else
      foundy = .false.
    endif

    if (.not. foundy) then
      !write(*,*) "Searching for y=",tgty
      if (yisascending) then
        jlow=1
        jhigh=ny_src
      else
        jlow=ny_src
        jhigh=1
      endif
      do while (abs(jhigh - jlow) > 1)
        jmid = (jlow + jhigh)/2
        if (abs(tgty - ysrc(jmid)) < atol) then
          jlow = jmid
          jhigh = jmid
        else if (tgty > ysrc(jmid)) then
          jlow = jmid
        else if (tgty < ysrc(jmid)) then
          jhigh = jmid
        endif
      end do
    endif

    !write(*,*) "Found y=",tgty," jlow=",jlow," and jhigh=",jhigh

    !xmax1 = 360 - 360./redpts_src(jlow) + xsrc(jlow)
    !xmax2 = 360 - 360./redpts_src(jhigh) + xsrc(jlow)
    ilow1 = mod(floor((tgtx-xsrc(jlow))/360.*redpts_src(jlow)),redpts_src(jlow))+1
    if (get_dist(tgtx, xsrc(jlow)+360/redpts_src(jlow)*(ilow1-1))<atol) then
      ihigh1 = ilow1
    else
      ihigh1 = mod(ilow1,redpts_src(jlow))+1
    endif
    ilow2 = mod(floor((tgtx-xsrc(jhigh))/360.*redpts_src(jhigh)),redpts_src(jhigh))+1
    if (get_dist(tgtx, xsrc(jhigh)+360/redpts_src(jhigh)*(ilow2-1))<atol) then
      ihigh2 = ilow2
    else
      ihigh2 = mod(ilow2,redpts_src(jhigh))+1
    endif

    ! A11 A12 A21 A22
    xyidxsrc(:) = 0
    do jj=1,jlow-1
      xyidxsrc(1) = xyidxsrc(1) + redpts_src(jj)
      xyidxsrc(2) = xyidxsrc(2) + redpts_src(jj)
    end do
    xyidxsrc(1) = xyidxsrc(1) + ilow1
    xyidxsrc(2) = xyidxsrc(2) + ihigh1
    do jj=1,jhigh-1
      xyidxsrc(3) = xyidxsrc(3) + redpts_src(jj)
      xyidxsrc(4) = xyidxsrc(4) + redpts_src(jj)
    end do
    xyidxsrc(3) = xyidxsrc(3) + ilow2
    xyidxsrc(4) = xyidxsrc(4) + ihigh2
    
    ! x11 x12 x21 x22 y11 y12 y21 y22
    xysrc(1) = xsrc(jlow) +360./redpts_src(jlow )*(ilow1 -1)
    xysrc(2) = xsrc(jlow) +360./redpts_src(jlow )*(ihigh1-1)
    xysrc(3) = xsrc(jhigh)+360./redpts_src(jhigh)*(ilow2 -1)
    xysrc(4) = xsrc(jhigh)+360./redpts_src(jhigh)*(ihigh2-1)
                
    xysrc(5) = ysrc(jlow)
    xysrc(6) = ysrc(jlow)
    xysrc(7) = ysrc(jhigh)
    xysrc(8) = ysrc(jhigh)

  END SUBROUTINE get_interp_coeffs_reduced

  SUBROUTINE get_interp_coeffs_lonlat(xsrc, ysrc, nx_src, ny_src, tgtx, tgty, &
                                    & xysrc, xyidxsrc, lat_row, atol_in)
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan

    integer(intk), intent(in), value :: ny_src, nx_src

    real(realk), intent(in), dimension(nx_src) :: xsrc  
    real(realk), intent(in), dimension(ny_src) :: ysrc

    real(realk), intent(in) :: tgtx, tgty
    real(realk), intent(out), dimension(8) :: xysrc
    integer(intk), intent(out), dimension(4) :: xyidxsrc

    logical, intent(in), value :: lat_row

    real(realk), intent(in), optional :: atol_in
    real(realk) :: atol = 1.e-6

    integer(intk) :: i_min, i_max, j_min, j_max

    integer(intk) :: imid, ilow, ihigh, jmid, jlow, jhigh
    logical :: foundx, foundy
    logical :: xisascending, yisascending
    logical :: xisperiodic

    ! Characterise grid ordering 
    xisascending = (xsrc(nx_src) > xsrc(1))
    yisascending = (ysrc(ny_src) > ysrc(1))

    ! Characterise periodicity
    xisperiodic = lat_row
    !yisperiodic = .not. lat_row

    if (present(atol_in)) then
      atol = atol_in
    endif

    foundx = .true.
    foundy = .true.

    ! Look at boundaries
    if (xisascending) then
      i_min = 1
      i_max = nx_src
    else
      i_min = nx_src
      i_max = 1
    endif

    if (yisascending) then
      j_min = 1
      j_max = ny_src
    else
      j_min = ny_src
      j_max = 1
    endif

    ! Triage
    if (xisperiodic) then ! x is longitude (periodic), y is latitude (clip)
      if ( (tgtx < xsrc(i_min)) .or. (tgtx > xsrc(i_max)) ) then
        ilow = i_max
        ihigh = i_min
      else
        foundx = .false.
      endif

      if (tgty < ysrc(j_min)) then
        jlow = j_min
        jhigh = j_min
      else if (tgty > ysrc(j_max)) then
        jlow = j_max
        jhigh = j_max
      else
        foundy = .false.
      endif
    else ! y is longitude (periodic)
      if (tgtx < xsrc(i_min)) then
        ilow = i_min
        ihigh = i_min
      else if (tgtx > xsrc(i_max)) then
        ilow = i_max
        ihigh = i_max
      else
        foundx = .false.
      endif

      if ( (tgty < ysrc(j_min)) .or. (tgty > ysrc(j_max)) ) then
        jlow = j_max
        jhigh = j_min
      else
        foundy = .false.
      endif
    endif

    ! Binary search inside domain

    if (.not. foundx) then
      if (xisascending) then
        ilow=1
        ihigh=nx_src
      else
        ilow=nx_src
        ihigh=1
      endif
      do while (abs(ihigh - ilow) > 1)
        imid = (ilow + ihigh)/2

        if (abs(tgtx - xsrc(imid)) < atol) then
          ilow = imid
          ihigh = imid
        else if (tgtx > xsrc(imid)) then
          ilow = imid
        else if (tgtx < xsrc(imid)) then
          ihigh = imid
        endif
      end do
    endif

    if (.not. foundy) then
      if (yisascending) then
        jlow=1
        jhigh=ny_src
      else
        jlow=ny_src
        jhigh=1
      endif
      do while (abs(jhigh - jlow) > 1)
        jmid = (jlow + jhigh)/2
        if (abs(tgty - ysrc(jmid)) < atol) then
          jlow = jmid
          jhigh = jmid
        else if (tgty > ysrc(jmid)) then
          jlow = jmid
        else if (tgty < ysrc(jmid)) then
          jhigh = jmid
        endif
      end do
    endif
    

  !write(*,*) "ilow=",ilow," ihigh=",ihigh
  !write(*,*) "jlow=",jlow," jhigh=",jhigh
  !yx 11
  if (lat_row) then
    xyidxsrc(1) = (jlow-1)*nx_src+ilow
    ! yx 12
    xyidxsrc(2) = (jlow-1)*nx_src+ihigh
    ! yx 21
    xyidxsrc(3) = (jhigh-1)*nx_src+ilow
    ! yx 22
    xyidxsrc(4) = (jhigh-1)*nx_src+ihigh
  else
    xyidxsrc(1) = (ilow-1)*ny_src+jlow
    xyidxsrc(2) = (ihigh-1)*ny_src+jlow
    xyidxsrc(3) = (ilow-1)*ny_src+jhigh
    xyidxsrc(4) = (ihigh-1)*ny_src+jhigh
  endif

  ! indices are (y,x)
  ! A11 A12 A21 A22
  !xysrc(1:4) = (/xsrc(ilow), xsrc(ihigh), xsrc(ilow), xsrc(ihigh)/)
  !xysrc(5:8) = (/ysrc(jlow), ysrc(jlow), ysrc(jhigh), ysrc(jhigh)/)
  xysrc(1) = xsrc(ilow)
  xysrc(2) = xsrc(ihigh)
  xysrc(3) = xsrc(ilow)
  xysrc(4) = xsrc(ihigh)
  xysrc(5) = ysrc(jlow)
  xysrc(6) = ysrc(jlow)
  xysrc(7) = ysrc(jhigh)
  xysrc(8) = ysrc(jhigh)

  END SUBROUTINE get_interp_coeffs_lonlat

  SUBROUTINE interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                     & nsrc, &
                     & ny_src, nx_src, nxy_src, &
                     & ny_dst, nx_dst, nxy_dst, &
                     & typ_src, typ_dst, &
                     & redpts_src, redpts_dst, &
                     & lat_row, &
                     & chunk_size_in, abs_tolerance_in)

    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
    use omp_lib
    
    ! Number of fields to interpolate
    integer(intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    
    ! Grids types
    integer(intk), intent(in), value  :: typ_src, typ_dst

    ! Fields source and interpolated
    real(realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(realk), intent(in), dimension(nx_src) :: xsrc
    real(realk), intent(in), dimension(ny_src) :: ysrc
    real(realk), intent(in), dimension(nx_dst) :: xdst
    real(realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_
    integer(intk), optional, intent(in), dimension(ny_src)  :: redpts_src
    integer(intk), optional, intent(in), dimension(ny_dst)  :: redpts_dst

    ! Whether the dimension ordering is (Fortran writing, column-major)
    ! Fortran column-major: (x,y) (true) or (y,x) (false)
    ! C row-major: (y,x) (true) or (x,y) (false)
    ! for the flattened fields on regular grid. For reduced
    ! and unstructured grids, lat_row is ignored
    logical, intent(in), value :: lat_row
    ! Specify chunk size
    integer(intk), intent(in), value, optional :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(realk), intent(in), value, optional :: abs_tolerance_in
    

    ! Default parameters for chunk size and tolerance
    integer(intk) :: actual_chunk_size
    real(realk) :: atol = 1.e-3

    ! variables for OMP
    integer(intk) :: num_threads
    integer(intk) :: ndimomp
    ! loop variables
    integer(intk) :: jdst, idst, idst_min, idst_max, jj
    
    real(realk) :: tgtx, tgty
    integer(intk) :: tgtidx
    real(realk), dimension(8) :: xysrc
    integer(intk), dimension(4) :: xyidxsrc

    real(realk) :: x11,x12,x21,x22
    real(realk) :: y11,y12,y21,y22
    integer(intk) :: xyidxsrc11, xyidxsrc12, xyidxsrc21, xyidxsrc22 

    real(realk) :: timer
    
    ! typ_ = 
    ! 1 = lonlat
    ! 2 = reduced
    ! 3 = unstructured

    
    ! rectangular
    ! present(redpts_) is false
    ! nxy_ = nx*ny
 
    ! unstructured
    ! present(redpts_) is false
    !
    ! nx_ = ny_ = nxy_
    ! case(3) assume ordered by lat,lon
    ! case(4) completely unstructured

    ! reduced gaussian
    ! present(redpts_) is true
    ! nx_ = ny_ != nxy_
    ! x_ indicates the first longitude for each y_
    select case(typ_src)
      case(1)
        if (nxy_src /= (nx_src*ny_src)) then
          call throw_error("nxy = nx*ny needed for rectangular grids!")
        endif
        !write(*,*) "Src grid is rectangular"
      case(2)
        if ((nx_src /= ny_src) .or. (.not. present(redpts_src))) then
          call throw_error("nx_src must be equal to ny_src "&
                       & //"redpts_src must be provided for reduced grids!")
        endif
        !write(*,*) "Src grid is reduced"
      case(3)
      case(4)
        if ((nx_src /= ny_src) .or. (ny_src/= nxy_src)) then
          call throw_error("nx_src = ny_src = nxy_src required for unstructured grids!")
        endif
        !write(*,*) "Src grid is unstructured"
      case default
        call throw_error("Could not recognize typ_src")
    end select

    select case(typ_dst)
      case(1)
        if (nxy_dst /= (nx_dst*ny_dst)) then
          call throw_error("nxy = nx*ny needed for rectangular grids!")
        endif
        !write(*,*) "tgt grid is rectangular"
      case(2)
        if ((nx_dst /= ny_dst) .or. (.not. present(redpts_dst))) then
          call throw_error("nx_dst must have the same length of ny_dst "&
                       & //"redpts_dst must be provided for reduced grids!")
        endif
        !write(*,*) "tgt grid is reduced"
      case(3)
      case(4)
        if ((nx_dst /= ny_dst) .or. (ny_dst/= nxy_dst)) then
          call throw_error("nx_dst = ny_dst = nxy_dst required for unstructured grids!")
        endif
        !write(*,*) "tgt grid is unstructured"
      case default
        call throw_error("Could not recognize typ_dst")
    end select

    timer = real(omp_get_wtime(), realk)

    if (present(abs_tolerance_in)) then
      atol = abs_tolerance_in
    endif


    ! Master thread decides
    num_threads = int(omp_get_max_threads(), intk)
    write(*,'(A,I0,A)') "Interp_2d using ",num_threads," threads"

    ! Parallelize on y
    if (nxy_dst > num_threads) then
        ndimomp = ny_dst
    else
        ndimomp = 1
    endif

    ! Compute chunk size
    if (present(chunk_size_in)) then
      actual_chunk_size = chunk_size_in
    elseif (typ_dst == 1) then
      ! Regular grids: minimum 16 iterations/chunk to reduce scheduling overhead
      actual_chunk_size = max(2, ndimomp / (num_threads * 4))
    elseif (typ_dst == 2) then
      ! This is because the workload for latitude ydst(j) 
      ! is determined by redpts_dst(j)
      actual_chunk_size = 1
    elseif (typ_dst == 3 .or. typ_dst == 4) then
      ! This is because the workload for point idst is determined by the distance to source points
      actual_chunk_size = max(16, ndimomp / (num_threads * 4))
    else
      actual_chunk_size = 1000
    endif

    !$OMP PARALLEL DO SCHEDULE(GUIDED, actual_chunk_size) &
!    !$OMP DEFAULT(NONE) &
    !$OMP PRIVATE(jdst, jj, idst_min, idst_max, idst, &
    !$OMP tgty, tgtx, tgtidx, xysrc, xyidxsrc, &
    !$OMP xyidxsrc11, xyidxsrc12, xyidxsrc21, xyidxsrc22, &
    !$OMP x11, x12, x21, x22, y11, y12, y21, y22)

    ! Along y (lats)        
    do jdst=1,ndimomp
      ! Find x range
      select case(typ_dst)
        case (1)
          idst_min=1
          idst_max=nx_dst
        case(2)
          idst_min=1
          idst_max=redpts_dst(jdst)
        case(3:4)
          idst_min=jdst
          idst_max=jdst
      end select
      tgty = ydst(jdst)
      
      do idst=idst_min,idst_max
        !write(*,*) "Thread ",thread_id," (jdst,idst)=(",jdst,",",idst,")"
        select case(typ_dst)
          !  (y,x) flattening assumed
          case (1)
            if (lat_row) then
              tgtidx = (jdst-1)*nx_dst + idst
            else
              tgtidx = (idst-1)*ny_dst + jdst
            endif
            tgtx = xdst(idst)
          case (2)
            tgtidx = 0
              do jj=1,jdst-1
                tgtidx = tgtidx + redpts_dst(jj)
              end do
            tgtidx = tgtidx + idst
            tgtx = xdst(jdst) + (idst - 1)*360.0/idst_max
          case (3:4)
            tgtidx = idst
            tgtx = xdst(idst)
        end select
        !write(*,*) " Type dst: ",typ_dst
        !write(*,*) " Type src: ",typ_src

        select case(typ_src)
          case (1)
            !write(*,*) "Call interp coeffs lonlat"
            call get_interp_coeffs_lonlat(xsrc, ysrc, nx_src, ny_src, tgtx, tgty, &
                                        & xysrc, xyidxsrc, lat_row, atol_in=atol)
          case (2)
            !write(*,*) "Call interp coeffs reduced"
            call get_interp_coeffs_reduced(ysrc, xsrc, ny_src, redpts_src, &
            & tgtx, tgty, xysrc, xyidxsrc, atol_in=atol)
          !case (3)
          !  call get_interp_coeffs_unstruct(xsrc, ysrc, tgtx, tgty, xysrc(1:8), xyidxsrc(1:4), .true.)
          !case (4)
          !  call get_interp_coeffs_unstruct(xsrc, ysrc, tgtx, tgty, xysrc(1:8), xyidxsrc(1:4), .false.)
        end select
        
        ! indices are (y=lat,x=lon) lat_row
        ! or (y=lon,x=lat) .not. lat_row
        ! A11 A12 A21 A22
        x11 = xysrc(1)
        x12 = xysrc(2)
        x21 = xysrc(3)
        x22 = xysrc(4)

        y11 = xysrc(5)
        y12 = xysrc(6)
        y21 = xysrc(7)
        y22 = xysrc(8)

        xyidxsrc11 = xyidxsrc(1)
        xyidxsrc12 = xyidxsrc(2)
        xyidxsrc21 = xyidxsrc(3)
        xyidxsrc22 = xyidxsrc(4)


        if (lat_row) then !y is lat, x is lon
          fdst(tgtidx,:) = bilinear_lonlat(tgty, tgtx, &
                                        & x11, x12, x21, x22, &
                                        & y11, y12, y21, y22, &
                                        & fsrc(xyidxsrc11,:), fsrc(xyidxsrc12, :), &
                                        & fsrc(xyidxsrc21,:), fsrc(xyidxsrc22, :), &
                                        & nsrc, atol)
        else !x is lat, y is lon
          fdst(tgtidx,:) = bilinear_lonlat(tgtx, tgty, &
                                         & y11, y21, y12, y22, &
                                         & x11, x21, x12, x22, &
                                         & fsrc(xyidxsrc11,:), fsrc(xyidxsrc21, :), &
                                         & fsrc(xyidxsrc12,:), fsrc(xyidxsrc22, :), &                                         
                                         & nsrc, atol)
        !write(*,'(A,I0,A,F6.2,A,F6.2)') "Thread ",thread_id,"tgty=",tgty," tgtx=",tgtx
        endif
        

      enddo        
    enddo
    !$OMP END PARALLEL DO

  !write(*,*)  "OMP time: ",omp_get_wtime()-timer,"s"

  END SUBROUTINE interp_2d

  FUNCTION bilinear_lonlat(tgtlat, tgtlon, &
                         & lon11, lon12, lon21, lon22, &
                         & lat11, lat12, lat21, lat22, &
                         & fsrc11, fsrc12, &
                         & fsrc21, fsrc22, &
                         & nsrc, atol) RESULT(fdst)
    ! Handle bilinear interpolation lon-lat gnostic
    real(realk), intent(in) :: tgtlat, tgtlon
    real(realk), intent(in) :: lon11, lon12, lon21, lon22
    real(realk), intent(in) :: lat11, lat12, lat21, lat22
    real(realk), intent(in), dimension(nsrc) :: fsrc11, fsrc12, fsrc21, fsrc22
    integer(intk), intent(in) :: nsrc
    real(realk), intent(in) :: atol

    real(realk), dimension(nsrc) :: fdst, fxalat, fxblat
    real(realk) :: w, wa, wb, lona, lonb

    ! (lat, lon)
    !                                       . A22(lat22,lon22)
    !      . A21(lat21,lon21)
    !   . Aalat(tgtlat,lona)      . (tgtlat,tgtlon)    . Ablat(tgtlat,lonb)
    ! . A11(lat11,lon11)
    !
    !                                          . A12(lat12,lon12)
        
    if ((get_dist(lon22, lon12)>atol) .and. ((lat22-lat12)>atol)) then ! They are latitude-separated
      w = (lat22-tgtlat)/(lat22-lat12)
      if (lon22-lon12>180) then !The cut line (0 meridian or date line) sits in between
        lonb = w*(lon12+360)+(1-w)*lon22
      else if (lon22-lon12<-180) then !The cut line (0 meridian or date line) sits in between
        lonb = w*lon12+(1-w)*(lon22+360)
      else
        lonb = w*lon12+(1-w)*lon22
      endif
      
      if (lonb > 360) then
        lonb = modulo(lonb,360.)
      endif
    else ! they are the same point
      lonb = lon12
    endif

    if ((get_dist(lonb,lon12)>atol)) then
      wb = sqrt((get_dist(lonb,lon12)**2 + (tgtlat-lat12)**2)/(get_dist(lon22,lon12)**2+(lat22-lat12)**2))
      !wb = get_geodist((tgtlat-lat12), (tgtlat+lat12)/2., get_dist(lonb,lon12))/get_geodist((lat22-lat12), (lat22+lat12)/2., get_dist(lon22,lon12))
    else if (lat22-lat12>atol) then
      wb = (tgtlat-lat12)/(lat22-lat12)
    else
      wb = 0._realk
    endif
    wb = max(0._realk, min(1._realk, wb))

    fxblat(:) = wb*fsrc22(:) + (1-wb)*fsrc12(:)

    ! Find lona and interpolate left side
    if ((get_dist(lon21,lon11)>atol) .and. (lat21-lat11>atol)) then
      w = (lat21-tgtlat)/(lat21-lat11)
      if (lon21-lon11>180) then
        lona = w*(lon11+360)+(1-w)*lon21
      else if (lon21-lon11<-180) then
        lona = w*lon11+(1-w)*(lon21+360)
      else
        lona = w*lon11+(1-w)*lon21
      endif
      if (lona > 360) then
        lona = modulo(lona,360.)
      endif
    else
      lona = lon11
    endif
    if (get_dist(lona, lon11)>atol) then
      wa = sqrt((get_dist(lona,lon11)**2 + (tgtlat-lat11)**2)/(get_dist(lon21,lon11)**2+(lat21-lat11)**2))
      !wa = get_geodist((tgtlat-lat11), (tgtlat+lat11)/2., get_dist(lona,lon11))/get_geodist((lat21-lat11), (lat21+lat11)/2., get_dist(lon21,lon11))
    else if (lat21-lat11>atol) then
      wa = (tgtlat-lat11)/(lat21-lat11)
    else
      wa = 0._realk
    endif
    wa = max(0._realk, min(1._realk, wa))
    fxalat(:) = wa*fsrc21(:) + (1-wa)*fsrc11(:)
    
    if (get_dist(lonb,lona)>atol) then
      w = get_dist(tgtlon,lonb)/get_dist(lonb,lona)
      w = max(0._realk, min(1._realk, w))
      fdst(:) = w*fxalat(:) + (1-w)*fxblat(:)
    else
      fdst(:) = fxalat(:)
    endif

  END FUNCTION bilinear_lonlat


  SUBROUTINE f_interp_2d_rec2rec(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & lat_row, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_rec2rec")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_

    ! row/col source flattening order switch
    integer(c_int), intent(in), value :: lat_row
    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=1
    integer(c_intk) :: typ_dst=1

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                 & nsrc, &
                 & ny_src, nx_src, nxy_src, &
                 & ny_dst, nx_dst, nxy_dst, &
                 & typ_src, typ_dst, &
                 & lat_row=(lat_row /= 0), &
                 & chunk_size_in=chunk_size_in, abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_rec2rec

  SUBROUTINE f_interp_2d_rec2red(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & redpts_dst, &
                               & lat_row, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_rec2red")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_
    integer(c_intk), intent(in), dimension(ny_dst)  :: redpts_dst

    ! row/col source flattening order switch
    integer(c_int), intent(in), value :: lat_row
    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=1
    integer(c_intk) :: typ_dst=2

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                 & nsrc, &
                 & ny_src, nx_src, nxy_src, &
                 & ny_dst, nx_dst, nxy_dst, &
                 & typ_src, typ_dst, &
                 & redpts_dst=redpts_dst, &
                 & lat_row=(lat_row /= 0),  &
                 & chunk_size_in=chunk_size_in, abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_rec2red
  
  SUBROUTINE f_interp_2d_red2rec(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & redpts_src, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_red2rec")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_
    integer(c_intk), intent(in), dimension(ny_src)  :: redpts_src

    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=2
    integer(c_intk) :: typ_dst=1

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                 & nsrc, &
                 & ny_src, nx_src, nxy_src, &
                 & ny_dst, nx_dst, nxy_dst, &
                 & typ_src, typ_dst, &
                 & redpts_src=redpts_src, &
                 & lat_row=.true.,  &
                 & chunk_size_in=chunk_size_in, abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_red2rec

  SUBROUTINE f_interp_2d_red2red(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & redpts_src, redpts_dst, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_red2red")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_
    integer(c_intk), intent(in), dimension(ny_src)  :: redpts_src
    integer(c_intk), intent(in), dimension(ny_dst)  :: redpts_dst

    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=2
    integer(c_intk) :: typ_dst=2

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                  & nsrc, &
                  & ny_src, nx_src, nxy_src, &
                  & ny_dst, nx_dst, nxy_dst, &
                  & typ_src, typ_dst, &
                  & redpts_src=redpts_src, redpts_dst=redpts_dst, &
                  & lat_row=.true.,  chunk_size_in=chunk_size_in, abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_red2red

  SUBROUTINE f_interp_2d_red2uns(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & redpts_src, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_red2uns")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst
    
    ! Only for reduced grids - x_ indicates the first longitude for each y_
    ! and redpts_ defines how many longitude points between x_ and 360. 
    ! per each y_. x_ <= 360./redpts_
    integer(c_intk), intent(in), dimension(ny_src)  :: redpts_src

    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=2
    integer(c_intk) :: typ_dst=3

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                  & nsrc, &
                  & ny_src, nx_src, nxy_src, &
                  & ny_dst, nx_dst, nxy_dst, &
                  & typ_src, typ_dst, &
                  & redpts_src=redpts_src, &
                  & lat_row=.true.,  chunk_size_in=chunk_size_in, abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_red2uns

  SUBROUTINE f_interp_2d_rec2uns(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                               & nsrc, &
                               & ny_src, nx_src, nxy_src, &
                               & ny_dst, nx_dst, nxy_dst, &
                               & lat_row, &
                               & chunk_size_in, abs_tolerance_in) &
            & bind(C, name="interp_2d_rec2uns")
    
    ! Number of fields to interpolate
    integer(c_intk), intent(in), value  :: nsrc

    ! Grids shapes
    integer(c_intk), intent(in), value  :: ny_src, nx_src, nxy_src
    integer(c_intk), intent(in), value  :: ny_dst, nx_dst, nxy_dst
    

    ! Fields source and interpolated
    real(c_realk), intent(in), dimension(nxy_src,nsrc) :: fsrc
    real(c_realk), intent(out), dimension(nxy_dst,nsrc) :: fdst
    
    ! Grid points
    real(c_realk), intent(in), dimension(nx_src) :: xsrc
    real(c_realk), intent(in), dimension(ny_src) :: ysrc
    real(c_realk), intent(in), dimension(nx_dst) :: xdst
    real(c_realk), intent(in), dimension(ny_dst) :: ydst

    ! row/col flattening order switch
    integer(c_int), intent(in), value :: lat_row 
    ! Specify chunk size
    integer(c_intk), intent(in), value :: chunk_size_in
    ! Specify absolute tolerance (default is 10^-3 degrees)
    real(c_realk), intent(in), value :: abs_tolerance_in
    
    ! Grids types
    integer(c_intk) :: typ_src=1
    integer(c_intk) :: typ_dst=3

    call interp_2d(fsrc, fdst, ysrc, xsrc, ydst, xdst, &
                 & nsrc, &
                 & ny_src, nx_src, nxy_src, &
                 & ny_dst, nx_dst, nxy_dst, &
                 & typ_src, typ_dst, &
                 & lat_row=(lat_row /= 0), &
                 & chunk_size_in=chunk_size_in, &
                 & abs_tolerance_in=abs_tolerance_in)
  
  END SUBROUTINE f_interp_2d_rec2uns
  
END MODULE mo_interp2d
