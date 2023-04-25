module tetrahedral_interpolation

use grid

implicit none

type point
  real :: x, y, z
end type point

type polygon
  type(point), dimension(:), allocatable :: vertex
end type polygon

CONTAINS

!-------------------------------------

subroutine calculate_surface(num_vals, values, mask, &
                             u_area_frac, v_area_frac, w_area_frac)
  integer :: num_vals
  real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, num_vals) :: values
  logical mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real u_area_frac(dimx1_u:dimx2_u, dimy1_u:dimy2_u, nzm)
  real v_area_frac(dimx1_v:dimx2_v, dimy1_v:dimy2_v, nzm)
  real w_area_frac(dimx1_w:dimx2_w, dimy1_w:dimy2_w, nz )

  real, dimension(3, 3, 3, num_vals) :: temp_values
  logical temp_mask(3, 3, 3)
  integer i, j, k

  u_area_frac(:, :, 1) = 0.
  u_area_frac(:, :, nz) = 0.
  v_area_frac(:, :, 1) = 0.
  v_area_frac(:, :, nz) = 0.
  w_area_frac(:, :, 1) = 0.
  w_area_frac(:, :, nz) = 0.

  do k=2,nzm-1
    do j=1,ny+1
      do i=1,nx+1
        if ( any( all(values(i-1:i+1, j-1:j+1, k-1:k+1, :) > 0., dim=4)) ) then
          temp_values = average_cgrid( values(i-1:i+1, j-1:j+1, k-1:k+1, :) )
          temp_mask = all(temp_values > 0., dim=4)
          u_area_frac(i, j, k) = calc_area( temp_values(1, :, :, :), &
                                            temp_mask(1, :, :) )
          v_area_frac(i, j, k) = calc_area( temp_values(:, 1, :, :), &
                                            temp_mask(:, 1, :) )
          w_area_frac(i, j, k) = calc_area( temp_values(:, :, 1, :), &
                                            temp_mask(:, :, 1) )
        else
          u_area_frac(i, j, k) = 0.
          v_area_frac(i, j, k) = 0.
          w_area_frac(i, j, k) = 0.
        end if
      end do
    end do
  end do

end subroutine calculate_surface

!-------------------------------------

function calc_area( values, mask )
  real calc_area
  real, dimension(:, :, :) :: values
  logical, dimension(3, 3) :: mask

  real, dimension(9) :: x_vals = (/ 0., 0., 0., .5, 1., 1., 1., .5, 0./)
  real, dimension(9) :: y_vals = (/ 0., .5, 1., 1., 1., .5, 0., 0., 0./)
  integer, dimension(9) :: i_vals = (/ 1, 1, 1, 2, 3, 3, 3, 2, 1/)
  integer, dimension(9) :: j_vals = (/ 1, 2, 3, 3, 3, 2, 1, 1, 1/)

  real, dimension(18) :: x, y
  real x1, x2, frac
  real delta1, delta2

  integer n, m

  m = 1
  do n=1,8
    frac = calculate_position(mask(2, 2), &
                              values(2, 2, :), &
                              mask(i_vals(n), j_vals(n)), &
                              values(i_vals(n), j_vals(n), :))

    x(m) = x_vals(n)*frac + .5*(1.-frac)
    y(m) = y_vals(n)*frac + .5*(1.-frac)
    m = m + 1
    
    if (mask(i_vals(n), j_vals(n)) .neqv. mask(i_vals(n+1), j_vals(n+1))) then
      frac = calculate_position(mask(i_vals(n), j_vals(n)), &
                                values(i_vals(n), j_vals(n), :), &
                                mask(i_vals(n+1), j_vals(n+1)), &
                                values(i_vals(n+1), j_vals(n+1), :))

      x(m) = x_vals(n+1)*frac + x_vals(n)*(1.-frac)
      y(m) = y_vals(n+1)*frac + y_vals(n)*(1.-frac)
      m = m + 1
    end if
    
  end do

  frac = calculate_position(mask(2, 2), &
                            values(2, 2, :), &
                            mask(i_vals(9), j_vals(9)), &
                            values(i_vals(9), j_vals(9), :))

  x(m) = x_vals(9)*frac + .5*(1.-frac)
  y(m) = y_vals(9)*frac + .5*(1.-frac)

  x(m+1) = -2.
  y(m+1) = -2.

  calc_area = 0.
  n = 1
  do while (x(n) > -1.)
    calc_area = calc_area - x(n)*y(n+1)
    calc_area = calc_area + y(n)*x(n+1)
    n = n + 1
  end do

  calc_area = 0.5*calc_area

  if (.not. mask(2, 2)) then
    calc_area = 1. - calc_area
  end if

end function calc_area

!---------------------------------

function average_cgrid( values )
  real, dimension(:, :, :, :) :: values
  real, dimension(3, 3, 3, size(values, dim=4)) :: average_cgrid
  integer n

  do n=1,size(values, dim=4)
    average_cgrid(1, 1, 1, n) = sum(values(1:2, 1:2, 1:2, n))/8.
    average_cgrid(1, 1, 2, n) = sum(values(1:2, 1:2, 2, n))/4.
    average_cgrid(1, 1, 3, n) = sum(values(1:2, 1:2, 2:3, n))/8.
    average_cgrid(1, 2, 1, n) = sum(values(1:2, 2, 1:2, n))/4.
    average_cgrid(1, 2, 2, n) = sum(values(1:2, 2, 2, n))/2.
    average_cgrid(1, 2, 3, n) = sum(values(1:2, 2, 2:3, n))/4.
    average_cgrid(1, 3, 1, n) = sum(values(1:2, 2:3, 1:2, n))/8.
    average_cgrid(1, 3, 2, n) = sum(values(1:2, 2:3, 2, n))/4.
    average_cgrid(1, 3, 3, n) = sum(values(1:2, 2:3, 2:3, n))/8.
    average_cgrid(2, 1, 1, n) = sum(values(2, 1:2, 1:2, n))/4.
    average_cgrid(2, 1, 2, n) = sum(values(2, 1:2, 2, n))/2.
    average_cgrid(2, 1, 3, n) = sum(values(2, 1:2, 2:3, n))/4.
    average_cgrid(2, 2, 1, n) = sum(values(2, 2, 1:2, n))/2.
    average_cgrid(2, 2, 2, n) = values(2, 2, 2, n)
    average_cgrid(2, 2, 3, n) = sum(values(2, 2, 2:3, n))/2.
    average_cgrid(2, 3, 1, n) = sum(values(2, 2:3, 1:2, n))/4.
    average_cgrid(2, 3, 2, n) = sum(values(2, 2:3, 2, n))/2.
    average_cgrid(2, 3, 3, n) = sum(values(2, 2:3, 2:3, n))/4.
    average_cgrid(3, 1, 1, n) = sum(values(2:3, 1:2, 1:2, n))/8.
    average_cgrid(3, 1, 2, n) = sum(values(2:3, 1:2, 2, n))/4.
    average_cgrid(3, 1, 3, n) = sum(values(2:3, 1:2, 2:3, n))/8.
    average_cgrid(3, 2, 1, n) = sum(values(2:3, 2, 1:2, n))/4.
    average_cgrid(3, 2, 2, n) = sum(values(2:3, 2, 2, n))/2.
    average_cgrid(3, 2, 3, n) = sum(values(2:3, 2, 2:3, n))/4.
    average_cgrid(3, 3, 1, n) = sum(values(2:3, 2:3, 1:2, n))/8.
    average_cgrid(3, 3, 2, n) = sum(values(2:3, 2:3, 2, n))/4.
    average_cgrid(3, 3, 3, n) = sum(values(2:3, 2:3, 2:3, n))/8.
  end do
  
end function average_cgrid

!---------------------------------

real function calculate_position(mask1, values1, &
                                 mask2, values2)
  logical mask1, mask2
  real, dimension(:) :: values1
  real, dimension(:) :: values2

  real x1, x2
  
  integer n
  n = 1
  
  if (mask1 .neqv. mask2) then
    if (mask1) then
      ! all values1 > 0
      x1 = 1.
      do while (n .le. size(values1))
        if (values2(n) > 0.) then
          x2 = 1.
        else
          x2 = values1(n) / (values1(n) - values2(n))
        end if
        x1 = min(x1, x2)
        n = n + 1
      end do
      
    else
      ! all values2 > 0    
      x1 = 0.
      do while (n .le. size(values1))
        if (values1(n) > 0.) then
          x2 = 0.
        else
          x2 = values1(n) / (values1(n) - values2(n))
        end if
        x1 = max(x1, x2)
        n = n + 1
      end do
    end if
    
    calculate_position = x1
    
  else
    calculate_position = 1.
  end if

end function calculate_position

!===============================

subroutine calculate_volume(num_vals, values, mask, volume)
  ! Calculates volume of cloud in each grid cell
  ! by dividing each cell into 6 pyramids
  ! All volumes are a nondimensional fraction of the cell volume.
  integer num_vals
  real, dimension(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm, num_vals) :: values
  logical mask(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)
  real volume(dimx1_s:dimx2_s, dimy1_s:dimy2_s, nzm)

  real, dimension(3, 3, 3, num_vals) :: temp_values
  logical temp_mask(3, 3, 3)

  integer i, j, k, n

  do k=2,nzm-1
    do j=1,ny
      do i=1,nx
        if (all(mask(i-1:i+1, j-1:j+1, k-1:k+1))) then
          volume(i, j, k) = 1.
          cycle
        else if (.not. any(mask(i-1:i+1, j-1:j+1, k-1:k+1))) then
          volume(i, j, k) = 0.
          cycle
        else        
          temp_values = average_cgrid( values(i-1:i+1, j-1:j+1, k-1:k+1, :) )
          temp_mask = all(temp_values > 0., dim=4)

          if (all(temp_mask)) then
            volume(i,j,k) = 1.
          else if (.not. any(temp_mask)) then
            volume(i,j,k) = 0.
          else
            call calc_tet_vol(temp_values, temp_mask, volume(i, j, k))
          end if
        end if
      end do
    end do
  end do
  
end subroutine calculate_volume

!---------------------------------

subroutine calc_tet_vol(values, mask, volume)
  real, dimension(:, :, :, :) :: values
  logical, dimension(3, 3, 3) :: mask
  real :: volume

  integer n

  real, dimension(size(values, dim=4)) :: core_values, face_values, edge_values, corner_values
  logical :: core_mask, face_mask, edge_mask, corner_mask

  integer, dimension(48) :: i_face, j_face, k_face
  integer, dimension(48) :: i_edge, j_edge, k_edge
  integer, dimension(48) :: i_corner, j_corner, k_corner

  i_face = (/1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2 /)
  j_face = (/2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2 /)
  k_face = (/2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             2, 2, 2, 2, 2, 2, 2, 2, &
             1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3 /)

  i_edge = (/1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3, &
             1, 1, 3, 3, 2, 2, 2, 2, &
             1, 1, 3, 3, 2, 2, 2, 2, &
             1, 1, 3, 3, 2, 2, 2, 2, &
             1, 1, 3, 3, 2, 2, 2, 2 /)
  j_edge = (/1, 1, 3, 3, 2, 2, 2, 2, &
             1, 1, 3, 3, 2, 2, 2, 2, &
             1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3, &
             2, 2, 2, 2, 1, 1, 3, 3, &
             2, 2, 2, 2, 1, 1, 3, 3 /)
  k_edge = (/2, 2, 2, 2, 1, 1, 3, 3, &
             2, 2, 2, 2, 1, 1, 3, 3, &
             2, 2, 2, 2, 1, 1, 3, 3, &
             2, 2, 2, 2, 1, 1, 3, 3, &
             1, 1, 1, 1, 1, 1, 1, 1, &
             3, 3, 3, 3, 3, 3, 3, 3 /)

  i_corner = (/1, 1, 1, 1, 1, 1, 1, 1, &
               3, 3, 3, 3, 3, 3, 3, 3, &
               1, 1, 3, 3, 1, 3, 1, 3, &
               1, 1, 3, 3, 1, 3, 1, 3, &
               1, 1, 3, 3, 1, 3, 1, 3, &
               1, 1, 3, 3, 1, 3, 1, 3 /)
  j_corner = (/1, 1, 3, 3, 1, 3, 1, 3, &
               1, 1, 3, 3, 1, 3, 1, 3, &
               1, 1, 1, 1, 1, 1, 1, 1, &
               3, 3, 3, 3, 3, 3, 3, 3, &
               1, 3, 1, 3, 1, 1, 3, 3, &
               1, 3, 1, 3, 1, 1, 3, 3 /)
  k_corner = (/1, 3, 1, 3, 1, 1, 3, 3, &
               1, 3, 1, 3, 1, 1, 3, 3, &
               1, 3, 1, 3, 1, 1, 3, 3, &
               1, 3, 1, 3, 1, 1, 3, 3, &
               1, 1, 1, 1, 1, 1, 1, 1, &
               3, 3, 3, 3, 3, 3, 3, 3 /)

  volume = 0.

  core_values = values(2, 2, 2, :)
  core_mask = mask(2, 2, 2)

  do n=1,48
    face_values = values(i_face(n), j_face(n), k_face(n), :)
    face_mask = mask(i_face(n), j_face(n), k_face(n))

    edge_values = values(i_edge(n), j_edge(n), k_edge(n), :)
    edge_mask = mask(i_edge(n), j_edge(n),  k_edge(n))

    corner_values = values(i_corner(n), j_corner(n), k_corner(n), :)
    corner_mask = mask(i_corner(n), j_corner(n), k_corner(n))
        
    volume = volume + calc_vol(core_values, core_mask, &
                               face_values, face_mask, &
                               edge_values, edge_mask, &
                               corner_values, corner_mask)
  end do

  if (.not. mask(2, 2, 2)) then
     volume = 1. - volume
  end if

end subroutine calc_tet_vol

!-------------------------------------

real function calc_vol(core_values, core_mask, &
                       face_values, face_mask, &
                       edge_values, edge_mask, &
                       corner_values, corner_mask)
  real, dimension(:) :: core_values, face_values, edge_values, corner_values
  logical :: core_mask, face_mask, edge_mask, corner_mask

  integer :: tetra_type, n
  type(polygon), dimension(:), allocatable :: poly_array

  type(point) :: core, face, edge, corner
  type(point) :: core_face, core_edge, core_corner
  type(point) :: face_edge, face_corner
  type(point) :: edge_corner
    
  core = point(0., 0., 0.)
  face = point(0., 0., .5)
  edge = point(0., .5, .5)
  corner = point(.5, .5, .5)

  calc_vol = 0.
  tetra_type = 0
  if (core_mask .eqv. face_mask) then
    tetra_type = tetra_type + 1
  end if
  if (core_mask .eqv. edge_mask) then 
    tetra_type = tetra_type + 2
  end if
  if (core_mask .eqv. corner_mask) then 
    tetra_type = tetra_type + 4
  end if

  ! Note: any polygon with the origin in its plane doesn't add to the volume calculation, 
  ! so we ignore those surfaces in the following calculations.
  select case(tetra_type)
    case (0)
      ! No points same as Core
      core_face = interpolate_point(core, core_values, core_mask, &
                                    face, face_values, face_mask)
      core_edge = interpolate_point(core, core_values, core_mask, &
                                    edge, edge_values, edge_mask)
      core_corner = interpolate_point(core, core_values, core_mask, &
                                      corner, corner_values, corner_mask)
      
      allocate(poly_array(1))

      allocate(poly_array(1)%vertex(3))
      poly_array(1)%vertex(1) = core_face
      poly_array(1)%vertex(2) = core_corner
      poly_array(1)%vertex(3) = core_edge

    case (1)
      ! Face point same as Core
      core_edge = interpolate_point(core, core_values, core_mask, &
                                    edge, edge_values, edge_mask)
      core_corner = interpolate_point(core, core_values, core_mask, &
                                      corner, corner_values, corner_mask)
      face_edge = interpolate_point(face, face_values, face_mask, &
                                    edge, edge_values, edge_mask)
      face_corner = interpolate_point(face, face_values, face_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(3))

      allocate(poly_array(1)%vertex(3))
      poly_array(1)%vertex(1) = face
      poly_array(1)%vertex(2) = face_corner
      poly_array(1)%vertex(3) = face_edge

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_edge
      poly_array(2)%vertex(2) = face_edge
      poly_array(2)%vertex(3) = core_corner

      allocate(poly_array(3)%vertex(3))
      poly_array(3)%vertex(1) = face_edge
      poly_array(3)%vertex(2) = face_corner
      poly_array(3)%vertex(3) = core_corner

    case (2)
      ! Edge point same as Core
      core_face = interpolate_point(core, core_values, core_mask, &
                                    face, face_values, face_mask)
      core_corner = interpolate_point(core, core_values, core_mask, &
                                      corner, corner_values, corner_mask)
      face_edge = interpolate_point(face, face_values, face_mask, &
                                    edge, edge_values, edge_mask)
      edge_corner = interpolate_point(edge, edge_values, edge_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(3))

      allocate(poly_array(1)%vertex(3))
      poly_array(1)%vertex(1) = edge
      poly_array(1)%vertex(2) = face_edge
      poly_array(1)%vertex(3) = edge_corner

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_face
      poly_array(2)%vertex(2) = core_corner
      poly_array(2)%vertex(3) = face_edge

      allocate(poly_array(3)%vertex(3))
      poly_array(3)%vertex(1) = face_edge
      poly_array(3)%vertex(2) = core_corner
      poly_array(3)%vertex(3) = edge_corner

    case (3)
      ! Face and Edge points same as Core
      core_corner = interpolate_point(core, core_values, core_mask, &
                                      corner, corner_values, corner_mask)
      face_corner = interpolate_point(face, face_values, face_mask, &
                                      corner, corner_values, corner_mask)
      edge_corner = interpolate_point(edge, edge_values, edge_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(2))

      allocate(poly_array(1)%vertex(4))
      poly_array(1)%vertex(1) = face
      poly_array(1)%vertex(2) = face_corner
      poly_array(1)%vertex(3) = edge_corner
      poly_array(1)%vertex(4) = edge

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_corner
      poly_array(2)%vertex(2) = edge_corner
      poly_array(2)%vertex(3) = face_corner

    case (4)
      ! Corner point same as Core
      core_face = interpolate_point(core, core_values, core_mask, &
                                    face, face_values, face_mask)
      core_edge = interpolate_point(core, core_values, core_mask, &
                                    edge, edge_values, edge_mask)
      edge_corner = interpolate_point(edge, edge_values, edge_mask, & 
                                      corner, corner_values, corner_mask)
      face_corner = interpolate_point(face, face_values, face_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(3))

      allocate(poly_array(1)%vertex(3))
      poly_array(1)%vertex(1) = corner
      poly_array(1)%vertex(2) = edge_corner
      poly_array(1)%vertex(3) = face_corner

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_face
      poly_array(2)%vertex(2) = face_corner
      poly_array(2)%vertex(3) = edge_corner

      allocate(poly_array(3)%vertex(3))
      poly_array(3)%vertex(1) = core_edge
      poly_array(3)%vertex(2) = core_face
      poly_array(3)%vertex(3) = edge_corner

    case (5)
      ! Corner and Face points same as Core
      core_edge = interpolate_point(core, core_values, core_mask, &
                                    edge, edge_values, edge_mask)
      face_edge = interpolate_point(face, face_values, face_mask, &
                                    edge, edge_values, edge_mask)
      edge_corner = interpolate_point(edge, edge_values, edge_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(2))

      allocate(poly_array(1)%vertex(4))
      poly_array(1)%vertex(1) = face
      poly_array(1)%vertex(2) = corner
      poly_array(1)%vertex(3) = edge_corner
      poly_array(1)%vertex(4) = face_edge

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_edge
      poly_array(2)%vertex(2) = face_edge
      poly_array(2)%vertex(3) = edge_corner

    case (6)
      ! Corner and Edge points same as Core
      core_face = interpolate_point(core, core_values, core_mask, &
                                    face, face_values, face_mask)
      face_edge = interpolate_point(face, face_values, face_mask, &
                                    edge, edge_values, edge_mask)
      face_corner = interpolate_point(face, face_values, face_mask, &
                                      corner, corner_values, corner_mask)

      allocate(poly_array(2))

      allocate(poly_array(1)%vertex(4))
      poly_array(1)%vertex(1) = corner
      poly_array(1)%vertex(2) = edge
      poly_array(1)%vertex(3) = face_edge
      poly_array(1)%vertex(4) = face_corner

      allocate(poly_array(2)%vertex(3))
      poly_array(2)%vertex(1) = core_face
      poly_array(2)%vertex(2) = face_corner
      poly_array(2)%vertex(3) = face_edge

    case (7)
      ! Corner and Edge and Face points same as Core        
      calc_vol = 1./48.
      return
  end select

  do n=1,size(poly_array)
    calc_vol = calc_vol + polygon_volume(poly_array(n))

    deallocate(poly_array(n)%vertex)
  end do
  
  deallocate(poly_array)

end function calc_vol

!-------------------------------

function interpolate_point(point1, values1, mask1, &
                           point2, values2, mask2)
  type(point) interpolate_point
  type(point) point1, point2
  real values1(:)
  real values2(:)
  logical mask1, mask2

  real frac

  frac = calculate_position(mask1, values1, mask2, values2)

  interpolate_point = point( point2%x*frac + point1%x*(1.-frac), &
                             point2%y*frac + point1%y*(1.-frac), &
                             point2%z*frac + point1%z*(1.-frac) )

end function interpolate_point

!-------------------------------

real function polygon_volume(poly)
  type(polygon) :: poly

  real ax, ay, az
  real bx, by, bz
  real cx, cy, cz
  integer n

  polygon_volume = 0.
  
  ax = poly%vertex(1)%x
  ay = poly%vertex(1)%y
  az = poly%vertex(1)%z
  do n=2,(size(poly%vertex)-1)
    bx = poly%vertex(n)%x
    by = poly%vertex(n)%y
    bz = poly%vertex(n)%z

    cx = poly%vertex(n+1)%x
    cy = poly%vertex(n+1)%y
    cz = poly%vertex(n+1)%z

    polygon_volume = polygon_volume + ax*(by*cz - bz*cy)
    polygon_volume = polygon_volume + ay*(bz*cx - bx*cz)
    polygon_volume = polygon_volume + az*(bx*cy - by*cx)
  end do

  polygon_volume = polygon_volume/6.

end function polygon_volume

end module tetrahedral_interpolation
