module soputils
implicit none

contains

subroutine rho(carray, ndat, g_ref, e_ref, ndim, gdim, tdim, rms)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: ndim
    integer, intent(in) :: ndat
    integer, dimension(ndim), intent(in) :: gdim
    integer, dimension(ndim), intent(in) :: tdim
    real(dp), dimension(ndat), intent(in) :: e_ref
    real(dp), dimension(ndat, ndim), intent(in) ::g_ref 
    real(dp), dimension(:), intent(in) :: carray
    real(dp), intent(out) :: rms

    integer :: i, idx, j, k, jkappa, tkappa, unt
    integer :: cheblength, nchebs
    real(dp), dimension(:), allocatable :: u_vects
    real(dp), dimension(:), allocatable :: core
    real(dp), dimension(:), allocatable :: coef_u_vects
    real(dp), dimension(ndat) :: e_sop
    real(dp) :: serieval

    ! Allocate parameter arrays

    nchebs = sum(gdim * tdim)
    allocate(coef_u_vects(nchebs))
    allocate(core(product(gdim)))

    coef_u_vects = carray(:nchebs)
    core = carray(nchebs + 1:)

    ! Allocate Chebyshev vectors

    cheblength = sum(gdim)
    allocate(u_vects(cheblength))

    ! Main loop to compute SOP energy
    
    e_sop = 0.0_dp

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(core, ndim, gdim, g_ref, tdim, coef_u_vects, e_sop)
    do idx = 1, ndat
      
    u_vects = 0.0_dp

    ! Initialize factor vectors

    jkappa = 0
    tkappa = 0
    do i = 1, ndim
       do j = 1, gdim(i)
          jkappa = jkappa + 1
          serieval = 0.0_dp
          do k = 1, tdim(i)
             tkappa = tkappa + 1
             serieval = serieval + coef_u_vects(tkappa) * chebpoly(g_ref(idx, i), k - 1)
          enddo
          u_vects(jkappa) = u_vects(jkappa) + serieval
       enddo
    enddo

    e_sop(idx) = tucker_vect(ndim, gdim, u_vects, core)
 
    enddo
    !$OMP END PARALLEL DO
 
    ! Compute RMS

    rms = sqrt(sum((e_sop - e_ref) ** 2) / ndat)
 
    ! Deallocate resources
 
    deallocate(u_vects)
    deallocate(core)
    deallocate(coef_u_vects)

    ! Write RMSE to file

    open(newunit=unt, file="rmse", position="append", status="unknown")
    write(unt, *) rms
    close(unt)

end subroutine rho

subroutine rho_cheb(coef_u_vects, core, ndat, g_ref, e_ref, ndim, gdim, tdim, rms)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: ndim
    integer, intent(in) :: ndat
    integer, dimension(ndim), intent(in) :: gdim
    integer, dimension(ndim), intent(in) :: tdim
    real(dp), dimension(ndat), intent(in) :: e_ref
    real(dp), dimension(ndat, ndim), intent(in) ::g_ref 
    real(dp), dimension(:), intent(in) :: coef_u_vects
    real(dp), dimension(:), intent(in) :: core
    real(dp), intent(out) :: rms

    integer :: i, idx, j, k, jkappa, tkappa, unt
    integer :: cheblength
    real(dp), dimension(:), allocatable :: u_vects
    real(dp), dimension(ndat) :: e_sop
    real(dp) :: serieval

    ! Allocate Chebyshev vectors

    cheblength = sum(gdim)
    allocate(u_vects(cheblength))

    ! Main loop to compute SOP energy
    
    e_sop = 0.0_dp

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(core, ndim, gdim, g_ref, tdim, coef_u_vects, e_sop)
    do idx = 1, ndat
      
    u_vects = 0.0_dp

    ! Initialize factor vectors

    jkappa = 0
    tkappa = 0
    do i = 1, ndim
       do j = 1, gdim(i)
          jkappa = jkappa + 1
          serieval = 0.0_dp
          do k = 1, tdim(i)
             tkappa = tkappa + 1
             serieval = serieval + coef_u_vects(tkappa) * chebpoly(g_ref(idx, i), k - 1)
          enddo
          u_vects(jkappa) = u_vects(jkappa) + serieval
       enddo
    enddo

    e_sop(idx) = tucker_vect(ndim, gdim, u_vects, core)
 
    enddo
    !$OMP END PARALLEL DO
 
    ! Compute RMS

    rms = sqrt(sum((e_sop - e_ref) ** 2) / ndat)
 
    ! Deallocate resources
 
    deallocate(u_vects)

    ! Write RMSE to file

    open(newunit=unt, file="rmse", position="append", status="unknown")
    write(unt, *) rms
    close(unt)

end subroutine rho_cheb

subroutine rho_core(core, coef_u_vects, ndat, g_ref, e_ref, ndim, gdim, tdim, rms)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: ndim
    integer, intent(in) :: ndat
    integer, dimension(ndim), intent(in) :: gdim
    integer, dimension(ndim), intent(in) :: tdim
    real(dp), dimension(ndat), intent(in) :: e_ref
    real(dp), dimension(ndat, ndim), intent(in) ::g_ref 
    real(dp), dimension(:), intent(in) :: coef_u_vects
    real(dp), dimension(:), intent(in) :: core
    real(dp), intent(out) :: rms

    integer :: i, idx, j, k, jkappa, tkappa, unt
    integer :: cheblength
    real(dp), dimension(:), allocatable :: u_vects
    real(dp), dimension(ndat) :: e_sop
    real(dp) :: serieval

    ! Allocate Chebyshev vectors

    cheblength = sum(gdim)
    allocate(u_vects(cheblength))

    ! Main loop to compute SOP energy
    
    e_sop = 0.0_dp

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(core, ndim, gdim, g_ref, tdim, coef_u_vects, e_sop)
    do idx = 1, ndat
      
    u_vects = 0.0_dp

    ! Initialize factor vectors

    jkappa = 0
    tkappa = 0
    do i = 1, ndim
       do j = 1, gdim(i)
          jkappa = jkappa + 1
          serieval = 0.0_dp
          do k = 1, tdim(i)
             tkappa = tkappa + 1
             serieval = serieval + coef_u_vects(tkappa) * chebpoly(g_ref(idx, i), k - 1)
          enddo
          u_vects(jkappa) = u_vects(jkappa) + serieval
       enddo
    enddo

    e_sop(idx) = tucker_vect(ndim, gdim, u_vects, core)
 
    enddo
    !$OMP END PARALLEL DO
 
    ! Compute RMS

    rms = sqrt(sum((e_sop - e_ref) ** 2) / ndat)
 
    ! Deallocate resources
 
    deallocate(u_vects)

    ! Write RMSE to file

    open(newunit=unt, file="rmse", position="append", status="unknown")
    write(unt, *) rms
    close(unt)

end subroutine rho_core


function tucker_vect(ndim, gdim, u_vects, core) result(prodoned)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(in) :: gdim
    real(dp), dimension(:), intent(in) :: u_vects
    real(dp), dimension(:), intent(in) :: core
    real(dp) :: prodoned

    integer :: mode, newsize, chebslice
    real(dp), dimension(:), allocatable :: tensor_holder
    real(dp), dimension(:), allocatable :: tensor_prod

    ! Compute tensor n-mode product with all vectors

    if (allocated(tensor_holder)) deallocate(tensor_holder)
    allocate(tensor_holder(product(gdim)))

    tensor_holder = core
    chebslice = 1
 
    do mode = 1, ndim
 
       if (allocated(tensor_prod)) deallocate(tensor_prod)
 
       newsize = product(gdim(mode:)) / gdim(mode)
       allocate(tensor_prod(newsize))
 
       tensor_prod = n_mode(mode, gdim, newsize, &
                     tensor_holder, u_vects(chebslice:(chebslice + gdim(mode) - 1)))
 
       deallocate(tensor_holder)
       allocate(tensor_holder(newsize))
       tensor_holder = tensor_prod
 
       chebslice = chebslice + gdim(mode)
       
    enddo

    prodoned = tensor_prod(1)

    deallocate(tensor_holder)
    deallocate(tensor_prod)
    
end function tucker_vect

recursive function chebpoly(x_point, degree) result(polval)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    real(dp), intent(in) :: x_point
    integer, intent(in) :: degree
    real(dp) :: polval

    if (degree == 0) then
        polval = 1.0_dp
    else if (degree == 1) then
        polval = x_point
    else
        polval = 2.0_dp * x_point * chebpoly(x_point, degree - 1) - chebpoly(x_point, degree - 2)
    endif
end function chebpoly

function n_mode(mode, gdim, newsize, core, spp) result(tenprod)
    implicit none
    integer, parameter :: dp = kind(1.d0)
    integer, intent(in) :: mode
    integer, intent(in) :: newsize
    integer, dimension(:), intent(in) :: gdim
    real(dp), dimension(:), intent(in) :: spp 
    real(dp), dimension(:), intent(in) :: core 
    real(dp), dimension(newsize) :: tenprod

    integer, dimension(2) :: newshape

    if (size(gdim(mode:)) == 1) then
        tenprod = dot_product(core, spp)
    else
        newshape = (/newsize, gdim(mode)/)
        tenprod = matmul(reshape(core, newshape), spp)
    end if

end function n_mode

end module soputils
