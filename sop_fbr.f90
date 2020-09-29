module soputils
implicit none

contains

subroutine rho(carray, ndat, g_ref, e_ref, ndim, gdim, tdim, rms)
    implicit none
    integer, intent(in) :: ndim
    integer, intent(in) :: ndat
    integer, intent(in) :: tdim
    integer, dimension(ndim), intent(in) :: gdim
    real, dimension(ndat), intent(in) :: e_ref
    real, dimension(ndat, ndim), intent(in) ::g_ref 
    real, dimension(:), intent(in) :: carray
    real, intent(out) :: rms

    integer :: i, idx, j, k, jkappa, tkappa
    integer :: cheblength, nchebs
    real, dimension(:), allocatable :: u_vects
    real, dimension(:), allocatable :: core
    real, dimension(:), allocatable :: coef_u_vects
    real, dimension(ndat) :: e_sop
    real :: serieval

    ! Allocate parameter arrays

    nchebs = sum(gdim) * tdim
    allocate(coef_u_vects(nchebs))
    allocate(core(product(gdim)))

    coef_u_vects = carray(:nchebs)
    core = carray(nchebs + 1:)

    ! Allocate Chebyshev vectors

    cheblength = sum(gdim)
    allocate(u_vects(cheblength))

    ! Main loop to compute SOP energy
    
    e_sop = 0.0

    !$OMP PARALLEL DO DEFAULT(PRIVATE) SHARED(core, ndim, gdim, g_ref, tdim, coef_u_vects, e_sop)
    do idx = 1, ndat
      
    u_vects = 0.0

    ! Initialize factor vectors

    jkappa = 0
    tkappa = 0
    do i = 1, ndim
       do j = 1, gdim(i)
          jkappa = jkappa + 1
          serieval = 0.0
          do k = 1, tdim
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

end subroutine rho

function tucker_vect(ndim, gdim, u_vects, core) result(prodoned)
    implicit none
    integer, intent(in) :: ndim
    integer, dimension(ndim), intent(in) :: gdim
    real, dimension(:), intent(in) :: u_vects
    real, dimension(:), intent(in) :: core
    real :: prodoned

    integer :: mode, newsize, chebslice
    real, dimension(:), allocatable :: tensor_holder
    real, dimension(:), allocatable :: tensor_prod

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

real recursive function chebpoly(x_point, degree) result(polval)
    real, intent(in) :: x_point
    integer, intent(in) :: degree

    if (degree == 0) then
        polval = 1
    else if (degree == 1) then
        polval = x_point
    else
        polval = 2 * x_point * chebpoly(x_point, degree - 1) - chebpoly(x_point, degree - 2)
    endif
end function chebpoly

function n_mode(mode, gdim, newsize, core, spp) result(tenprod)
    integer, intent(in) :: mode
    integer, intent(in) :: newsize
    integer, dimension(:), intent(in) :: gdim
    real, dimension(:), intent(in) :: spp 
    real, dimension(:), intent(in) :: core 
    real, dimension(newsize) :: tenprod

    integer, dimension(2) :: newshape

    if (size(gdim(mode:)) == 1) then
        tenprod = dot_product(core, spp)
    else
        newshape = (/newsize, gdim(mode)/)
        tenprod = matmul(reshape(core, newshape), spp)
    end if

end function n_mode

end module soputils
