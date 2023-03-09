! Compile with
! gfortran -o TestODEIntegrator.x TestODEIntegrator.f90  MatrixExponential.o MatrixInverse.o Precision.o -llapack
! and then execute ./TestODEIntegrator.x
program TestODEIntegrator
    use Precision
    use LinearOperator
    use BandedMatrix
    implicit none
    
    integer :: n ! Number of ODE variables
    real(wp), allocatable :: A(:,:)
    
    type, extends(BandedMatrix) :: TestType
        integer :: test_comp
    contains
        procedure :: destroy=>TestProc  
    end type
            
    !call TestMatrixRoutines(n=2) ! Tests linear algebra in BandedMatrix class
    
    n=2
    allocate(A(n,n))
    
    lambda=1.0_wp ! Non-stiff equation
    A=A_lambda(n,lambda)
    
    
contains

subroutine TestProc(this,comp)
    class(TestType), intent(inout) :: this
    this%test_comp=comp
    call this%BandMat%destroy()
end subroutine

function A_lambda(n,lambda) result(A)
    integer, intent(in) :: n
    real(wp), intent(in) :: lambda
    real(wp) :: A(n,n) 
    
    A(1,1)=0.222222e-1_wp+0.888889_wp*lambda
    A(1,2)=0.628539e-1_wp-0.314269_wp*lambda
    A(2,1)=0.628536e-1_wp-0.314269_wp*lambda
    A(2,2)=0.177777_wp+0.111111_wp*lambda
    !write(*,*) "A=", A
    
end function

subroutine TestMatrixRoutines(n)
    integer, intent(in) :: n
 
    real(wp), dimension(n,n) :: A, test_mat, expA, invA
    real(wp), dimension(n) :: x, b, b_tmp
    real(wp) :: lambda
 
    type(BandMat) :: band_mat 
    
    band_mat = BandMat(n=n, m=n, n_bands=[-1,1], periodic=.false., id=1)

    lambda=1.0_wp           
    A = A_lambda(n,lambda)
    
    ! From Maple:
    invA(1,1) = 0.144444337752028318e1_wp
    invA(1,2) = 0.125707793940648971e1_wp
    invA(2,1) = 0.125707793940648971e1_wp
    invA(2,2) = 0.455556565221083165e1_wp
    
    call band_mat%create(n=n)
    
    call band_mat%init(matrix=A, tol=1e-6_wp, t=1.0_wp)

    !--------------------
    ! Solvers
    !--------------------       
    write(*,*) "Error in inverse = ", band_mat%invA-invA
    
    call random_number(b)
    call band_mat%solve(b=b, x=x)
    write(*,*) "Error in A^(-1)*b = ", x-matmul(invA,b)
    call band_mat%apply(x=x, b=b_tmp)
    write(*,*) "Error in residual = ", b-b_tmp
       
    !--------------------
    ! Exponentiation
    !--------------------       
    ! From Maple:
    expA(1,1) = 0.255196124690734960e1_wp
    expA(1,2) = -0.470422068844958519_wp
    expA(2,1) = -0.470422068844958519_wp
    expA(2,2) = 0.138772112224129263e1_wp

    write(*,*) "Error in exponential = ", band_mat%expA-expA

    call random_number(x)
    call band_mat%exp_t(x=x, b=b)
    write(*,*) "Error in exp(At)*b = ", b-matmul(expA,x)

    call band_mat%destroy()
    
end subroutine       

end program

