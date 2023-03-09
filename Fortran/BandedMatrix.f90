module BandedMatrix
    use, intrinsic :: ISO_C_BINDING ! For interfacing with C
    use Precision
    use LinearOperator
    use MatrixInverse
    use TridiagonalSolver
    use MatrixExponential
    implicit none
    private

    ! The BandedMatrix class implements a LinearOperator interface for banded matrices
    ! In this example code only fully dense, tridiagonal, and diagonal matrices are supported,
    ! but full support for all banded matrices is easy to add using an interface to LAPACK
    
    ! Abstract representation of a linear map/matrix/operator A from R^n->R^m
    type, extends(LinOp), public :: BandMat
        integer  :: n_bands(2)=-1 ! Range of band below and above the diagonal, negative for fully dense
            ! [-1,-1] fully dense, [-1:0] lower triangular, [0:-1] upper triagonal, 
            ! [0,0] means diagonal, [1:1] means tridiagonal, [2:2] pentadiagonal, etc.
        real(wp), dimension(:,:), allocatable :: A, invA, expA ! Matrix A and A^(-1) and exp(A*t)
            ! In more serious codes this would be instead an LU/Cholesky factorization or alike
        logical(c_bool) :: periodic=.false. ! Does this banded matrix come from periodic BCs in 1D    
        logical, private :: ready=.false. ! Must call constructor first
    contains
        ! Overloaded procedures from parent class LinOp
        procedure, pass(this) :: create => CreateBandMat
        procedure, pass(this) :: init => InitBandMat
        procedure, pass(this) :: destroy => DestroyBandMat
        procedure, pass(this) :: apply => ApplyBandMat
        procedure, pass(this) :: solve => SolveBandSys
        procedure, pass(this) :: exp_t => ApplyBandExp
    end type

    interface BandMat
        module procedure ConstructBandMat
    end interface
    
contains
    
function ConstructBandMat(n,m,n_bands,periodic,id) result(this) ! Create (allocate memory, open files, etc.)
   type(BandMat) :: this ! Constructed object
   integer, intent(in) :: n ! Must be supplied
   integer, intent(in), optional :: m ! Defaults to n (square matrix)
   integer, intent(in), optional :: n_bands(2) ! Number of bands, defaults to [-1,-1] 
   logical, intent(in), optional :: periodic          
   integer, intent(in), optional :: id
   
   if(present(m)) then
       if(.not.(m==n)) stop "Only square matrices (n=m) supported in CreateBandMat"
   end if
   
   ! call this%LinOp%create(n=n,m=n,id=id) ! This seems to be illegal in Fortran, not sure why
   call CreateLinOp(this,n=n,m=n,id=id) ! Use the utility routine

   if(present(n_bands)) this%n_bands=n_bands 
   if(present(periodic)) this%periodic=periodic       
   this%ready=.true.
         
end function

subroutine CreateBandMat(this,n,m,id) ! Create (allocate memory, open files, etc.)
   class(BandMat), intent(inout) :: this
   integer, intent(in) :: n ! Must be supplied
   integer, intent(in), optional :: m ! Defaults to n (square matrix)            
   integer, intent(in), optional :: id
    
   if(.not.this%ready) then
       stop "Banded matrix must be created using BandMat constructor"
   end if
   
   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       ! We only need to store the diagonal
       allocate(this%A(n,1), this%invA(n,1), this%expA(n,1))
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       allocate(this%A(n,3)) ! Only need to store 3 diagonals
   else ! dense matrix by default
       allocate(this%A(n,n), this%invA(n,n), this%expA(n,n))
   end if

   write(*,*) "Created banded matrix with id=", this%id
  
end subroutine

subroutine InitBandMat(this,matrix,tol,t) ! Initialize (set values, factorize, etc.)
   class(BandMat), intent(inout) :: this
   real(wp), intent(in), optional :: matrix(:,:) ! Data used to create the matrix   
   real(wp), intent(in), optional :: tol ! Maximum relative error tolerance
   real(wp), intent(in), optional :: t ! If present, assume ApplyExpOp is called with the same t, and precompute exp(A*t)
   
   integer :: n   
   n=this%n

   call InitLinOp(this, matrix=matrix, tol=tol, t=t)

   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       if(.not.(size(matrix,2)==1)) write(*,*) "WARNING: matrix should be of size (n,1) for diagonal matrices"
       this%A(:,1) = matrix(:,1)
       this%invA(:,1)=1.0_wp/this%A(:,1)
       this%expA(:,1)=exp(this%A(:,1)*t)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       if(.not.(size(matrix,2)==3)) write(*,*) "WARNING: matrix should be of size (n,3) for diagonal matrices"
       this%A(:,1:3) = matrix(:,1:3)
       ! No need to precompute inverse for tridiagonal
       ! For now we don't have a special function for exponential of a tridiagonal matrix
   else ! dense matrix by default
       if(.not.(size(matrix,2)==n)) write(*,*) "WARNING: matrix should be of size (n,n) for dense matrices" 
       this%A = matrix
       this%invA=mat_inv(this%A)
       if(present(t)) call mat_expm_matlab(n=n, a=this%A, e=this%expA)
   end if
   
end subroutine

subroutine ApplyBandMat(this,x,b) ! Compute b=A*x
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(out) :: b(:) ! Actual argument: b(m)

   integer :: n   
   n=this%n

   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       b=this%A(:,1)*x(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       stop "tridiagonal matrix times vector not yet implemented"
   else ! dense matrix by default
       b=matmul(this%A,x) ! Use built in intrinsic function matmul
   end if
   
end subroutine

! Child classes are free to restrict the A^(-1) and exp(A*t) operations
! to n=m, or, to compute something like a least squares type solution

subroutine SolveBandSys(this,b,x) ! Solve b=A*x to tolerance
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: b(:) ! Actual argument: b(m)
   real(wp), intent(out) :: x(:)  ! Actual argument: x(n)

   integer :: n   
   n=this%n
   
   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       x=this%invA(:,1)*b(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       call solve_tridiag(a=this%A(:,1), b=this%A(:,2), c=this%A(:,3), &
                          r=b, u=x, n=n, periodic=this%periodic)
   else ! dense matrix by default
       x=matmul(this%invA,b) ! x=A^(-1)*b
   end if
   
end subroutine

subroutine ApplyBandExp(this,x,b,t) ! Compute b=exp(A*t)*x to tolerance
   class(BandMat), intent(inout) :: this
   real(wp), intent(in) :: x(:)  ! Actual argument: x(n)
   real(wp), intent(out) :: b(:) ! Actual argument: b(m)
   real(wp), intent(in), optional :: t ! If not present, use the t set when initialized

   integer :: n   
   n=this%n

   if(all(this%n_bands==[0,0])) then ! Diagonal matrix
       b=this%expA(:,1)*x(:)
   else if((n>2).and.all(this%n_bands==[1,1])) then ! Tridiagonal matrix
       stop "Exponential of a tridiagonal matrix not yet implemented"
   else ! dense matrix by default
       if(present(t)) then
           if(.not.(t==this%t)) then
               write(*,*) "Recomputing dense expA for BandedMatrix with id=", this%id
               call mat_expm_matlab ( n=n, a=this%A, e=this%expA) ! Recompute exp(A*t)
           end if     
       end if
       b=matmul(this%expA,x) ! b=exp(A*t)*x               
   end if
   
end subroutine
        
subroutine DestroyBandMat(this) ! Does nothing but print a message to confirm destruction
  class(BandMat), intent(inout) :: this

  write(*,*) "Destroying banded matrix with id=", this%id
  
  if(allocated(this%A)) deallocate(this%A)
  if(allocated(this%invA)) deallocate(this%invA)
  if(allocated(this%expA)) deallocate(this%expA)
  
  this%ready = .false.
  
  !call this%LinOp%destroy() ! This seems to illegal in Fortran, not sure why
  call DestroyLinOp(this)

end subroutine

end module BandedMatrix

