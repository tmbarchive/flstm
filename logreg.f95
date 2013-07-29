module logregnet
    use abstract_network, only : network

    implicit none
    real :: dsigmoid_floor = 1e-6
    integer :: verbose = 1
    integer :: nmax = 5000

    type, public, extends(network) :: logreg
        integer nh,no
        real, allocatable :: zs(:,:),inputs(:,:)
        real, allocatable :: w(:,:),dw(:,:)
    contains
        procedure :: forward => forward
        procedure :: backward => backward
        procedure :: update => update
        procedure :: save => save
        procedure :: load => load
    end type logreg

    interface logreg
        module procedure init
    end interface logreg

contains

    elemental real function sigmoid(x)
        real, intent(in) :: x
        sigmoid = 1.0 / (1.0+exp(-min(max(x,-200.0),200.0)))
    end function sigmoid

    elemental real function dsigmoidy(y)
        real, intent(in) :: y
        dsigmoidy = max(y*(1-y),dsigmoid_floor)
    end function dsigmoidy

    function init(no,nh) result(net)
        type(logreg) net
        integer no,nh
        net%eta = 1e-4
        net%no = no
        net%nh = nh
        allocate(net%w(no,nh+1))
        call random_number(net%w)
        net%w = 1e-4*(net%w-0.5)
        allocate(net%dw(no,nh+1))
        net%dw = 0
        allocate(net%zs(nmax,no))
        allocate(net%inputs(nmax,nh))
    end function init

    subroutine forward(net,z,x)
        class(logreg) net
        real,intent(out) :: z(:,:)
        real,intent(in) :: x(:,:)
        integer t
        net%inputs(:,2:) = x
        net%inputs(:,1) = 1
        do t=1,size(z,1)
            z(t,:) = sigmoid(matmul(net%w,net%inputs(t,:)))
        end do
        net%zs(:size(z,1),:) = z
    end subroutine forward

    subroutine sumouter(out,a,b)
        real :: out(:,:),a(:,:),b(:,:)
        integer i,j
        do i=1,size(a,2)
            do j=1,size(b,2)
                out(i,j) = sum(a(:,i)*b(:,j))
            end do
        end do
    end subroutine sumouter

    subroutine backward(net,deltas,pdeltas)
        class(logreg) net
        real, intent(in) :: deltas(:,:)
        real, intent(out) :: pdeltas(:,:)
        integer :: n,t
        n = size(deltas,1)
        do t=n,1,-1
            associate (v => matmul(deltas(t,:)*dsigmoidy(net%zs(t,:)),net%w))
              pdeltas(t,:) = v(1:)
            end associate
        end do
        call sumouter(net%DW,deltas,net%inputs)
    end subroutine backward

    subroutine update(net)
        class(logreg) net
        real eta
        eta = net%eta
        net%W = net%W + eta * net%DW
        net%DW = 0
    end subroutine update

    subroutine save(net,u)
        class(logreg) net
        integer u
    end subroutine save

    subroutine load(net,u)
        class(logreg) net
        integer u
    end subroutine load

end module logregnet
