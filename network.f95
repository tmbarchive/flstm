module abstract_network
    implicit none

    type, public :: network
        real :: eta = 1e-4
    contains
        procedure :: forward => forward
        procedure :: backward => backward
        procedure :: update => update
        procedure :: save => save
        procedure :: load => load
    end type network

contains
    subroutine forward(net,z,x)
        class(network) net
        real, intent(out) :: z(:,:)
        real, intent(in) :: x(:,:)
        stop "unimplemented"
    end subroutine forward

    subroutine backward(net,deltas,pdeltas)
        class(network) net
        real, intent(in) :: deltas(:,:)
        real, intent(out) :: pdeltas(:,:)
        stop "unimplemented"
    end subroutine backward

    subroutine update(net)
        class(network) net
        real eta
        stop "unimplemented"
    end subroutine update

    subroutine save(net,u)
        class(network) net
        integer u
        stop "unimplemented"
    end subroutine save

    subroutine load(net,u)
        class(network) net
        integer u
        stop "unimplemented"
    end subroutine load

end module abstract_network
