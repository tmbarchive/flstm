module ftester

    use abstract_network, only: network
    class(network), allocatable :: net
    integer :: nforward = 0
    integer :: nbackward = 0

contains

    subroutine init_logistic(no,ni)
        use logregnet, only: logreg, init
        integer no,ni
        allocate(net,source=logreg(no,ni))
        print *,"init_logistic",no,ni
        ntrain = 0; ntest = 0; prederr = 0
    end subroutine init_logistic

    subroutine forward(z,x)
        real z(:,:),x(:,:)
        call net%forward(z,x)
        nforward = nforward + 1
        if (mod(nforward,1000)==0) print *,"nforward",nforward
    end subroutine forward

    subroutine backward(deltas,pdeltas)
        real deltas(:,:),pdeltas(:,:)
        call net%backward(deltas,pdeltas)
        nbackward = nbackward + 1
        if (mod(nbackward,1000)==0) print *,"nbackward",nbackward
    end subroutine backward

end module ftester
