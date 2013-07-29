module lstmnet
    use abstract_network, only : network

    implicit none
    real :: dsigmoid_floor = 1e-6
    integer :: verbose = 1
    integer :: nmax = 5000

    type, public, extends(network) :: lstm
        integer ni,ns,na
        real, allocatable :: WGI(:,:),WGF(:,:),WGO(:,:),WCI(:,:),WIP(:),WFP(:),WOP(:)
        real, allocatable :: DWGI(:,:),DWGF(:,:),DWGO(:,:),DWCI(:,:),DWIP(:),DWFP(:),DWOP(:)
        real, allocatable :: cix(:,:),ci(:,:),gix(:,:),gi(:,:),gox(:,:),go(:,:),gfx(:,:),gf(:,:)
        real, allocatable :: source(:,:),state(:,:),output(:,:),outerr(:,:)
        real, allocatable :: gierr(:,:),gferr(:,:),goerr(:,:),cierr(:,:),stateerr(:,:),sourceerr(:,:)
    contains
        procedure :: forward => forward
        procedure :: backward => backward
        procedure :: update => update
        procedure :: save => save
        procedure :: load => load
    end type lstm

contains

    elemental real function sigmoid(x)
        real, intent(in) :: x
        sigmoid = 1.0 / (1.0+exp(-min(max(x,-200.0),200.0)))
    end function sigmoid

    elemental real function dsigmoidy(y)
        real, intent(in) :: y
        dsigmoidy = max(y*(1-y),dsigmoid_floor)
    end function dsigmoidy

    elemental real function ffunc(x)
        real, intent(in) :: x
        ffunc = sigmoid(x)
    end function ffunc
    elemental real function fprime(y)
        real, intent(in) :: y
        fprime = y*y
    end function fprime
    elemental real function gfunc(x)
        real, intent(in) :: x
        gfunc = sigmoid(x)
    end function gfunc
    elemental real function gprime(y)
        real, intent(in) :: y
        gprime = y*y
    end function gprime
    elemental real function hfunc(x)
        real, intent(in) :: x
        hfunc = sigmoid(x)
    end function hfunc
    elemental real function hprime(y)
        real, intent(in) :: y
        hprime = y*y
    end function hprime

    subroutine init(net,ni,ns)
        class(lstm) net
        integer ni,ns,na
        na = 1+ni+ns
        net%ni = ni
        net%ns = ns
        net%na = na
#define F(a,b) allocate(net%a(ns,na)); call random_number(net%a); allocate(net%b(ns,na))
        F(WGI,DWGI)
        F(WGF,DWGF)
        F(WGO,DWGO)
        F(WCI,DWCI)
#undef F
#define F(a,b) allocate(net%a(ns)); call random_number(net%a); allocate(net%b(ns))
        F(WIP,DWIP)
        F(WFP,DWFP)
        F(WOP,DWOP)
#undef F
#define F(a) allocate(net%a(nmax,ns)); call random_number(net%a)
        F(cix); F(ci)
        F(gix); F(gi)
        F(gox); F(go)
        F(gfx); F(gf)
#undef F
#define F(a) allocate(net%a(nmax,na)); call random_number(net%a)
        F(source)
        F(sourceerr)
#undef F
    end subroutine init

    subroutine forward(net,z,x)
        class(lstm) net
        real, intent(out) :: z(:,:)
        real, intent(in) :: x(:,:)
        real :: prev(net%ns)
        integer :: n,t
        n = size(x,1)
        do t=1,n
            prev(:) = 0
            if (t>1) prev(:) = net%output(t-1,:)
            net%source(t,1) = 1
            net%source(t,2:2+net%ni) = x(t,:)
            net%source(t,3+net%ni:) = prev(:)
            net%gix(t,:) = matmul(net%WGI,net%source(t,:))
            net%gfx(t,:) = matmul(net%WGF,net%source(t,:))
            net%gox(t,:) = matmul(net%WGO,net%source(t,:))
            net%cix(t,:) = matmul(net%WCI,net%source(t,:))
            if (t>1) then
                net%gix(t,:) = net%gix(t,:) + net%WIP(:)*net%state(t-1,:)
                net%gfx(t,:) = net%gfx(t,:) + net%WFP(:)*net%state(t-1,:)
            end if
            net%gi(t,:) = ffunc(net%gix(t,:))
            net%gf(t,:) = ffunc(net%gfx(t,:))
            net%ci(t,:) = gfunc(net%cix(t,:))
            net%state(t,:) = net%ci(t,:)*net%gi(t,:)
            if (t>1) then
                net%state(t,:) = net%state(t,:)+net%gf(t,:)*net%state(t-1,:)
                net%gox(t,:) = net%gox(t,:)+net%WOP(:)*net%state(t,:)
            end if
            net%go(t,:) = ffunc(net%gox(t,:))
            net%output(t,:) = hfunc(net%state(t,:))*net%go(t,:)
        end do
    end subroutine forward

    subroutine sumprod(out,a,b)
        real :: out(:),a(:,:),b(:,:)
        out = sum(a*b,dim=1)
    end subroutine sumprod

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
        class(lstm) net
        real, intent(in) :: deltas(:,:)
        real, intent(out) :: pdeltas(:,:)
        real :: prev(net%ns)
        integer :: n,t
        n = size(deltas,1)
        do t=n,1,-1
            net%outerr(t,:) = deltas(t,:)
            if (t<n) net%outerr(t,:) = net%outerr(t,:)+net%sourceerr(t+1,9999999)
            net%goerr(t,:) = fprime(net%go(t,:))*hfunc(net%state(t,:))*net%outerr(t,:)
            net%stateerr(t,:) = hprime(net%state(t,:)) * net%go(t,:) * net%outerr(t,:)+ net%goerr(t,:) * net%WOP(:)
            if (t<n) net%stateerr(t,:) = net%stateerr(t,:) + net%gferr(t+1,:)*net%WFP(:) + net%gierr(t+1,:)*net%WIP(:) &
                &+ net%stateerr(t+1,:)*net%gf(t+1,:)
            if (t>1) net%sourceerr(t,:) = fprime(net%gf(t,:)) * net%stateerr(t,:)*net%state(t-1,:)
            net%gierr(t,:) = fprime(net%gi(t,:))*net%stateerr(t,:)*net%ci(t,:)
            net%cierr(t,:) = gprime(net%ci(t,:))*net%stateerr(t,:)*net%gi(t,:)
            net%sourceerr(t,:) = matmul(net%gierr(t,:),net%WGI)
            if (t>0) net%sourceerr(t,:) = net%sourceerr(t,:)+matmul(net%gferr(t,:),net%WGF)
            net%sourceerr(t,:) = net%sourceerr(t,:) + matmul(net%goerr(t,:),net%WGF) + matmul(net%cierr(t,:),net%WCI)
        end do
        call sumprod(net%DWIP,net%gierr(2:n,:),net%state(1:n-1,:))
        call sumprod(net%DWFP,net%gferr(2:n,:),net%state(1:n-1,:))
        call sumprod(net%DWOP,net%goerr(2:n,:),net%state(1:n-1,:))
        call sumouter(net%DWGI,net%gierr(1:n,:),net%source(1:n,:))
        call sumouter(net%DWGF,net%gferr(1:n,:),net%source(1:n,:))
        call sumouter(net%DWGO,net%goerr(1:n,:),net%source(1:n,:))
        call sumouter(net%DWCI,net%cierr(1:n,:),net%source(1:n,:))
    end subroutine backward

    subroutine update(net)
        class(lstm) net
        real eta
        eta = net%eta
        net%WIP = net%WIP + eta * net%DWIP
        net%WFP = net%WFP + eta * net%DWFP
        net%WOP = net%WOP + eta * net%DWOP
        net%WGI = net%WGI + eta * net%DWGI
        net%WGF = net%WGF + eta * net%DWGF
        net%WGO = net%WGO + eta * net%DWGO
        net%WCI = net%WCI + eta * net%DWCI
        net%DWIP = 0; net%DWFP = 0; net%DWOP = 0
        net%DWGI = 0; net%DWGF = 0; net%DWGO = 0; net%DWCI = 0
    end subroutine update

    subroutine save(net,u)
        class(lstm) net
        integer u
    end subroutine save

    subroutine load(net,u)
        class(lstm) net
        integer u
    end subroutine load

end module lstmnet
