!!! reorder hidden units by average activation (to make growing/shrinking work better
!!! handle automatic class mapping... sparse vectors?

module mlps

    use utils
    use dataio
    use quicksort

    implicit none

    type mlp
        integer ninput,nhidden,noutput
        real, allocatable :: b1(:), w1(:,:), b2(:), w2(:,:)
        integer updates
        real delta1,delta2
        real eta
    end type mlp

    type autoparams
        integer :: stagesize = 200000
        integer :: nstages = 50
        integer :: nnets = 16
        real :: hidden_lo = 15.0
        real :: hidden_hi = 150.0
        real :: hidden_max = 1000.0
        real :: cv_frac = 0.1
        real :: sizecost = 0.01
    end type autoparams

    real :: dsigmoid_floor = 0.0

    integer :: verbose = 1
contains

    subroutine mlp_write(net,u)
        type(mlp) net
        integer u
        write(unit=u) net%ninput,net%nhidden,net%noutput
        write(unit=u) net%b1,net%w1,net%b2,net%w2
    end subroutine mlp_write

    subroutine mlp_read(net,u)
        type(mlp) net
        integer u
        read(unit=u) net%ninput,net%nhidden,net%noutput
        read(unit=u) net%b1,net%w1,net%b2,net%w2
    end subroutine mlp_read

    pure real function sigmoid(x)
        real, intent(in) :: x
        sigmoid = 1.0 / (1.0+exp(-min(max(x,-200.0),200.0)))
    end function sigmoid

    pure real function dsigmoidy(y)
        real, intent(in) :: y
        dsigmoidy = max(y*(1-y),dsigmoid_floor)
    end function dsigmoidy

    subroutine mlp_init(net,ninput,nhidden,noutput)
        type(mlp) net
        integer ninput,nhidden,noutput
        if (verbose>=10) print *,"mlp nin",ninput,"nhid",nhidden,"nout",noutput
        net%ninput = ninput
        net%nhidden = nhidden
        net%noutput = noutput
        allocate(net%w1(nhidden,ninput))
        allocate(net%b1(nhidden))
        allocate(net%w2(noutput,nhidden))
        allocate(net%b2(noutput))
        call random_number(net%w1)
        call random_number(net%b1)
        call random_number(net%w2)
        call random_number(net%b2)
        net%w1 = 2*net%w1 - 1
        net%b1 = 2*net%b1 - 1
        net%w2 = 2*net%w2 - 1
        net%b2 = 2*net%b2 - 1
        net%eta = -1.0
        net%updates = 0
    end subroutine mlp_init

    subroutine mlp_check(net)
        type(mlp) net
        if (net%ninput<2 .or. net%ninput>1000000) stop "bad net ninput"
        if (net%nhidden<2 .or. net%nhidden>1000000) stop "bad net nhidden"
        if (net%noutput<2 .or. net%noutput>1000000) stop "bad net noutput"
        if (size(net%w1,1)/=net%nhidden) stop "net oops"
        if (size(net%b1)/=net%nhidden) stop "net oops"
        if (size(net%w1,2)/=net%ninput) stop "net oops"
        if (size(net%w2,1)/=net%noutput) stop "net oops"
        if (size(net%w2,2)/=net%nhidden) stop "net oops"
        if (size(net%b2)/=net%noutput) stop "net oops"
    end subroutine mlp_check

    subroutine data_check(classes,inputs,net)
        integer, intent(in) :: classes(:)
        real, intent(in) :: inputs(:,:)
        type(mlp),optional :: net

        if (maxval(abs(inputs))>100.0) stop "input not properly normalized"
        if (maxval(abs(inputs))<1-2) stop "input not properly normalized"
        if (size(inputs,1)/=size(classes)) stop "input size mismatch"
        if (minval(classes)<1) stop "classes numbered from zero???"
        if (present(net)) then
            if (size(inputs,2)/=net%ninput) stop "input vector size mismatch"
            if (maxval(classes)>net%noutput) stop "output class larger than output vector"
        end if
    end subroutine data_check

    subroutine mlp_adjust_hidden0(net,nhidden)
        type(mlp) net,anet
        real scale1,scale2
        integer nhidden,i
        call mlp_init(anet,net%ninput,nhidden,net%noutput)
        anet%eta = net%eta
        scale1 = sum(abs(net%w1))/size(net%w1)
        scale2 = sum(abs(net%w2))/size(net%w2)
        do i=1,nhidden
            anet%w1(i,:) = net%w1(imod(i,net%nhidden),:)
            anet%b1(i) = net%b1(imod(i,net%nhidden))
            anet%w2(:,i) = net%w2(:,imod(i,net%nhidden))
            if (i>net%nhidden) then
                call perturb(anet%w1(i,:),1e-2*scale1)
                anet%w2(:,i) = randunif()*1e-2*scale2
            end if
        end do
        net = anet
    end subroutine mlp_adjust_hidden0

    subroutine mlp_adjust_hidden(net,nhidden)
        type(mlp) net,anet
        real scale1,scale2
        integer nhidden,i
        call mlp_init(anet,net%ninput,nhidden,net%noutput)
        anet%eta = net%eta
        scale1 = sum(abs(net%w1))/size(net%w1)
        scale2 = sum(abs(net%w2))/size(net%w2)
        do i=1,nhidden
            if (i<=net%nhidden) then
                anet%w1(i,:) = net%w1(i,:)
                anet%b1(i) = net%b1(i)
                anet%w2(:,i) = net%w2(:,i)
            else
                call random_number(anet%w1(i,:))
                call random_number(anet%b1(i))
                call random_number(anet%w2(:,i))
                anet%w1(i,:) = scale1 * (2*anet%w1(i,:) - 1.0)
                anet%b1(i) = scale1 * (2*anet%b1(i) - 1.0)
                anet%w2(:,i) = scale2 * 1e-2 * (2*anet%w2(:,i) - 1.0)
            end if
        end do
        net = anet
    end subroutine mlp_adjust_hidden

    subroutine mlp_forward(net,z,x)
        type(mlp) net
        real, intent(out) :: z(:)
        real, intent(in) :: x(:)
        real y(net%nhidden)
        integer i
        y = matmul(net%w1,x) + net%b1
        forall (i=1:net%nhidden) y(i) = sigmoid(y(i))
        z = matmul(net%w2,y) + net%b2
        forall (i=1:net%noutput) z(i) = sigmoid(z(i))
    end subroutine mlp_forward

    function mlp_error(net,classes,inputs) result(errs)
        type(mlp) net
        integer, intent(in) :: classes(:)
        real, intent(in) :: inputs(:,:)
        real :: z(net%noutput)
        integer :: targt
        real :: errs
        integer i,row

        if (verbose>=10) print *,"mlp_error",size(inputs,1),"samples",&
            &size(classes),maxval(classes),size(inputs,2)

        call mlp_check(net)
        call data_check(classes,inputs)

        errs = 0

        do row=1,size(inputs,1)
            call mlp_forward(net,z,inputs(row,:))
            if (maxloc(z,1) /= classes(row)) errs = errs + 1
        end do
    end function mlp_error

    subroutine mlp_confusion(cm,net,targets,inputs)
        type(mlp) net
        real, allocatable :: cm(:,:)
        integer, intent(in) :: targets(:)
        real, intent(in) :: inputs(:,:)
        real :: z(net%noutput)
        real :: targt(net%noutput)
        integer i,row,nclasses
        logical classmode

        call mlp_check(net)
        nclasses = maxval(targets)
        allocate(cm(nclasses,nclasses))
        cm = 0

        do row=1,size(inputs,1)
            call mlp_forward(net,z,inputs(row,:))
            cm(targets(row),maxloc(z,1)) = cm(targets(row),maxloc(z,1)) + 1
        end do
    end subroutine mlp_confusion

    subroutine print_confusion(cm)
        real :: cm(:,:)
        integer i,j
        print *,"--- confusion matrix ---"
        do i=1,size(cm,1)
            write (*,fmt="(i4)",advance="no") (floor(cm(i,j)),j=1,size(cm,2))
            write (*,fmt="(a)") " "
        end do
        print *,"---"
    end subroutine print_confusion

    subroutine mlp_clear_info(net)
        type(mlp) net
        net%updates = 0
        net%delta1 = 0
        net%delta2 = 0
    end subroutine mlp_clear_info

    subroutine mlp_print_info(net)
        type(mlp) net
        print *,"layer1",minval(net%w1),maxval(net%w1),minval(net%b1),maxval(net%b1)
        print *,"layer2",minval(net%w2),maxval(net%w2),minval(net%b2),maxval(net%b2)
        print *,"update",net%updates,net%delta1/net%updates,net%delta2/net%updates
    end subroutine mlp_print_info

    subroutine mlp_train(net,classes,inputs,ntrials)
        type(mlp) net
        integer, intent(in), target :: classes(:)
        real, intent(in), target :: inputs(:,:)
        real :: weights(size(inputs,1)), weight
        real x(net%ninput),targt(net%noutput)
        real z(net%noutput)
        real y(net%nhidden)
        real delta2(net%noutput),delta1(net%nhidden)
        integer i,j,row,ntrials
        logical classmode

        if (verbose>=12) print *,"mlp_train",size(inputs,1),&
            &"samples",size(classes),size(inputs,2)

        call mlp_check(net)
        call data_check(classes,inputs,net)
        if (net%eta<1e-8) stop "bad eta"

        weights = 1

        do i=1,ntrials
            row = randint(size(inputs,1))
            weight = weights(row)
            x = inputs(row,:)
            targt = 0
            targt(classes(row)) = 1

            if (maxval(x)>10 .or. minval(x)<-10) stop "should normalize mlp inputs"

            ! forward propagation
            y = matmul(net%w1,x) + net%b1
            forall (i=1:net%nhidden) y(i) = sigmoid(y(i))
            z = matmul(net%w2,y) + net%b2
            forall (i=1:net%noutput) z(i) = sigmoid(z(i))

            ! backward propagation of deltas
            forall (i=1:net%noutput) delta2(i) = (z(i)-targt(i)) * dsigmoidy(z(i))
            delta1 = matmul(delta2,net%w2)
            forall (i=1:net%nhidden) delta1(i) = delta1(i) * dsigmoidy(y(i))

            ! weight update
            forall (i=1:net%noutput,j=1:net%nhidden)
                net%w2(i,j) = net%w2(i,j) - net%eta * delta2(i) * y(j) * weight
            end forall
            forall (i=1:net%noutput) net%b2(i) = net%b2(i) - net%eta * delta2(i) * weight
            forall (i=1:net%nhidden,j=1:net%ninput)
                net%w1(i,j) = net%w1(i,j) - net%eta * delta1(i) * x(j)
            end forall
            forall (i=1:net%nhidden) net%b1(i) = net%b1(i) - net%eta * delta1(i) * weight

            net%updates = net%updates + 1
            net%delta1 = net%delta1 + maxval(abs(delta1))
            net%delta2 = net%delta2 + maxval(abs(delta2))
        end do
    end subroutine mlp_train

    subroutine mlp_decay(net,delta)
        type(mlp) net
        real delta
        integer i,j
        if (delta<0.5 .or. delta>1.0) stop "bad delta value in mlp_decay"
        forall (j=1:net%noutput,i=1:net%nhidden)
            net%w2(i,j) = net%w2(i,j) * delta
        end forall
        forall (j=1:net%nhidden,i=1:net%ninput)
            net%w1(i,j) = net%w1(i,j) * delta
        end forall
    end subroutine mlp_decay

    subroutine mlp_autotrain0(net,targets,inputs)
        type(mlp) net
        real inputs(:,:)
        integer targets(:)
        integer epoch,n
        integer :: epochs = 30
        integer, parameter :: ntrials = 4
        integer trial,best
        type(mlp) trials(ntrials)
        real :: etas(ntrials) = [0.1,0.1,0.3,0.7]
        real :: errs(ntrials,2)
        ! FIXME do more here
        n = size(targets)
        !$omp parallel do shared(targets,inputs,errs)
        do trial=1,ntrials
            print *,"trial",trial,"starting"
            call mlp_init(trials(trial),size(inputs,2),20,maxval(targets))
            do epoch=1,epochs
                call mlp_train(trials(trial),targets,inputs,100000)
            end do
            errs(trial,:) = mlp_error(trials(trial),targets,inputs)
            print *,"trial",trial,"errs",errs(trial,2)/n
        end do
        !$omp end parallel do
        best = minloc(errs(:,2),1)
        print *,"best",best,errs(best,2)/n
        net = trials(best)
    end subroutine mlp_autotrain0

    subroutine mlp_autotrain(ps,net,targets,inputs,test_targets,test_inputs)
        type(autoparams) ps
        type(mlp) net
        real inputs(:,:)
        integer targets(:)
        real test_inputs(:,:)
        integer test_targets(:)

        type(mlp) nets(ps%nnets)
        integer nsamples,nclasses,nfeat,stage,nn,i,j,nhidden
        real errs(ps%nnets)
        integer order(ps%nnets),n
        real :: best

        errs = randnormal()
        nfeat = size(inputs,2)
        nsamples = size(targets)
        nclasses = maxval(targets)

        print *,"autotrain params",ps
        print *,"autotrain nfeat",nfeat
        print *,"autotrain nsamples",nsamples
        print *,"autotrain nclasses",nclasses
        print *,"autotrain ntest",size(test_targets)

        ! initialize the networks

        best = huge(1.0)
        call mlp_init(net,1,1,1)

        do nn=1,ps%nnets
            nhidden = logspace(nn,ps%nnets,ps%hidden_lo,ps%hidden_hi)
            nhidden = min(nhidden,floor(ps%hidden_max))
            call mlp_init(nets(nn),nfeat,nhidden,nclasses)
            nets(nn)%eta = exp(randnormal()+log(0.3))
            print *,"init",nn,nhidden,nets(nn)%eta
        end do

        do stage=1,ps%nstages
            errs = -1
            !$omp parallel do shared(targets,inputs,errs)
            do nn=1,ps%nnets
                call mlp_train(nets(nn),&
                    &targets,inputs,ps%stagesize)
                errs(nn) = mlp_error(nets(nn),&
                    &test_targets,test_inputs)
                write (*,fmt="(99i6)") floor(errs(:))
            end do
            !$omp end parallel do

            do nn=1,ps%nnets
                errs(nn) = errs(nn) + nets(nn)%nhidden * ps%sizecost
            end do
            call quicksort1(order,errs(:),ps%nnets)

            print *,"stage",stage
            print *,"errs+c",(errs(order(j)),j=1,ps%nnets)
            print *,"etas  ",(nets(order(j))%eta,j=1,ps%nnets)
            print *,"nhid  ",(nets(order(j))%nhidden,j=1,ps%nnets)
            print *,"updt  ",(nets(order(j))%updates,j=1,ps%nnets)

            ! save the best network so far
            if (best<huge(1.0)) call mlp_check(net)
            if (errs(order(1))<best) then
                best = errs(order(1))
                call mlp_check(nets(order(1)))
                net = nets(order(1)) ! save the best net so far
            end if
            call mlp_check(net)

            ! replace the worst 50% with the best 50%
            do i=1,ps%nnets/2
                j = order(i+ps%nnets/2)
                nets(j) = nets(order(i))
                call mlp_check(net)
                nets(j)%eta = exp(randnormal()+log(nets(j)%eta))
                n = max(2,nets(j)%nhidden+floor(0.2*nets(j)%nhidden*(2.0*randunif()-1.0)+0.5))
                n = min(n,floor(ps%hidden_max))
                print *,"adjust",nets(j)%nhidden,n
                call mlp_adjust_hidden(nets(j),n)
                print *,"net",i,"nhid",nets(order(i))%nhidden,"eta",nets(order(i))%eta
                print *,"net",i+ps%nnets/2,"nhid",nets(j)%nhidden,"eta",nets(j)%eta
            end do
        end do
        call mlp_check(net)
    end subroutine mlp_autotrain

    subroutine mlp_autotrain_cv(ps,net,targets,inputs)
        type(autoparams) ps
        type(mlp) net
        real inputs(:,:)
        integer targets(:)
        integer selection(size(inputs,1))
        real, allocatable :: test_inputs(:,:),train_inputs(:,:)
        integer, allocatable :: test_targets(:), train_targets(:)
        integer nt,ndim,n

        n = size(selection)
        nt = floor(ps%cv_frac*n)
        ndim = size(inputs,2)

        call permutation1(selection,size(selection))
        allocate(test_inputs(nt,ndim))
        allocate(test_targets(nt))
        allocate(train_inputs(n-nt,ndim))
        allocate(train_targets(n-nt))
        test_inputs = inputs(:nt,:)
        test_targets = targets(:nt)
        train_inputs = inputs(nt+1:,:)
        train_targets = targets(nt+1:)
        call mlp_autotrain(ps,net,train_targets,train_inputs,test_targets,test_inputs)
        call mlp_check(net)
    end subroutine mlp_autotrain_cv

end module mlps
