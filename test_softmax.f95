module helpers
    implicit none
contains
    elemental real function H(x)
        real, intent(in) :: x
        H = 0.
        if (x>0.) H = 1.
    end function H
end module helpers

program test_softmax
    use helpers
    use softmaxnet, only: softmax
    implicit none
    type(softmax) net
    integer trial
    real inputs(100,4),outputs(100,1),targets(100,1),deltas(100,1),pdeltas(100,4)
    do trial=1,100000
        call random_number(inputs)
        targets(:,1) = H(sum(inputs,dim=2)-2.0)
        call net%forward(outputs,inputs)
        deltas = targets-outputs
        if (mod(trial,1000)==0) print *,trial,sum(abs(deltas))
        call net%backward(deltas,pdeltas)
        call net%update()
    end do
end program test_softmax
