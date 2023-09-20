#############################################################################
#############################################################################
#
# This file contains units tests for multiplication through the CRT
#
#############################################################################
#############################################################################

function CRT_mult_test(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for i in 1:N
        t1 = rand(PolynomialSparse128)
        t2 = rand(PolynomialSparse128)
        @assert t1*t2 == multiplication(t1,t2)
    end
    println("CRT Basic Multiplication Test - PASSED")
end

function CRT_commute_test(;N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for i in 1:N
        t1 = rand(PolynomialSparse128)
        t2 = rand(PolynomialSparse128)
        @assert multiplication(t1,t2) == multiplication(t2,t1)
    end
    print("CRT commutativity test - PASSED")
end

CRT_mult_test()
CRT_commute_test()