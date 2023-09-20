#############################################################################
#############################################################################
#
# This file contains tests for multiplication through the CRT
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

function CRT_time_test()
    println("Time comparison for small numbers.")
    p1 = PolynomialSparse128([Term128(10,10), Term128(5,5)])
    p2 = PolynomialSparse128([Term128(8,8), Term128(6,7)])
    println("Original Multiplication")
    @time p1*p2
    println("CRT Multiplication")
    @time multiplication(p1,p2)
    println("Time comparison for large numbers.")
    p3 = PolynomialSparse128([Term128(10^15,10), Term128(9^12,17)])
    p4 = PolynomialSparse128([Term128(12^12,16), Term128(8^15,10)])
    println("Original Multiplication")
    @time p3*p4
    println("CRT Multiplication")
    @time multiplication(p3,p4)
end

CRT_time_test()