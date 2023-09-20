####################################################################################
####################################################################################
#
# This file contains units tests for polynomial factorization, both sparse and dense
#                                                                               
####################################################################################
####################################################################################
using Random

"""
Test factorization of polynomialdense
"""
function factor_test_polydense(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialDense)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_polydense - PASSED")
end


"""
Test factorization of polynomialsparses
"""
function factor_test_polysparse(;N::Int = 10, seed::Int = 0, primes::Vector{Int} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end

    println("\nfactor_test_polysparse - PASSED")
end

"""
Test factorization of polynomialsparse128
"""
function factor_test_polysparse128(;N::Int = 10, seed::Int = 0, primes::Vector{} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = rand(PolynomialSparse128)
            factorization = factor(p, prime)
            pr = mod(expand_factorization(factorization),prime)
            @assert mod(p-pr,prime) == 0 
        end
    end
    println("\nfactor_test_polysparse128 - PASSED")
end

"""
Test factorization of polynomialmodp
"""
function factor_test_polymodp(;N::Int = 10, seed::Int = 0, primes::Vector{} = [5,17,19])
    Random.seed!(seed)
    for prime in primes
        print("\ndoing prime = $prime \t")
        for _ in 1:N
            print(".")
            p = PolynomialModP(rand(PolynomialSparse), prime)
            factorization = factor(p)
            pr = expand_factorization(factorization)
            @assert iszero(p-pr) == true
        end
    end
    println("\nfactor_test_polymodp - PASSED")
end

factor_test_polymodp()