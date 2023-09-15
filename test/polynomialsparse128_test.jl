#############################################################################
#############################################################################
#
# This file contains units tests for polynomial sparse operations
#                                                                               
#############################################################################
#############################################################################

###########################################
# POLYNOMIALSPARSE128 FUNCTIONALITY TESTS #
###########################################
"""
Test product of polynomialsparse128
"""
function prod_test_polysparse128(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparse128([Term128(1,0)])
        for _ in 1:N_prods
            p = rand(PolynomialSparse128)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_polysparse128 - PASSED")
end


"""
Test derivative of polynomialsparse128 (as well as product).
"""
function prod_derivative_test_polysparse128(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_polysparse128 - PASSED")
end


"""
Test division of polynomialsparse128 modulo p.
"""
function division_test_polysparse128(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        p_prod = p1*p2
        q, r = PolynomialSparse128(), PolynomialSparse128()
        try
            q, r = divide(p_prod, p2)(prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert mod(p2, prime) == 0
            else
                throw(e)
            end
        end
        @assert iszero( mod(q*p2+r - p_prod, prime) )
    end
    println("division_test_polysparse128 - PASSED")
end


"""
Test the extended euclid algorithm for polynomialsparse128 modulo p.
"""
function ext_euclid_test_polysparse128(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse128)
        p2 = rand(PolynomialSparse128)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_polysparse128 - PASSED")
end


""" 
Test excess zeroes are not stored in polynomialsparse128
"""
function test_excess_zeros_128()
    p1 = PolynomialSparse128([Term128(1,2), Term128(2,3), Term128(2,300)])
    @assert length(p1) == 3
    println("test_excess_zeros_128 - PASSED")
end

""" 
Test that adding the zero term does not change the length of the polynomialsparse128
"""
function test_zero_addition_128()
    p1 = rand(PolynomialSparse128)
    p1plus = p1 + Term128(0,0)
    @assert length(p1plus) == length(p1)
    println("test_zero_addition_128 - PASSED")
end

"""
Test that when coefficients are sufficiently large, PolynomialSparse will overflow and yield incorrect answers, while PolynomialSparse128 will not
"""
function overflow_test()
    sparsepoly = PolynomialSparse([Term(100,2), Term(100,1), Term(100,0)])
    sparsepoly128 = PolynomialSparse128([Term128(100,2), Term128(100,1), Term128(100,0)])
    # the coefficients of both should be a multiple of 10, when raised to any power. check if polynomialsparse has overflown
    if iszero(mod(sparsepoly^10, 10))
        println("PolynomialSparse has not overflown, insufficient example.")
    else
    # check that the coefficients of the polynomialsparse128 are in fact multiples of ten, when raised to the tenth power
        if iszero(mod(sparsepoly128^10, 10))
        else
            throw(DomainError("OVERFLOW_TEST - FAILED"))
        end
    end
    # the following coefficient is the maximum number included in Int64 - show that it does not overflow when we add the unit term
    maxpolysparse = PolynomialSparse128([Term128(9223372036854775807,1)])
    maxpolysparse += Term128(1,1)
    # show that Int64 would normally overflow at this point
    @assert Int64(9223372036854775807) + Int64(1) != Int128(9223372036854775808)
    # show that since Int128 does not overflow at this point, PolynomialSparse128 also does not overflow
    @assert maxpolysparse.terms[1] == Term128(9223372036854775808,1)
    println("OVERFLOW_TEST - PASSED.")
end

prod_test_polysparse128()
prod_derivative_test_polysparse128()
division_test_polysparse128()
ext_euclid_test_polysparse128()
test_excess_zeros_128()
test_zero_addition_128()
overflow_test()

