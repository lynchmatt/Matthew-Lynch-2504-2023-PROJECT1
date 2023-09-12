#############################################################################
#############################################################################
#
# This file contains units tests for polynomial sparse operations
#                                                                               
#############################################################################
#############################################################################


"""
Test product of polynomialsparses
"""
function prod_test_polysparse(;N::Int = 10^3, N_prods::Int = 20, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        prod = p1*p2
        @assert leading(prod) == leading(p1)*leading(p2)
    end

    for _ in 1:N
        p_base = PolynomialSparse([Term(1,0)])
        for _ in 1:N_prods
            p = rand(PolynomialSparse)
            prod = p_base*p
            @assert leading(prod) == leading(p_base)*leading(p)
            p_base = prod
        end
    end
    println("prod_test_polysparse - PASSED")
end


"""
Test derivative of polynomialsparses (as well as product).
"""
function prod_derivative_test_polysparse(;N::Int = 10^2,  seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_polysparse - PASSED")
end


"""
Test division of polynomialsparse modulo p.
"""
function division_test_polysparse(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        p_prod = p1*p2
        q, r = PolynomialSparse(), PolynomialSparse()
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
    println("division_test_polysparse - PASSED")
end


"""
Test the extended euclid algorithm for polynomialsparse modulo p.
"""
function ext_euclid_test_polysparse(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialSparse)
        p2 = rand(PolynomialSparse)
        g, s, t = extended_euclid_alg(p1, p2, prime)
        @assert mod(s*p1 + t*p2 - g, prime) == 0
    end
    println("ext_euclid_test_polysparse - PASSED")
end


""" 
Test excess zeroes are not stored in sparse polynomials
"""
function test_excess_zeros()
    p1 = PolynomialSparse([Term(1,2), Term(2,3), Term(2,300)])
    @assert length(p1) == 3
    println("test_excess_zeros - PASSED")
end

""" 
Test that adding the zero term does not change the length of the polynomialsparse
"""
function test_zero_addition()
    p1 = rand(PolynomialSparse)
    p1plus = p1 + Term(0,0)
    @assert length(p1plus) == length(p1)
    println("test_excess_zeros - PASSED")
end

prod_test_polysparse()
prod_derivative_test_polysparse()
division_test_polysparse()
ext_euclid_test_polysparse()
test_excess_zeros()
test_zero_addition()