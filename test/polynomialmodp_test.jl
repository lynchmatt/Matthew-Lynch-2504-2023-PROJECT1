#############################################################################
#############################################################################
#
# This file contains tests for polynomial sparse operations
#                                                                               
#############################################################################
#############################################################################


"""
Test product of polynomials mod  p
"""
function prod_test_polymodp(;N::Int = 10, N_prods::Int = 20, seed::Int = 0, prime::Int = 101)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        prod = p1*p2
        @assert mod(leading(prod), prime) == mod(leading(p1)*leading(p2), prime)
    end
    # need to specify prime for rand(polynomialmodp)
    for _ in 1:N
        p_base = PolynomialModP(Term(1,0), prime)
        for _ in 1:N_prods
            p = rand(PolynomialModP, prime)
            prod = p_base*p
            @assert mod(leading(prod), prime) == mod(leading(p_base)*leading(p), prime)
            p_base = prod
        end
    end
    println("prod_test_polymodp - PASSED")
end

"""
Test derivative of polynomialsmodp (as well as product).
"""
function prod_derivative_test_polymodp(;N::Int = 10^2,  seed::Int = 0, prime::Int = 101)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        p1d = derivative(p1)
        p2d = derivative(p2)
        @assert (p1d*p2) + (p1*p2d) == derivative(p1*p2)
    end
    println("prod_derivative_test_polymodp - PASSED")
end

"""
Test division of polynomialmodp
"""
function division_test_polymodp(;prime::Int = 101, N::Int = 10^4, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        p_prod = p1*p2
        q, r = PolynomialModP(zero(Term), prime), PolynomialModP(zero(Term), prime)
        try
            q, r = divide(p_prod, p2)
            r = PolynomialModP(r, prime)
            q = PolynomialModP(q, prime)
            if (q, r) == (nothing,nothing)
                println("Unlucky prime: $p1 is reduced to $(p1 % prime) modulo $prime")
                continue
            end
        catch e
            if typeof(e) == DivideError
                @assert p2 == 0
            else
                throw(e)
            end
        end
        @assert iszero(q*p2+r - p_prod)
    end
    println("division_test_polymodp - PASSED")
end

"""
Test the extended euclid algorithm for polynomialsparse modulo p.
"""
function ext_euclid_test_polymodp(;prime::Int=101, N::Int = 10^3, seed::Int = 0)
    Random.seed!(seed)
    for _ in 1:N
        p1 = rand(PolynomialModP, prime)
        p2 = rand(PolynomialModP, prime)
        g, s, t = extended_euclid_alg(p1, p2)
        g, s, t = PolynomialModP(g, prime), PolynomialModP(s, prime), PolynomialModP(t, prime)
        @assert s*p1 + t*p2 - g == 0
    end
    println("ext_euclid_test_polymodp - PASSED")
end

""" 
Test excess zeroes are not stored in polynomialmodp where there are many terms with the zero coefficient
"""
function test_excess_zeros_modp(;prime::Int = 101)
    p1 = PolynomialSparse([Term(1,2), Term(2,3), Term(2,300)])
    p1 = PolynomialModP(p1, prime)
    @assert length(p1) == 3
    println("test_excess_zeros - PASSED")
end

""" 
Test that adding the zero term does not change the length of the polynomialmodp
"""
function test_zero_addition_modp(;prime::Int=101)
    p1 = rand(PolynomialModP, prime)
    p1plus = p1 + Term(0,0)
    @assert length(p1plus) == length(p1)
    println("test_excess_zeros - PASSED")
end

prod_test_polymodp()
prod_derivative_test_polymodp()
division_test_polymodp()
ext_euclid_test_polymodp()
test_excess_zeros_modp()
test_zero_addition_modp()