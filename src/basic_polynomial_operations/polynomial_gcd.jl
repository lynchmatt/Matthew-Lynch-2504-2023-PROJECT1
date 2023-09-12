#############################################################################
#############################################################################
#
# This file implements polynomial GCD 
#                                                                               
#############################################################################
#############################################################################

####################
# EUCLID ALGORITHM #
####################

"""
The extended euclid algorithm for polynomialdense modulo prime.
"""
function extended_euclid_alg(a::PolynomialDense, b::PolynomialDense, prime::Int)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialDense), zero(PolynomialDense)
    old_t, t = zero(PolynomialDense), one(PolynomialDense)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t 
end


"""
The extended euclid algorithm for polynomialsparses modulo prime.
"""
function extended_euclid_alg(a::PolynomialSparse, b::PolynomialSparse, prime::Int)
    empty = PolynomialSparse(zero(Term))
    delete_element!(empty.terms, empty.dict, 0)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialSparse), empty
    old_t, t = empty, one(PolynomialSparse)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

"""
The extended euclid algorithm for polynomialdense modulo prime.
"""
function extended_euclid_alg(a::PolynomialSparse128, b::PolynomialSparse128, prime::Int128)
    empty = PolynomialSparse128(zero(Term128))
    delete_element!(empty.terms, empty.dict, 0)
    old_r, r = mod(a, prime), mod(b, prime)
    old_s, s = one(PolynomialSparse128), empty
    old_t, t = empty, one(PolynomialSparse128)

    while !iszero(mod(r,prime))
        q = first(divide(old_r, r)(prime))
        old_r, r = r, mod(old_r - q*r, prime)
        old_s, s = s, mod(old_s - q*s, prime)
        old_t, t = t, mod(old_t - q*t, prime)
    end
    g, s, t = old_r, old_s, old_t
    @assert mod(s*a + t*b - g, prime) == 0
    return g, s, t  
end

#######
# GCD #
#######

"""
The GCD of two polynomials modulo prime.
"""
gcd(a::PolynomialDense, b::PolynomialDense, prime::Int) = extended_euclid_alg(a,b,prime) |> first

"""
The GCD of two polynomialsparses modulo prime.
"""
gcd(a::PolynomialSparse, b::PolynomialSparse, prime::Int) = extended_euclid_alg(a,b,prime) |> first

"""
The GCD of two polynomialsparses modulo prime.
"""
gcd(a::PolynomialSparse128, b::PolynomialSparse128, prime::Int128) = extended_euclid_alg(a,b,prime) |> first