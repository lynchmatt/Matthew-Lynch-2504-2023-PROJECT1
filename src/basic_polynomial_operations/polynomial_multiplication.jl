#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

##################
# MULTIPLICATION #
##################

"""
Multiply two polynomials.
"""
function *(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p_out = PolynomialDense()
    for t in p1
        new_summand = (t * p2)
        p_out = p_out + new_summand
    end
    return p_out
end


"""
Multiply two polynomialsparses.
"""
function *(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    p_out = PolynomialSparse(Term(0,0))
    delete_element!(p_out.terms, p_out.dict, 0) # create empty polysparse to 'fill'
    for t in p1.terms
        new_summand = (t * p2) # produces another polysparse
        p_out = p_out + new_summand # add polysparses
    end
    return p_out
end

"""
Multiply two polynomialsparse128s.
"""
function *(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128
    p_out = PolynomialSparse128(Term128(0,0))
    delete_element!(p_out.terms, p_out.dict, 0) # create empty polysparse to 'fill'
    for t in p1.terms
        new_summand = (t * p2) # produces another polysparse
        p_out = p_out + new_summand # add polysparses
    end
    return p_out
end

"""
Multiply two polynomialmodps, provided they are mod the same prime
"""
function *(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    return PolynomialModP((p1.polynomial*p2.polynomial), p1.prime)
end

"""
Multiply two polynomialmodp128s, provided they are mod the same prime
"""
function *(p1::PolynomialModP128, p2::PolynomialModP128)::PolynomialModP128
    @assert p1.prime == p2.prime
    return PolynomialModP128((p1.polynomial*p2.polynomial), p1.prime)
end

#########
# POWER #
#########

"""
Power of a polynomial.
"""
function ^(p::PolynomialDense, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

"""
Power of a polynomialsparse.
"""
function ^(p::PolynomialSparse, n::Int)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

"""
Original PolySparse128 power
"""
function ^(p::PolynomialSparse128, n::Integer)
    n < 0 && error("No negative power")
    out = one(p)
    for _ in 1:n
        out *= p
    end
    return out
end

"""
Original power of a polynomialmodp
"""
function ^(p::PolynomialModP, n::Int)
    return PolynomialModP((p.polynomial^n), p.prime)
end

"""
Implements the repeated squares method of powers for PolynomialSparse128
"""
function repsq_power(p::PolynomialSparse128, n::Integer)
    # find max binary power needed to reach exponent
    maxpower = Int(trunc(log2(n)))
    exponents = [2^i for i in 0:maxpower]
    binary_string = digits(n, base=2, pad=maxpower) # convert to binary string
    outpoly = one(Term128)
    for i in 1:length(exponents)
        if binary_string[i] == 1
            outpoly *= ^(p, exponents[i])
        else
            nothing
        end
    end
    return outpoly
end

"""
Implements the repeated squares method of powers (mod p) for PolynomialSparse128
"""
function repsq_pow_mod(p::PolynomialSparse128, n::Integer, prime::Integer)
    # find max binary power needed to reach exponent
    maxpower = Int(trunc(log2(n)))
    exponents = [2^i for i in 0:maxpower]
    binary_string = digits(n, base=2, pad=maxpower) # convert to binary string
    outpoly = one(Term128)
    for i in 1:length(exponents)
        if binary_string[i] == 1
            outpoly *= pow_mod(p, exponents[i], prime)
        else
            nothing
        end
    end
    return outpoly
end

"""
Implements the repeated squares method of powers for PolynomialModP
"""
function repsq_power(p::PolynomialModP, n::Integer)
    # find max binary power needed to reach exponent
    maxpower = Int(trunc(log2(n)))
    exponents = [2^i for i in 0:maxpower]
    binary_string = digits(n, base=2, pad=maxpower) # convert to binary string
    outpoly = one(Term)
    for i in 1:length(exponents)
        if binary_string[i] == 1
            outpoly *= ^(p, exponents[i])
        else
            nothing
        end
    end
    return outpoly
end