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
Power of a polynomialsparse128.
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
Power of a polynomialmodp
"""
function ^(p::PolynomialModP, n::Int)
    return PolynomialModP((p.polynomial^n), p.prime)
end