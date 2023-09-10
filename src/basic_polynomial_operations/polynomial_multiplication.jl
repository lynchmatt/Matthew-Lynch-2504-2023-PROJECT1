#############################################################################
#############################################################################
#
# This file implements polynomial multiplication 
#                                                                               
#############################################################################
#############################################################################

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
    delete_element!(p_out.terms, p_out.dict, 0)
    for t in p1.terms
        new_summand = (t * p2) # produces another polysparse
        p_out = p_out + new_summand # add polysparse
    end
    return p_out
end


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