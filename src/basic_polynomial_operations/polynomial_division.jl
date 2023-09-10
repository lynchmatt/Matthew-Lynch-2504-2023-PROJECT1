#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################

"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::Polynomial, den::Polynomial)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p) # numerator mod p and denominator mod p
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = Polynomial() #empty polynomial
        prev_degree = degree(f) # prev degree is degree of numerator mod p
        while degree(f) ≥ degree(g) # while num mod p is >= denom mod p
            h = Polynomial( (leading(f) ÷ leading(g))(p) )  # here dividing leading term by leading term, multiplyng by p, putting into polynomial
            f = mod((f - h*g), p) # numerator minus h*g
            q = mod((q + h), p)  # empty plus h, mod p
            prev_degree == degree(f) && break 
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end


"""  Modular algorithm for PolynomialSparse
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::PolynomialSparse, den::PolynomialSparse)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p) # numerator mod p and denominator mod p
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparse(Term(0,0)) #empty polynomial
        delete_element!(q.terms, q.dict, 0)
        prev_degree = degree(f) # prev degree is degree of numerator mod p
        while degree(f) ≥ degree(g) # while num mod p is >= denom mod p
            h = PolynomialSparse( (leading(f) ÷ leading(g))(p) )  # here dividing leading term by leading term, multiplyng by p, putting into polynomial
            f = mod((f - h*g), p) # numerator minus h*g
            q = mod((q + h), p)  # empty plus h, mod p
            prev_degree == degree(f) && break 
            prev_degree = degree(f)
        end
        @assert iszero( mod((num  - (q*g + f)),p))
        return q, f
    end
    return division_function
end



"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::Polynomial, den::Polynomial)  = (p::Int) -> first(divide(num,den)(p))


"""
The quotient from polynomialsparse division. Returns a function of an integer.
"""
÷(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> first(divide(num,den)(p))


"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::Polynomial, den::Polynomial)  = (p::Int) -> last(divide(num,den)(p))

"""
The remainder from polynomialsparse division. Returns a function of an integer.
"""
rem(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> last(divide(num,den)(p))