#############################################################################
#############################################################################
#
# This file implements polynomial division 
#                                                                               
#############################################################################
#############################################################################


####################
# DIVIDE ALGORITHM #
####################


"""  Modular algorithm.
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::PolynomialDense, den::PolynomialDense)
    function division_function(p::Int)
        f, g = mod(num,p), mod(den,p) # numerator mod p and denominator mod p
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialDense() #empty polynomialdense
        prev_degree = degree(f) # prev degree is degree of numerator mod p
        while degree(f) ≥ degree(g) # while num mod p is >= denom mod p
            h = PolynomialDense( (leading(f) ÷ leading(g))(p) )  # here dividing leading term by leading term, multiplyng by p, putting into polynomial
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

"""  Modular algorithm for PolynomialSparse128
f divide by g

f = q*g + r

p is a prime
"""
function divide(num::PolynomialSparse128, den::PolynomialSparse128)
    function division_function(p::Integer)
        f, g = mod(num,p), mod(den,p) # numerator mod p and denominator mod p
        degree(f) < degree(num) && return nothing 
        iszero(g) && throw(DivideError())
        q = PolynomialSparse128(Term128(0,0)) #empty polynomial
        delete_element!(q.terms, q.dict, 0)
        prev_degree = degree(f) # prev degree is degree of numerator mod p
        while degree(f) ≥ degree(g) # while num mod p is >= denom mod p
            h = PolynomialSparse128( (leading(f) ÷ leading(g))(p) )  # here dividing leading term by leading term, multiplyng by p, putting into polynomial
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

#############
# QUOTIENT  #
#############

"""
The quotient from polynomial division. Returns a function of an integer.
"""
÷(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> first(divide(num,den)(p))


"""
The quotient from polynomialsparse division. Returns a function of an integer.
"""
÷(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> first(divide(num,den)(p))

"""
The quotient from polynomialsparse division. Returns a function of an integer.
"""
÷(num::PolynomialSparse128, den::PolynomialSparse128)  = (p::Integer) -> first(divide(num,den)(p))

##############
# REMAINDER  #
##############

"""
The remainder from polynomial division. Returns a function of an integer.
"""
rem(num::PolynomialDense, den::PolynomialDense)  = (p::Int) -> last(divide(num,den)(p))

"""
The remainder from polynomialsparse division. Returns a function of an integer.
"""
rem(num::PolynomialSparse, den::PolynomialSparse)  = (p::Int) -> last(divide(num,den)(p))

"""
The remainder from polynomialsparse128 division. Returns a function of an integer.
"""
rem(num::PolynomialSparse128, den::PolynomialSparse128)  = (p::Integer) -> last(divide(num,den)(p))