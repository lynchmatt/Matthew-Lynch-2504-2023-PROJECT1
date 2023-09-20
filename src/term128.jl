#############################################################################################
#############################################################################################
#
# This file defines the Term128 type with several operations, for use in PolynomialSparse128
#                                                                               
#############################################################################################
#############################################################################################

#################################
# Term128 type and construction #
#################################

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

"""
A Term128.
"""
struct Term128  #structs are immutable by default
    coeff::Int128
    degree::Integer
    function Term128(coeff::Integer, degree::Integer)
        coeff = Int128(coeff)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

# if coefficient is zero, automatically makes the Term128 into Term128(0,0)

"""
Creates the zero Term128.
"""
zero(::Type{Term128})::Term128 = Term128(Int128(0),Int128(0))

"""
Creates the unit Term128.
"""
one(::Type{Term128})::Term128 = Term128(Int128(1),Int128(0))

###########
# Display #
###########

function super(s)
    super_chars = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
    res = ""
    for c in s
        if c >= '0' && c <= '9'
            res *= super_chars[c - '0' + 1]
        end
    end
    return res
end

"""
Show a Term128.
"""
function show(io::IO, t::Term128)
    #define the coefficient and degree parts separately, then combine into one with io
    #Term128(1,0)
    # IF ZERO POLY PRINT 
        #first determine if constant
        coefficient, xdegree = 0,0
        if t.degree == 0
            xdegree = "" # removes the x^0
            # now checking if coefficient is one. if so, print +1 or -1 when it's the constant term
            if t.coeff == -1 
                coefficient = "-1"
            elseif t.coeff == 1
                coefficient = "1"
            else
                coefficient = t.coeff
            end
        elseif t.degree == 1 #non-constant scenario, degree is one. still need to remove coefficients if one
            xdegree = "x"
            if t.coeff == -1 
                coefficient = "-"
            elseif t.coeff == 1
                coefficient = ""
            else
                coefficient = t.coeff
            end
        else # nonconstant scenario, degree is not one or zero. still need to remove coefficients if one
            xdegree = "x$(super(string(t.degree)))"
            if t.coeff == -1 
                coefficient = "-"
            elseif t.coeff == 1
                coefficient = ""
            else
                coefficient = t.coeff
            end
        end 
        print(io, "$(coefficient)$xdegree")
    end



"""
Check if a term is 0.
"""
iszero(t::Term128)::Bool = iszero(t.coeff)


"""
Compare two Term128s.
"""
isless(t1::Term128,t2::Term128)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a Term128 at a point x.
"""
evaluate(t::Term128, x::T) where T <: Number = t.coeff * x^t.degree

#############################
# Operations with a Term128 #
#############################

"""
Add two Term128s of the same degree.
"""
function +(t1::Term128,t2::Term128)::Term128
    @assert t1.degree == t2.degree
    Term128(Int128(t1.coeff + t2.coeff), t1.degree)
end

"""
Negate a Term128.
"""
-(t::Term128,) = Term128(Int128(-t.coeff),Int128(t.degree))  

"""
Subtract two Term128s with the same degree.
"""
-(t1::Term128, t2::Term128)::Term128 = t1 + (-t2) 

"""
Multiply two Term128s.
"""
*(t1::Term128, t2::Term128)::Term128 = Term128(Int128(t1.coeff * t2.coeff), Int128(t1.degree + t2.degree))

"""
Compute the symmetric mod of a Term128 with an integer.
"""
mod(t::Term128, p::Integer) = Term128(Int128(mod(t.coeff,p)), Int128(t.degree)) # take term and what we're taking as mod. then do the coeffecient mod p, and keep the degree.

"""
Compute the derivative of a Term128.
"""
derivative(t::Term128) = Term128(Int128(t.coeff*t.degree),Int128(max(t.degree-1,0)))

"""
Divide two Term128s. Returns a function of an integer.
"""
function ÷(t1::Term128,t2::Term128) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Integer)::Term128 = Term128(Int128(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p)), Int128(t1.degree - t2.degree))
end

"""
Integer divide a Term128 by an integer.
"""
÷(t::Term128, n::Integer) = t ÷ Term128(Int128(n),Int128(0))
