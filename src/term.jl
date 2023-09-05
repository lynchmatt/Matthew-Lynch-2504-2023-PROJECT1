#############################################################################
#############################################################################
#
# This file defines the Term type with several operations 
#                                                                               
#############################################################################
#############################################################################

##############################
# Term type and construction #
##############################

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one

"""
A term.
"""
struct Term  #structs are immutable by default
    coeff::Int
    degree::Int
    function Term(coeff::Int, degree::Int)
        degree < 0 && error("Degree must be non-negative")
        coeff != 0 ? new(coeff,degree) : new(coeff,0)
    end
end

"""
Creates the zero term.
"""
zero(::Type{Term})::Term = Term(0,0)

"""
Creates the unit term.
"""
one(::Type{Term})::Term = Term(1,0)

###########
# Display #
###########

"""
Show a term.
"""
function show(io::IO, t::Term)
#define the coefficient and degree parts separately, then combine into one with io
#term(1,0)
    #first determine if constant
    coefficient, xdegree = 0,0
    if t.degree == 0
        xdegree = "" # removes the x^0
        # now checking if coefficient is one. if so, print +1 or -1
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
        xdegree = "x^$(t.degree)"
        if t.coeff == -1 
            coefficient = "-"
        elseif t.coeff == 1
            coefficient = ""
        else
            coefficient = t.coeff
        end
    end 
    print(io, "$coefficient$xdegree")
end

super_chars = ['⁰', '¹', '²', '³', '⁴', '⁵', '⁶', '⁷', '⁸', '⁹']
function super(s)
    res = ""
    for c in s
        if c >= '0' && c <= '9'
            res *= super_chars[c - '0' + 1]
        end
    end
    return res
end

function show2(io::IO, t::Term)
    #define the coefficient and degree parts separately, then combine into one with io
    #term(1,0)
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
            xdegree = ""
            if t.coeff == -1 
                coefficient = "-"
            elseif t.coeff == 1
                coefficient = ""
            else
                coefficient = t.coeff
            end
        else # nonconstant scenario, degree is not one or zero. still need to remove coefficients if one
            xdegree = "$(super(string(t.degree)))"
            if t.coeff == -1 
                coefficient = "-"
            elseif t.coeff == 1
                coefficient = ""
            else
                coefficient = t.coeff
            end
        end 
        print(io, "$(coefficient)x$xdegree")
    end


#####################################
# Updated Display - Pretty Printing #
#####################################

# function show(io::IO, t::Term)
#     #define the coefficient and degree parts separately, then combine into one with io
#     #first determine if constant
#     coefficient, xdegree = 0,0
#     if t.degree == 0
#         xdegree = "" # removes the x^0
#         # now checking if coefficient is one. if so, print +1 or -1
#         if t.coeff == -1 
#             coefficient = "-1"
#         elseif t.coeff == 1
#             coefficient = "1"
#         else
#             coefficient = t.coeff
#         end
#     elseif t.degree == 1 #non-constant scenario, degree is one. still need to remove coefficients if one
#         xdegree = "x"
#         if t.coeff == -1 
#             coefficient = "-"
#         elseif t.coeff == 1
#             coefficient = ""
#         else
#             coefficient = t.coeff
#         end
#     else # nonconstant scenario, degree is not one or zero. still need to remove coefficients if one
#         xdegree = "x^$(t.degree)"
#         if t.coeff == -1 
#             coefficient = "-"
#         elseif t.coeff == 1
#             coefficient = ""
#         else
#             coefficient = t.coeff
#         end
#     end 
#     print(io, "$coefficient$xdegree")
# end 


## version without plusminus one working when degree is zero
# function show(io::IO, t::Term)
#     #define the coefficient and degree parts separately, then combine into one with io
#     coefficient, xdegree = 0,0
#     if t.coeff == -1
#         coefficient = "-"
#     elseif t.coeff == 1
#         coefficient = ""
#     else
#         coefficient = t.coeff
#     end
#     if t.degree == 0
#         xdegree = ""
#     elseif t.degree == 1
#         xdegree = "x"
#     else
#         xdegree = "x^$(t.degree)"
#     end  
#     print(io, "$coefficient$xdegree")
# end 


########################
# Queries about a term #
########################

"""
Check if a term is 0.
"""
iszero(t::Term)::Bool = iszero(t.coeff)

iszero(Term(0,0))

"""
Compare two terms.
"""
isless(t1::Term,t2::Term)::Bool =  t1.degree == t2.degree ? (t1.coeff < t2.coeff) : (t1.degree < t2.degree)  

"""
Evaluate a term at a point x.
"""
evaluate(t::Term, x::T) where T <: Number = t.coeff * x^t.degree

##########################
# Operations with a term #
##########################

"""
Add two terms of the same degree.
"""
function +(t1::Term,t2::Term)::Term
    @assert t1.degree == t2.degree
    Term(t1.coeff + t2.coeff, t1.degree)
end

"""
Negate a term.
"""
-(t::Term,) = Term(-t.coeff,t.degree)  

"""
Subtract two terms with the same degree.
"""
-(t1::Term, t2::Term)::Term = t1 + (-t2) 

"""
Multiply two terms.
"""
*(t1::Term, t2::Term)::Term = Term(t1.coeff * t2.coeff, t1.degree + t2.degree)


"""
Compute the symmetric mod of a term with an integer.
"""
mod(t::Term, p::Int) = Term(mod(t.coeff,p), t.degree)

"""
Compute the derivative of a term.
"""
derivative(t::Term) = Term(t.coeff*t.degree,max(t.degree-1,0))

"""
Divide two terms. Returns a function of an integer.
"""
function ÷(t1::Term,t2::Term) #\div + [TAB]
    @assert t1.degree ≥ t2.degree
    f(p::Int)::Term = Term(mod((t1.coeff * int_inverse_mod(t2.coeff, p)), p), t1.degree - t2.degree)
end

"""
Integer divide a term by an integer.
"""
÷(t::Term, n::Int) = t ÷ Term(n,0)
