#############################################################################
#############################################################################
#
# This file defines the polynomialmodp type with several operations 
#                                                                               
#############################################################################
#############################################################################

##########################################
# Polynomialsparse type and construction #
##########################################
"""
PolynomialModP type - takes a PolynomialSparse and a prime as fields, and does operations mod prime
"""

struct PolynomialModP
    # fields
    polynomial::PolynomialSparse
    prime::Int
    # no zero case?
    # inner constructor
    function PolynomialModP(poly::PolynomialSparse, p::Int)
        @assert isprime(p)
        if isempty(poly)
            poly = PolynomialSparse()
        end
        return new(mod(poly,p),p)
    end
end

"""
Construct a polynomialmod p  with a single term and prime
"""
PolynomialModP(t::Term, p::Int) = PolynomialModP(PolynomialSparse(t), p)

"""
Construct a polynomialmodp of the form x^n-x.
"""
cyclotonic_polynomialmodp(n::Int, p::Int) = PolynomialModP(PolynomialSparse([Term(1,n), Term(-1,1)]), p)


"""
Construct a polynomialmodp of the form x-n.
"""
linear_monic_polynomialmodp(n::Int, p::Int) = PolynomialModP(PolynomialSparse([Term(1,1), Term(-n,0)]), p)

"""
Construct a polynomial of the form x, mod p
"""
x_polymodp(p::Int) = PolynomialModP(PolynomialSparse(Term(1,1)), p)

"""
Creates the zero polynomial, mod p
"""
zero(::Type{PolynomialModP}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse(), p)

"""
Creates the unit polynomial with p attached
"""
one(::Type{PolynomialModP}, p::Int)::PolynomialModP = PolynomialModP(PolynomialSparse(one(Term)), p)
one(p::PolynomialModP) = one(typeof(p),p.prime)

"""
Generates a random polynomial with a random p
"""
rand(::Type{PolynomialModP}, p::Int)::PolynomialModP = PolynomialModP(rand(PolynomialSparse), p) 


###########
# Display #
###########

"""
Show a polynomial mod p
"""
function show(io::IO, p::PolynomialModP)
    show(io, p.polynomial)
    print(io, " (mod $(p.prime))")
end


##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial mod p. This implements the iteration interface.
"""
iterate(p::PolynomialModP, state=1) = iterate(p.polynomial, state)



##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialModP) = length(p.polynomial)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialModP)::Term = leading(p.polynomial) 

"""
Returns a vector of the coefficients of the polynomial.
"""
coeffs(p::PolynomialModP)::Vector{Int} = coeffs(p.polynomial)

"""
The degree of the polynomial.
"""
degree(p::PolynomialModP)::Int = degree(p.polynomial)

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialModP)::Int = euclid_alg(coeffs(p.polynomial))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialModP, x::T) where T <: Number = evaluate(f.polynomial, x)

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialModP)::Bool = iszero(p.polynomial) 

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialModP) = PolynomialModP(-(p.polynomial), p.prime)

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialModP)::PolynomialModP
    return PolynomialModP(derivative(p.polynomial), p.prime)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialModP) = ÷(p, content(p))


"""
A square free polynomial, using the prime built into polynomialmodp
"""
square_free(p::PolynomialModP)::PolynomialModP = PolynomialModP(square_free(p.polynomial, p.prime), p.prime)


#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same (in both terms and prime)
"""
==(p1::PolynomialModP, p2::PolynomialModP)::Bool = p1.polynomial == p2.polynomial && p1.prime == p2.prime


"""
Check if a polynomial is equal to 0. 
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
==(p::PolynomialModP, n::T) where T <: Real = iszero(p.polynomial)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials mod p.
"""
-(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP = p1 + (-p2)


"""
Multiplication of polynomialmodp and term.
"""
function *(t::Term, p1::PolynomialModP)::PolynomialModP  # = iszero(t) ? PolynomialSparse() : Polynomial(map((pt)->t*pt, p1.terms))
    return PolynomialModP(t*p1.polynomial, p1.prime)
end

*(p1::PolynomialModP, t::Term)::PolynomialModP = t*p1

"""
Multiplication of polynomialmodp and an integer.
"""
*(n::Int, p::PolynomialModP)::PolynomialModP = p*Term(n,0)
*(p::PolynomialModP, n::Int)::PolynomialModP = n*p

"""
Integer division of a polynomial by an integer, mod p 

Warning this may not make sense if n does not divide all the coefficients of the polynomialmodp.
"""
÷(p::PolynomialModP, n::Int)::PolynomialModP = PolynomialModP((÷(p.polynomial,n)(p.prime)), p.prime)


"""
Power of a polynomialsparse mod prime using repeated squares (for Task 6).
"""
function pow_mod(p::PolynomialModP, n::Int)
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