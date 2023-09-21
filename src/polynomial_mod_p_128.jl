#############################################################################
#############################################################################
#
# This file defines the polynomialmodp type with several operations 
#                                                                               
#############################################################################
#############################################################################

##########################################
# Polynomialmodp type and construction   #
##########################################
"""
PolynomialModP128 type - takes a PolynomialSparse128 and a prime as fields, and does operations mod prime
"""
struct PolynomialModP128
    # fields
    polynomial::PolynomialSparse128
    prime::Integer
    # inner constructor
    function PolynomialModP128(poly::PolynomialSparse128, p::Integer)
        if iszero(poly)
            poly = PolynomialSparse128()
        end
        return new(mod(poly,p),p)
    end
end

"""
Construct a polynomialmod p  with a single term and prime
"""
PolynomialModP128(t::Term128, p::Integer) = PolynomialModP128(PolynomialSparse128(t), p)

"""
Construct a polynomial of the form x, mod p
"""
x_polymodp128(p::Integer) = PolynomialModP128(PolynomialSparse128(Term128(1,1)), p)

"""
Creates the zero polynomial, mod p
"""
zero(::Type{PolynomialModP128}, p::Integer)::PolynomialModP128 = PolynomialModP128(PolynomialSparse128(), p)

"""
Creates the unit polynomial with p attached
"""
one(::Type{PolynomialModP128}, p::Integer)::PolynomialModP128 = PolynomialModP128(PolynomialSparse128(one(Term128)), p)
one(p::PolynomialModP128) = one(typeof(p),p.prime)

"""
Generates a random polynomial with a random p
"""
rand(::Type{PolynomialModP128}, p::Integer)::PolynomialModP128 = PolynomialModP128(rand(PolynomialSparse128), p) 


###########
# Display #
###########

"""
Show a polynomial mod p
"""
function show(io::IO, p::PolynomialModP128)
    show(io, p.polynomial)
    print(io, " (mod $(p.prime))")
end


##############################################
# Iteration over the terms of the polynomial #
##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial mod p. This implements the iteration Interface.
"""
iterate(p::PolynomialModP128, state=1) = iterate(p.polynomial, state)



##############################
# Queries about a polynomial #
##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialModP128) = length(p.polynomial)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialModP128)::Term128 = leading(p.polynomial) 

"""
Returns a vector of the coefficients of the polynomial.
"""
coeffs(p::PolynomialModP128)::Vector{Integer} = coeffs(p.polynomial)

"""
The degree of the polynomial.
"""
degree(p::PolynomialModP128)::Integer = degree(p.polynomial)

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialModP128)::Integer = euclid_alg(coeffs(p.polynomial))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialModP128, x::T) where T <: Number = evaluate(f.polynomial, x)

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialModP128)::Bool = iszero(p.polynomial) 

#################################################################
# Transformation of the polynomial to create another polynomial #
#################################################################

"""
The negative of a polynomial.
"""
-(p::PolynomialModP128) = PolynomialModP128(-(p.polynomial), p.prime)

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialModP128)::PolynomialModP128
    return PolynomialModP128(derivative(p.polynomial), p.prime)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialModP128) = ÷(p, content(p))


"""
A square free polynomial, using the prime built Into PolynomialModP128
"""
square_free(p::PolynomialModP128)::PolynomialModP128 = PolynomialModP128(square_free(p.polynomial, p.prime), p.prime)


#################################
# Queries about two polynomials #
#################################

"""
Check if two polynomials are the same (in both terms and prime)
"""
==(p1::PolynomialModP128, p2::PolynomialModP128)::Bool = p1.polynomial == p2.polynomial && p1.prime == p2.prime


"""
Check if a polynomial is equal to 0. 
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the Integer 2.
==(p::PolynomialModP128, n::T) where T <: Real = iszero(p.polynomial)

##################################################################
# Operations with two objects where at least one is a polynomial #
##################################################################

"""
Subtraction of two polynomials mod p.
"""
-(p1::PolynomialModP128, p2::PolynomialModP128)::PolynomialModP128 = p1 + (-p2)


"""
Multiplication of PolynomialModP128 and term.
"""
function *(t::Term128, p1::PolynomialModP128)::PolynomialModP128  # = iszero(t) ? PolynomialSparse128() : Polynomial(map((pt)->t*pt, p1.terms))
    return PolynomialModP128(t*p1.polynomial, p1.prime)
end

*(p1::PolynomialModP128, t::Term128)::PolynomialModP128 = t*p1

"""
Multiplication of PolynomialModP128 and an Integer
"""
*(n::Integer, p::PolynomialModP128)::PolynomialModP128 = p*Term128(n,0)
*(p::PolynomialModP128, n::Integer)::PolynomialModP128 = n*p

"""
Integer division of a polynomial by an Integer, mod p 

Warning this may not make sense if n does not divide all the coefficients of the PolynomialModP128.
"""
÷(p::PolynomialModP128, n::Integer)::PolynomialModP128 = PolynomialModP128((÷(p.polynomial,n)(p.prime)), p.prime)


"""
Power of a polynomialsparse mod prime using repeated squares (for Task 6).
"""
function pow_mod(p::PolynomialModP128, n::Integer)
    # find max binary power needed to reach exponent
    maxpower = Integer(trunc(log2(n)))
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