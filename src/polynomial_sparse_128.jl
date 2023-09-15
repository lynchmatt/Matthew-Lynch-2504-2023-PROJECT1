#############################################################################
#############################################################################
#
# This file defines the PolynomialSparse128 type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################
using DataStructures

"""
A PolynomialSparse128 type - designed to be for polynomials with larger integer coefficients and many nonzero terms
"""

struct PolynomialSparse128

    terms::MutableLinkedList{Term128} 
    dict::Dict{Integer, DataStructures.ListNode{Term128}} 

    #Inner constructor of the 0 polynomial
    PolynomialSparse128() = new(MutableLinkedList{Term128}(zero(Term128)), Dict{Integer, DataStructures.ListNode{Term128}}(0=>DataStructures.ListNode{Term128}(zero(Term128))))

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse128(terms::Vector{Term128})
        #Filter the vector so that there is not more than a single zero term
        terms = sort(filter((t)->!iszero(t), terms))
        if isempty(terms)
            terms = [zero(Term128)]
        end
        # initialise empty linked list and dictionary
        lst = MutableLinkedList{Term128}()
        dict = Dict{Int, DataStructures.ListNode{Term128}}()
        # create list and dictionary out of the vector of terms
        for t in terms
            insert_sorted!(lst, dict, t.degree, t)
        end

        return new(lst, dict)
    end
end

"""
This function maintains the invariant of the Polynomial type so that there are no zero terms beyond the highest
non-zero term. 
"""
function trim!(p::PolynomialSparse128)::PolynomialSparse128
    i = length(p.terms)
    while i > 0
        if iszero(p.terms[i])
            a = first(p.terms)
            a.degree
            delete_element!(p.terms, p.dict, a.degree)
        else
            nothing
        end
        i = i-1
    end
    return p
end


"""
Construct a polynomial with a single term.
"""
PolynomialSparse128(t::Term128) = PolynomialSparse128([t])

"""
Construct a polynomial of the form x^n-x.
"""
cyclotonic_polynomialsparse128(n::Integer) = PolynomialSparse128([Term128(1,n), Term128(-1,1)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomialsparse128(n::Integer) = PolynomialSparse128([Term128(1,1), Term128(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_polysparse128() = PolynomialSparse128(Term128(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialSparse128})::PolynomialSparse128 = PolynomialSparse128()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialSparse128})::PolynomialSparse128 = PolynomialSparse128(one(Term128))
one(p::PolynomialSparse128) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialSparse128} ; 
                degree::Integer = -1, 
                terms::Integer = -1, 
                max_coeff::Integer = 100, 
                mean_degree::Float64 = 5.0,
                prob_term::Float64  = 0.7,
                monic = false,
                condition = (p)->true)
        
    while true 
        _degree = degree == -1 ? rand(Poisson(mean_degree)) : degree
        _terms = terms == -1 ? rand(Binomial(_degree,prob_term)) : terms
        degrees = vcat(sort(sample(0:_degree-1,_terms,replace = false)),_degree)
        coeffs = rand(1:max_coeff,_terms+1)
        monic && (coeffs[end] = 1)
        p = PolynomialSparse128( [Term128(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
        condition(p) && return p
    end
end

# ###########
# # Display #
# ###########
lowest_to_highest = false

"""
Show a polynomial.
"""
function show(io::IO, p::PolynomialSparse128)
    if iszero(p)
        print(io, "0")
    else
    # make a local variable, false if lowest_to_highest is false or doesnt exist, and true otherwise
    global localorder = true
    if (@isdefined lowest_to_highest) == false
        localorder = false
    elseif lowest_to_highest == false
        localorder = false
    else
        localorder = true
    end
        for (i,t) in (localorder ? enumerate(p.terms) : enumerate(reverse((p.terms))))
            if !iszero(p.terms[i])
                if i == 1 # if first term, only print sign if negative
                    print(io, t.coeff > 0 ? "$(string(t)[1:end])" : "- $(string(t)[2:end])")
                else # if not first term, print plus or minus
                    print(io, t.coeff < 0 ? " - $(string(t)[2:end])" : " + $(string(t)[1:end])")
                end
            end
        end
    end
end



# ##############################################
# # Iteration over the terms of the polynomial #
# ##############################################

"""
Allows to do iteration over the non-zero terms of the polynomial. This implements the iteration interface.
"""
iterate(p::PolynomialSparse128, state=1) = iterate(p.terms, state)



# ##############################
# # Queries about a polynomial #
# ##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse128) = length(p.terms)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialSparse128)::Term128 = isempty(p.terms) ? zero(Term128) : maximum(p.terms) 

"""
Returns a vector of the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparse128)::Vector{Integer} = [t.coeff for t in p.terms]

"""
The degree of the polynomial.
"""
degree(p::PolynomialSparse128)::Integer = leading(p).degree

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse128)::Integer = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialSparse128, x::T) where T <: Number = sum(evaluate(t,x) for t in f.terms)

# ################################
# # Pushing and popping of terms #
# ################################

"""
Push a new term into the polynomial.
"""
#Note that ideally this would throw an error if pushing another term of degree that is already in the polynomial
# function push!(p::Polynomial, t::Term) 
#     if t.degree <= degree(p)
#         p.terms[t.degree + 1] = t
#     else
#         append!(p.terms, zeros(Term, t.degree - degree(p)-1))
#         push!(p.terms, t)
#     end
#     return p        
# end

# """
# Pop the leading term out of the polynomial. When polynomial is 0, keep popping out 0.
# """
# function pop!(p::Polynomial)::Term 
#     popped_term = pop!(p.terms) #last element popped is leading coefficient

#     while !isempty(p.terms) && iszero(last(p.terms))
#         pop!(p.terms)
#     end

#     if isempty(p.terms)
#         push!(p.terms, zero(Term))
#     end

#     return popped_term
# end

"""
Check if the polynomial is zero.
"""
iszero(p::PolynomialSparse128)::Bool = (p.terms == MutableLinkedList{Term128}(zero(Term128))) || (p.terms == MutableLinkedList{Term128}()) 

# #################################################################
# # Transformation of the polynomial to create another polynomial #
# #################################################################

"""
The negative of a polynomial.
"""

-(p::PolynomialSparse128) = PolynomialSparse128(-collect(p.terms))

"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparse128)::PolynomialSparse128
    deriv_vector = Vector{Term128}(undef, length(p.terms))
    for (i,t) in enumerate(p.terms)
        deriv_vector[i] = derivative(t)
    end
    return PolynomialSparse128(deriv_vector)
end

"""
The prim part (multiply a polynomial by the inverse of its content).
"""
prim_part(p::PolynomialSparse128) = ÷(p, content(p))


"""
A square free polynomial.
"""
square_free(p::PolynomialSparse128, prime::Integer)::PolynomialSparse128 = (÷(p, gcd(p,derivative(p),prime)))(prime)


# #################################
# # Queries about two polynomials #
# #################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialSparse128, p2::PolynomialSparse128)::Bool = p1.terms == p2.terms


"""
Check if a polynomial is equal to 0. 
"""
#Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2. HERE HERE
==(p::PolynomialSparse128, n::T) where T <: Real = iszero(p) == iszero(n)

# ##################################################################
# # Operations with two objects where at least one is a polynomial #
# ##################################################################

"""
Subtraction of two PolynomialSparse128s.
"""
-(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128 = p1 + (-p2)


"""
Multiplication of PolynomialSparse128 and term.
"""
function *(t::Term128, p1::PolynomialSparse128)::PolynomialSparse128  # = iszero(t) ? PolynomialSparse128() : Polynomial(map((pt)->t*pt, p1.terms))
    outpoly = PolynomialSparse128(Term128(0,0))
    delete_element!(outpoly.terms, outpoly.dict, 0)
    for term in p1.terms
        product = term*t
        iszero(product) ? nothing : insert_sorted!(outpoly.terms, outpoly.dict, product.degree, product)
    end
    outpoly
end

*(p1::PolynomialSparse128, t::Term128)::PolynomialSparse128 = t*p1

"""
Multiplication of PolynomialSparse128 and an integer.
"""
*(n::Integer, p::PolynomialSparse128)::PolynomialSparse128 = p*Term128(n,0)
*(p::PolynomialSparse128, n::Integer)::PolynomialSparse128 = n*p

"""
Integer division of a polynomial by an integer.

Warning this may not make sense if n does not divide all the coefficients of p.
"""
÷(p::PolynomialSparse128, n::Integer) = (prime)->PolynomialSparse128(map((pt)->((pt ÷ n)(prime)), collect(p.terms)))


"""
Take the mod of a PolynomialSparse128 with an integer.
"""
function mod(f::PolynomialSparse128, p::Integer)::PolynomialSparse128
    mod_vector = Vector{Term128}(undef, length(f.terms))
    for (i,t) in enumerate(f.terms)
        mod_vector[i] = mod(t, p)
    end
    return PolynomialSparse128(mod_vector)
end

"""
Power of a PolynomialSparse128 mod prime.
"""
function pow_mod(p::PolynomialSparse128, n::Integer, prime::Integer)
    n < 0 && error("No negative power")
    out = one(p) # unit polynomial
    for _ in 1:n # up the value of integer n
        out *= p # multiply 
        out = mod(out, prime)
    end
    return out
end