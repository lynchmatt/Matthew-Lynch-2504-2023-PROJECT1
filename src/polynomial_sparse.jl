#############################################################################
#############################################################################
#
# This file defines the polynomialsparse type with several operations 
#                                                                               
#############################################################################
#############################################################################

####################################
# Polynomial type and construction #
####################################
using DataStructures

"""
A PolynomialSparse type - designed to be for polynomials with integer coefficients and few nonzero terms
"""

struct PolynomialSparse

    terms::MutableLinkedList{Term} 
    dict::Dict{Int, DataStructures.ListNode{Term}} 

    #Inner constructor of the 0 polynomial
    PolynomialSparse() = new(MutableLinkedList{Term}(zero(Term)), Dict{Int, DataStructures.ListNode{Term}}(0=>DataStructures.ListNode{Term}(zero(Term))))
    

    #Inner constructor of polynomial based on arbitrary list of terms
    function PolynomialSparse(terms::Vector{Term})
        #Filter the vector so that there is not more than a single zero term
        terms = sort(filter((t)->!iszero(t), terms))
        if isempty(terms)
            terms = [zero(Term)]
        end
        # initialise empty linked list and dictionary
        lst = MutableLinkedList{Term}()
        dict = Dict{Int, DataStructures.ListNode{Term}}()
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
function trim!(p::PolynomialSparse)::PolynomialSparse
    i = length(p.terms)
    while i > 0
        if iszero(p.terms[i])
            println("zero found")
            @show a = first(p.terms)
            @show a.degree
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
PolynomialSparse(t::Term) = PolynomialSparse([t])

"""
Construct a polynomial of the form x^p-x.
"""
cyclotonic_polynomialsparse(p::Int) = PolynomialSparse([Term(1,p), Term(-1,1)])


"""
Construct a polynomial of the form x-n.
"""
linear_monic_polynomialsparse(n::Int) = PolynomialSparse([Term(1,1), Term(-n,0)])

"""
Construct a polynomial of the form x.
"""
x_polysparse() = PolynomialSparse(Term(1,1))

"""
Creates the zero polynomial.
"""
zero(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse()

"""
Creates the unit polynomial.
"""
one(::Type{PolynomialSparse})::PolynomialSparse = PolynomialSparse(one(Term))
one(p::PolynomialSparse) = one(typeof(p))

"""
Generates a random polynomial.
"""
function rand(::Type{PolynomialSparse} ; 
                degree::Int = -1, 
                terms::Int = -1, 
                max_coeff::Int = 100, 
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
        p = PolynomialSparse( [Term(coeffs[i],degrees[i]) for i in 1:length(degrees)] )
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
function show(io::IO, p::PolynomialSparse)
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
                # print the term. if next term's coefficient is negative, print minus, if not print plus.
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
iterate(p::PolynomialSparse, state=1) = iterate(p.terms, state)



# ##############################
# # Queries about a polynomial #
# ##############################

"""
The number of terms of the polynomial.
"""
length(p::PolynomialSparse) = length(p.terms)

"""
The leading term of the polynomial.
"""
leading(p::PolynomialSparse)::Term = isempty(p.terms) ? zero(Term) : maximum(p.terms) 

"""
Returns a vector of the coefficients of the polynomial.
"""
coeffs(p::PolynomialSparse)::Vector{Int} = [t.coeff for t in p.terms]

"""
The degree of the polynomial.
"""
degree(p::PolynomialSparse)::Int = leading(p).degree

"""
The content of the polynomial is the GCD of its coefficients.
"""
content(p::PolynomialSparse)::Int = euclid_alg(coeffs(p))

"""
Evaluate the polynomial at a point `x`.
"""
evaluate(f::PolynomialSparse, x::T) where T <: Number = sum(evaluate(t,x) for t in f.terms)

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
iszero(p::PolynomialSparse)::Bool = p.terms == MutableLinkedList{Term}(zero(Term))

# #################################################################
# # Transformation of the polynomial to create another polynomial #
# #################################################################

"""
The negative of a polynomial.
"""

-(p::PolynomialSparse) = PolynomialSparse(-collect(p.terms))

#-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms), map((pt)->-pt, values(l.dict)))



"""
Create a new polynomial which is the derivative of the polynomial.
"""
function derivative(p::PolynomialSparse)::PolynomialSparse
    der_p = PolynomialSparse(Term(0,0)) # zero polynomialsparse, has list and dict
    delete_element!(der_p.terms, der_p.dict, 0)
    for term in p.terms# will go from lowest to highest
        der_term = derivative(term)
        iszero(der_term) ? nothing : insert_sorted!(der_p.terms, der_p.dict, der_term.degree, der_term)
        # # if constant, derivative will be zero term and will have two terms of zero degree that conflict. will only potentially need to worry about this for first term, lowest
        # # delete the original zero term in der_p. can't have empty polynomial sparse, so insert something else as placeholder first THEN insert the new constant term
        #     delete_element!(der_p.terms, der_p.dict, 0)
        #     insert_sorted!(der_p.terms, der_p.dict, 0, der_term)
        #     delete_element!(der_p.terms, der_p.dict, 0)
        # # account for when constants derive into zero - insert_sorted doesnt like having this zero term, since der_p is initialised to only have the zero term
        # else
        #     insert_sorted!(der_p.terms, der_p.dict, der_term.degree, der_term)
        # end
        #use insert_sorted to add the key and value of the derived termdto the list and dict of der_p
    end
    return der_p
end

# """
# The prim part (multiply a polynomial by the inverse of its content).
# """
# prim_part(p::PolynomialSparse) = p ÷ content(p)


# """
# A square free polynomial.
# """
# square_free(p::Polynomial, prime::Int)::Polynomial = (p ÷ gcd(p,derivative(p),prime))(prime)

# #################################
# # Queries about two polynomials #
# #################################

"""
Check if two polynomials are the same
"""
==(p1::PolynomialSparse, p2::PolynomialSparse)::Bool = p1.terms == p2.terms


# """
# Check if a polynomial is equal to 0.
# """
# #Note that in principle there is a problem here. E.g The polynomial 3 will return true to equalling the integer 2.
# ==(p::Polynomial, n::T) where T <: Real = iszero(p) == iszero(n)

# ##################################################################
# # Operations with two objects where at least one is a polynomial #
# ##################################################################

# """
# Subtraction of two polynomials.
# """
# -(p1::Polynomial, p2::Polynomial)::Polynomial = p1 + (-p2)


# """
# Multiplication of polynomial and term.
# """
# *(t::Term, p1::Polynomial)::Polynomial = iszero(t) ? Polynomial() : Polynomial(map((pt)->t*pt, p1.terms))
# *(p1::Polynomial, t::Term)::Polynomial = t*p1

# """
# Multiplication of polynomial and an integer.
# """
# *(n::Int, p::Polynomial)::Polynomial = p*Term(n,0)
# *(p::Polynomial, n::Int)::Polynomial = n*p

# """
# Integer division of a polynomial by an integer.

# Warning this may not make sense if n does not divide all the coefficients of p.
# """
# ÷(p::Polynomial, n::Int) = (prime)->Polynomial(map((pt)->((pt ÷ n)(prime)), p.terms))

# """
# Take the mod of a polynomial with an integer.
# """
# function mod(f::Polynomial, p::Int)::Polynomial
#     f_out = deepcopy(f)
#     for i in 1:length(f_out.terms)
#         f_out.terms[i] = mod(f_out.terms[i], p)
#     end
#     return trim!(f_out)
        
#     # p_out = Polynomial()
#     # for t in f
#     #     new_term = mod(t, p)
#     #     @show new_term
#     #     push!(p_out, new_term)
#     # end
#     # return p_out
# end

# """
# Power of a polynomial mod prime.
# """
# function pow_mod(p::Polynomial, n::Int, prime::Int)
#     n < 0 && error("No negative power")
#     out = one(p)
#     for _ in 1:n
#         out *= p
#         out = mod(out, prime)
#     end
#     return out
# end