#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

##################
# POLY PLUS TERM #
##################

"""
Add a polynomialdense and a term.
"""
function +(p::PolynomialDense, t::Term)
    p = deepcopy(p)
    if t.degree > degree(p)
        push!(p, t)
    else
        if !iszero(p.terms[t.degree + 1]) #+1 is due to indexing
            p.terms[t.degree + 1] += t
        else
            p.terms[t.degree + 1] = t
        end
    end

    return trim!(p)
end
+(t::Term, p::PolynomialDense) = p + t


"""
Add a polynomialsparse and a term.
"""
function +(p::PolynomialSparse, t::Term)
    # include edge cases for zero term and zero polynomial
    if iszero(t) == true
        return p
    else
        nothing
    end
    p = deepcopy(p)
    checkelement = get_element(p.terms, p.dict, t.degree) # try get the element of the same degree as the term we're adding
    if isnothing(checkelement) # if doesnt have a term of that degree
        insert_sorted!(p.terms, p.dict, t.degree, t)
    else #case where we're adding the term to an existing term
        checkelement += t #term addition
        # if they cancel out, the coefficient being zero auto sets the degree to be zero as well. if this happens, remove that degree but don't add anything
        delete_element!(p.terms, p.dict, t.degree)
        if checkelement.coeff == 0
            nothing
        else
            insert_sorted!(p.terms, p.dict, checkelement.degree, checkelement)
        end
    end

    return p
end

+(t::Term, p::PolynomialSparse) = p + t


"""
Add a polynomialsparse128 and a term.
"""
function +(p::PolynomialSparse128, t::Term128)
    if iszero(t) == true
        return p
    else
        nothing
    end
    p = deepcopy(p)
    checkelement = get_element(p.terms, p.dict, t.degree) # try get the element of the same degree as the term we're adding
    if isnothing(checkelement) # if doesnt have a term of that degree
        insert_sorted!(p.terms, p.dict, t.degree, t)
    else #case where we're adding the term to an existing term
        checkelement += t #term addition
        # if they cancel out, the coefficient being zero auto sets the degree to be zero as well. if this happens, remove that degree but don't add anything
        delete_element!(p.terms, p.dict, t.degree)
        if checkelement.coeff == 0
            nothing
        else
            insert_sorted!(p.terms, p.dict, checkelement.degree, checkelement)
        end
    end
    return p
end

+(t::Term128, p::PolynomialSparse128) = p + t

"""
Methods for adding terms to PolynomialModP
"""
+(p::PolynomialModP, t::Term) = PolynomialModP((p.polynomial+mod(t,p.prime)), p.prime)
+(t::Term, p::PolynomialModP) = p + t

"""
Methods for adding terms to PolynomialModP128
"""
+(p::PolynomialModP128, t::Term128) = PolynomialModP128((p.polynomial+mod(t,p.prime)), p.prime)
+(t::Term128, p::PolynomialModP128) = p + t

##################
# POLY PLUS POLY #
##################

"""
Add two polynomialdenses.
"""
function +(p1::PolynomialDense, p2::PolynomialDense)::PolynomialDense
    p = deepcopy(p1)
    for t in p2.terms
        p += t
    end
    return p
end


"""
Add two polynomialsparses.
"""
function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    a = deepcopy(p1)
    for t in p2.terms
        t
        a += t
    end
    return a
end

"""
Add two polynomialsparse128s.
"""
function +(p1::PolynomialSparse128, p2::PolynomialSparse128)::PolynomialSparse128
    a = deepcopy(p1)
    for t in p2.terms
        t
        a += t
    end
    return a
end

"""
Add two polynomialmodps, provided they are mod the same prime
"""
function +(p1::PolynomialModP, p2::PolynomialModP)::PolynomialModP
    @assert p1.prime == p2.prime
    return PolynomialModP((p1.polynomial+p2.polynomial), p1.prime)
end

"""
Add two polynomialmodp128s, provided they are mod the same prime
"""
function +(p1::PolynomialModP128, p2::PolynomialModP128)::PolynomialModP128
    @assert p1.prime == p2.prime
    return PolynomialModP128((p1.polynomial+p2.polynomial), p1.prime)
end


##################
# POLY PLUS INT  #
##################

"""
Add a polynomialdense and an integer.
"""
+(p::PolynomialDense, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialDense) = p + Term(n,0)
"""
Add a polynomialsparse and an integer.
"""
+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)

"""
Add a polynomialsparse128 and an integer.
"""
+(p::PolynomialSparse128, n::Integer) = p + Term128(n,0)
+(n::Integer, p::PolynomialSparse128) = p + Term128(n,0)

"""
Add a polynomialmodp and an integer.
"""
+(p::PolynomialModP, n::Integer) = p + Term(n,0)
+(n::Integer, p::PolynomialModP) = p + Term(n,0)

"""
Add a polynomialmodp128 and an integer.
"""
+(p::PolynomialModP128, n::Integer) = p + Term128(n,0)
+(n::Integer, p::PolynomialModP128) = p + Term128(n,0)

##################
# POLY MINUS INT #
################## 

"""
Subtraction of an integer from a polynomialdense
"""
-(p1::PolynomialDense, n::Int)::PolynomialDense = p1 + -n
-(n::Int, p1::PolynomialDense)::PolynomialDense = -p1 + n


"""
Subtraction of an integer from a polynomialsparse
"""
-(p1::PolynomialSparse, n::Int)::PolynomialSparse = p1 + -n
-(n::Int, p1::PolynomialSparse)::PolynomialSparse = -p1 + n

"""
Subtraction of an integer from a polynomialsparse128
"""
-(p1::PolynomialSparse128, n::Integer)::PolynomialSparse128 = p1 + -n
-(n::Integer, p1::PolynomialSparse128)::PolynomialSparse128 = -p1 + n

"""
Subtraction of an integer from a polynomialmodp
"""
-(p1::PolynomialModP, n::Integer)::PolynomialModP = p1 + -n
-(n::Integer, p1::PolynomialModP)::PolynomialModP = -p1 + n

"""
Subtraction of an integer from a polynomialmodp
"""
-(p1::PolynomialModP128, n::Integer)::PolynomialModP128 = p1 + -n
-(n::Integer, p1::PolynomialModP128)::PolynomialModP128 = -p1 + n