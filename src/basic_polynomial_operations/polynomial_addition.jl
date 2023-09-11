#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

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
    p = deepcopy(p)
    tmp_element = get_element(p.terms, p.dict, t.degree) # try get the element of the same degree as the term we're adding
    if isnothing(tmp_element) # if doesnt have a term of that degree
        insert_sorted!(p.terms, p.dict, t.degree, t)
    else #case where we're adding the term to an existing term
        tmp_element += t #term addition
        # if they cancel out, the coefficient being zero auto sets the degree to be zero as well. if this happens, remove that degree but don't add anything
        delete_element!(p.terms, p.dict, t.degree)
        if tmp_element.coeff == 0
            nothing
        else
            insert_sorted!(p.terms, p.dict, tmp_element.degree, tmp_element)
        end
    end

    return p
end

+(t::Term, p::PolynomialSparse) = p + t

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
Add a polynomial and an integer.
"""
+(p::PolynomialDense, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialDense) = p + Term(n,0)

+(p::PolynomialSparse, n::Int) = p + Term(n,0)
+(n::Int, p::PolynomialSparse) = p + Term(n,0)

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