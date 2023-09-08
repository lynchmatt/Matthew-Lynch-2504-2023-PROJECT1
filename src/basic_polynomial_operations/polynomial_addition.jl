#############################################################################
#############################################################################
#
# This file implements polynomial addition 
#                                                                               
#############################################################################
#############################################################################

"""
Add a polynomial and a term.
"""
function +(p::Polynomial, t::Term)
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
+(t::Term, p::Polynomial) = p + t


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
        delete_element!(p.terms, p.dict, t.degree)
        insert_sorted!(p.terms, p.dict, tmp_element.degree, tmp_element)
    end

    return p
end

+(t::Term, p::Polynomial) = p + t

"""
Add two polynomialsparse.
"""
function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    a = deepcopy(p1)
    b = deepcopy(p2)
    for (i,t) in enumerate(a.dict[1])
    # need to check each key of a and compare with b. If b has a key of the same degree, add their coefficients together and make a key of that degree - if their added coefficients equal
    # zero, don't make that degree into a key. If one polynomial has a key that the other doesn't, preserve the original term and insert it into the linkedlist.
    # need to check both directions?

    for t in p2
        p += t
    end
    return p
end


"""
Add two polynomials.
"""
function +(p1::Polynomial, p2::Polynomial)::Polynomial
    p = deepcopy(p1)
    for t in p2.terms
        p += t
    end
    return p
end

"""
Add a polynomial and an integer.
"""
+(p::Polynomial, n::Int) = p + Term(n,0)
+(n::Int, p::Polynomial) = p + Term(n,0)