using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### SPARSE TESTING #####

l = PolynomialSparse([Term(9,2), Term(2,4), Term(-3, 3), Term(4,0)])
p = derivative(l) # THIS FUCKING WORKS

j = PolynomialSparse([Term(9,2), Term(2,4)])
delete_element!(j.terms, j.dict, 4)
j
j.terms

der_p = PolynomialSparse(Term(0,0)) # zero polynomialsparse, has list and dict
delete_element!(der_p.terms, der_p.dict, 0)


haskey(l.dict, 0)
derivative(l)
# if haskey zero to begin with, delete the first term at the very end, since they'll be a double-up?
insert_sorted!(l.terms, l.dict, 0, Term(9,0))
l.terms
derivative(l)
# will auto-convert Term(0,Int) to be the zero term, so it's now a constant and a constant already exists
a = Term(0,2)
a.degre
leading(l)
coeffs(l)
degree(l)
content(l)
evaluate(l, 0)
iszero(l)
trim!(l)
p = derivative(l)
p.terms
trim!(p)
iszero(p.terms[1])

p.terms

-collect(l.terms)


example = PolynomialSparse([Term(2,4), Term(6,5), Term(8,3), Term(2,0)])
-(p::PolynomialSparse) = PolynomialSparse(map((pt)->-pt, p.terms), map((pt)->-pt, values(p.dict)))
-example

values(l.dict)


### SORTING OUT VALUE ACCESS ###
for t in l.terms
    println(-t)
end

map((pt)->(-pt), values(l.dict))

j = map((pt)->-pt, values(l.dict))
k = get_element(l.terms, l.dict, 3) 
k.degree

i = PolynomialSparse()
i.dict
i.terms
i.terms == MutableLinkedList{Term}(zero(Term))
iszero(i)
evaluate(i, 2)
degree(i)
content(i)


##### DENSE TESTING #####


p = Polynomial([Term(3,4), Term(4,5), Term(2,3)])
p.terms
iterate(p, 1)
leading(p)
derivative(p)


for t in l.terms # ITERATE IN POLYNOMIALSPARSE.JL WORKS
    println(t)
end


q = MutableLinkedList{Int}(1,2,3)
for i in 1:q.len
    println(getindex(q,i))
end

sum(q)

l.terms[2] # access specific term in polynomial l

l.terms[2].coeff # access the coefficeint of a specific term of polynomial l1

l.dict[3] # access the term of degree 3 of polynomial l
l.dict[0]

a = get_element(l.terms, l.dict, 3)
a.degree # access the term of degree 3 of polynomial 1, but nicer

k = PolynomialSparse(Term(2,3)) # single, non-vector term works




function alphabetprinter()
    alphabet = ["a", "b", "c", "d"]
    # make local variable
    if (@isdefined out_of_order) == false
        localvar = false
    elseif out_of_order == false 
        localvar = false
    else
        localvar = true
    end
    for (i,t) in (localvar ? enumerate(reverse(alphabet)) : enumerate(alphabet)) 
        println(t)
    end
end


function +(p1::PolynomialSparse, p2::PolynomialSparse)::PolynomialSparse
    a = deepcopy(p1)
    b = deepcopy(p2)
    for (i,t) in enumerate(a.dict)
        println(t[1])
    end

end

