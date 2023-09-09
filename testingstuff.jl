using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### SPARSE TESTING #####


x = x_polysparse()

y1 = x^2 + 3x

y1 == a

t = Term(2,3)

a = PolynomialSparse([Term(6,4), Term(6,9)])

a^2

p1 = PolynomialSparse([Term(4,6), Term(8,3)])

p1*a


    # der_p = PolynomialSparse(Term(0,0)) # zero polynomialsparse, has list and dict
    # delete_element!(der_p.terms, der_p.dict, 0)
    # for term in p.terms# will go from lowest to highest
    #     der_term = derivative(term)
    #     iszero(der_term) ? nothing : insert_sorted!(der_p.terms, der_p.dict, der_term.degree, der_term)
    # end
    # return der_p

# emptypoly = PolynomialSparse(Term(0,0))
# delete_element!(emptypoly.terms, emptypoly.dict, 0)
# for vt in n.terms
#     a = vt*t
#     iszero(a) ? nothing : insert_sorted!(emptypoly.terms, emptypoly.dict, a.degree, a)
# end

emptypoly

PolynomialSparse(init)

m.terms.*t

k = PolynomialSparse()
iszero(k)
delete_element!(k.terms, k.dict, 0)











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
func = p ÷ 4
func(5)
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

