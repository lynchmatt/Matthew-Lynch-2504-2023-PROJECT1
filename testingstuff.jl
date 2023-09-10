using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### SPARSE TESTING #####
s = PolynomialSparse([Term(1,3), Term(8,7)])
q = PolynomialSparse([Term(2,1), Term(1,2)])
t = extended_euclid_alg(s,q,3)
typeof(t[2])
square_free(s, 2)

square_free(p::PolynomialSparse, prime::Int)::PolynomialSparse = (÷(p, gcd(p,derivative(p),prime)))(prime)

a = ÷(s, gcd(s,derivative(s),2))
a(2)


a = Polynomial([Term(1,3), Term(8,7)])
b = Polynomial([Term(2,1), Term(1,2)])
square_free(a, 1)
extended_euclid_alg(a, b, 3)
gcd(a,b, 3)
# returns the greatest common denominator, which can be expressed as a linear combination of a and b - the coordinates of the linear combo are the other numbers
# ie returns
÷(s,2)

empty = PolynomialSparse(zero(Term))
delete_element!(empty.terms, empty.dict, 0)


x = x_polysparse()

y1 = x^2 + 3x

y1 == a

t = Term(2,3)

a = PolynomialSparse([Term(6,4), Term(6,9)])

a^2

p1 = PolynomialSparse([Term(4,6), Term(8,3)])

p1*a


k = PolynomialSparse()
iszero(k)
delete_element!(k.terms, k.dict, 0)

j = PolynomialSparse(Term(0,0))
iszero(j)
delete_element!(j.terms, j.dict, 0)



p = derivative(l) # THIS FUCKING WORKS

j = PolynomialSparse([Term(9,2), Term(2,4)])
delete_element!(j.terms, j.dict, 4)
j
j.terms



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

e = Polynomial()
e.terms

p = Polynomial([Term(3,4), Term(4,5), Term(2,3)])
e + p

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

