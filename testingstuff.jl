using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### TESTING #####


s = PolynomialSparse([Term(1,8), Term(2,4)])
q = PolynomialSparse([Term(2,1), Term(1,2)])

d = PolynomialDense([Term(2,8), Term(2,4)])
e = PolynomialDense([Term(2,2), Term(1,4)])


mod(1,3)

x = x_polysparse()

x + 2 + 4

mod(1 - (2x+2)*0, 3)
a = (4 + 2x)
typeof(a)
+(a::PolynomialSparse,2::Int)

g,s,t = 1, 0, x + 1
a, b = 2x^2 + x, 4x + 1

mod((s*a + t*b) - g, 3)




d = PolynomialDense([Term(1,1), Term(2,2)])
e = PolynomialDense([Term(2,0), Term(1,1)])
square_free(d,3)
multiplicity(d,e,2)
square_free(d,5)

function multiplicity(f::PolynomialSparse, g::PolynomialSparse, prime::Int)::Int
    degree(gcd(f, g, prime)) == 0 && return 0
    return 1 + multiplicity((÷(f,g)(prime)), g, prime)
end

# issue is with the assertion error in gcd - refuses to acknowledge that it's zero. Think the cause is how i've defined the zero polynomialsparse. 
# gcd works for polysparse if i remove the assertion error.
# all of square free will work if i remove the assertion error - need to fix the zero thing probably

gcd(a,derivative(a),3)

square_sparse(s,2)
square_free(a,3)


a = PolynomialDense([Term(1,3), Term(8,7)])
b = PolynomialDense([Term(2,1), Term(1,2)])
square_free(a, 3)
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


##### DENSE TESTING #####

e = PolynomialDense()
x = x_polydense()
z = x^2 + 2x
e.terms

p = PolynomialDense([Term(3,4), Term(4,5), Term(2,3)])
r = PolynomialDense([Term(4,7), Term(9,1)])


q = MutableLinkedList{Int}(1,2,3)
for i in 1:q.len
    println(getindex(q,i))
end

sum(q)


a = get_element(l.terms, l.dict, 3)
a.degree # access the term of degree 3 of polynomial 1, but nicer




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

