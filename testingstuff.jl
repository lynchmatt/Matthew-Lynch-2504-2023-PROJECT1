using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

l = PolynomialSparse([Term(2,3), Term(3,4)])

l.terms[2] # access specific term in polynomial l

l.terms[2].coeff # access the coefficeint of a specific term of polynomial l1

l.dict[3] # access the term of degree 3 of polynomial l

get_element(l.terms, l.dict, 3) # access the term of degree 3 of polynomial 1, but nicer

k = PolynomialSparse([Term(2,3)])

x = x_poly()
typeof(x)

y = one(Polynomial)

h = Term(1,1)


show(Term(2,3)) # will display constants correctly when the coefficient is non zero, ie without x^0

y1 = -1x^1 + 2x^3 -3x^0 -4x^2  #will display without +/-1, without powers of zero, and without x^1 when put into polynomial, but not when using show?
y2 = 1x^0 + -3x^2 # need to print 1 as constant, currently won't
show(y1)

using DataStructures
l = MutableLinkedList{Int}(1,2,3)
k = PolynomialSparse([Term(2,3), Term(4,5), Term(1,1)])
for (i,t) in enumerate(k)
    println(t.coeff)
end

dict = Dict(1=> "a", 2=>"b", 4=> "d")

for (i,t) in enumerate(dict)
    println(t[2]) #use t[1] to get key, t[2] to get value. using only t will give each key-value pair
end

insert_sorted!(l, dict, 2, 2)





import Pkg; Pkg.add("Subscripts")
using Subscripts
super("1")


# bring up the base.show issue?


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

k = PolynomialSparse([Term(2,3)])
l = PolynomialSparse([Term(2,3), Term(3,4)])

+(l,k) 
