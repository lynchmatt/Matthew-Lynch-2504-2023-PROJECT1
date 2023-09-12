using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### TESTING #####

densetest = PolynomialDense([Term(7,3), Term(2,2), Term(8,1), Term(1,0)])
sparsetest = PolynomialSparse([Term(7,3), Term(2,2), Term(8,1), Term(1,0)])

prim_part(sparsetest)(7)

prim_part(densetest)(7)

parts = factor(densetest, 11)
sparseparts = factor(sparsetest, 11)

expand_factorization(parts)
f = expand_factorization(sparseparts)

### SORTING OUT VALUE ACCESS ###

for t in l.terms
    println(-t)
end


##### DENSE TESTING #####

p = PolynomialDense([Term(3,4), Term(4,5), Term(2,3)])
r = PolynomialDense([Term(4,7), Term(9,1)])


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

