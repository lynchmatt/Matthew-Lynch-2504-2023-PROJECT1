using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### TESTING #####

f = Term128(9,4)
g = Term128(3,2)
s1 = PolynomialSparse128([Term128(4,3), Term128(2,9), Term128(3,4)])
s2 = PolynomialSparse128([Term128(7,3), Term128(2,2), Term128(8,1)])
d1 = PolynomialDense([Term(4,3), Term(2,9), Term(3,4)])
d2 = PolynomialDense([Term(7,3), Term(2,2), Term(8,1)])
easy = PolynomialDense([Term(3,2), Term(2,1), Term(3,0)])


sparse = PolynomialSparse([Term(100,2), Term(100,1), Term(100,0)])
sparse128 = PolynomialSparse128([Term128(100,2), Term128(100,1), Term128(100,0)])

iszero(mod(sparse128, 10))
iszero(mod(sparse128^10, 10))

function overflow_test()
    sparsepoly = PolynomialSparse([Term(100,2), Term(100,1), Term(100,0)])
    sparsepoly128 = PolynomialSparse128([Term128(100,2), Term128(100,1), Term128(100,0)])
    # the coefficients of both should be a multiple of 10, when raised to any power. check if they are equal
    if sparsepoly^10 == sparsepoly128^10
        println("No overflow in PolynomialSparse")
        return
    end
    # check that the coefficients of the polynomialsparse128 are in fact multiples of ten, when raised to the tenth power
    if iszero(mod(sparsepoly128^10, 10))
    else
        throw(OverflowError("PolynomialSparse128 has overflown"))
    end
    # check if polynomialsparse has overflown
    if iszero(mod(sparsepoly^10, 10))
        println("PolynomialSparse has not overflown.")
    else
        println("PolynomialSparse has overflown.")
    end
end


maxpolysparse = PolynomialSparse([Term(92233720368547758071,1)])
maxpolysparse += Term128(1,1)

overflow_test()

overflow = sparse^10
typeof(overflow)

dontoverflow = sparse128^10



@show p1 = rand(PolynomialSparse128)
length(p1)
@show p1plus = p1 + Term128(0,0)

p1plus = p1 + Term128(0,0)
length(p1)

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

