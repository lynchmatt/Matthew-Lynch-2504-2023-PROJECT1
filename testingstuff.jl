using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### TESTING #####
f = Term(11,4)
g = Term(3,2)
s1 = PolynomialSparse([Term(5,3), Term(2,9), Term(3,4)])
s2 = PolynomialSparse([Term(7,3), Term(2,2), Term(8,1)])
s3  = PolynomialSparse([Term(3,2), Term(9,1), Term(3,0)])
d1 = PolynomialDense([Term(4,3), Term(2,9), Term(3,4)])

factor(s1,5)
factor(s1modp)
x = x_polysparse()

square_free(s2,5)
square_free(s2modp)


s1modp = PolynomialModP(s1,5)
s2modp = PolynomialModP(s2, 5)
s3modp = PolynomialModP(s3,11)
divide(s2modp, s1modp)
rem(s2modp, s1modp)
÷(s2modp, s1modp)
divide(s2modp, s3modp)
rem(s2modp, s3modp)

check = 2x^9 + 2x^3
mod(check^2,3)
s1modp^2
pow_mod(s1modp,2)

4 - s1modp

length(s1)
length(modp)
length(s2)
length(s2modp)
leading(s1modp)
leading(s2modp)
coeffs(s1modp)
coeffs(s2modp)
degree(s1)
degree(s1modp)
degree(s2modp)
evaluate(s1, 2)
evaluate(s2modp,2)
evaluate(check,2)


s1modp == s1modp2

iszero(s2modp)


m = PolynomialModP(s1, 3)
z = zero(PolynomialModP,3)
iszero(z)
a = one(PolynomialModP, 3)
one(m)
a.polynomial
m.polynomial
zero(PolynomialModP,3)
rand(Polynomial)

cyclotonic_polynomialmodp(2,3)

linear_monic_polynomialmodp(3,5)

x = x_polymodp(3)

mod(s1,3)

### SORTING OUT VALUE ACCESS ###

x = x_polymodp(3)
k = x^2 + 3x + 4x^6
for t in k.polynomial
    println(t)
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

