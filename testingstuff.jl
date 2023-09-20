using Pkg
Pkg.activate(".")


include("poly_factorization_project.jl")


##### TESTING #####
s1 = PolynomialSparse128([Term128(3,2), Term128(4,0)])
s2 = PolynomialSparse128([Term128(6,1), Term128(5,0)])
s3 = PolynomialSparse128([Term128(2,2)])
s4 = PolynomialSparse128([Term128(9223372036854775807,1), Term128(5^9,3), Term128(4,80)])
repsq_power(s2,2)
s2^2
@time s2*s4
@time multiplication(s2,s4)

t1 = PolynomialSparse([Term(3,2), Term(4,0)])
t2 = PolynomialSparse([Term(6,1), Term(5,0)])

CRT_int([0,4], [3,5])

s3  = PolynomialSparse([Term(3,2), Term(9,1), Term(3,0)])
d1 = PolynomialDense([Term(4,3), Term(2,9), Term(3,4)])
s1modp = PolynomialModP128(s1,5)
s2modp = PolynomialModP(s2, 7)
s3modp = PolynomialModP(s3,11)

"""
Implements the repeated squares method of powers for PolynomialModP
"""
function repsq_power(p::PolynomialModP, n::Integer)
    # find max binary power needed to reach exponent
    maxpower = Int(trunc(log2(n)))
    exponents = [2^i for i in 0:maxpower]
    binary_string = digits(n, base=2, pad=maxpower) # convert to binary string
    outpoly = one(Term)
    for i in 1:length(exponents)
        if binary_string[i] == 1
            outpoly *= ^(p, exponents[i])
        else
            nothing
        end
    end
    return outpoly
end




function l_pow(p::PolynomialModP, n::Int)
    max_power = floor(Int, log2(n))
    powers::Vector{typeof(p)} = [p] #p^1, p^2, p^4, p^8, p^16, ...
    for i in 1:max_power
        powers = push!(powers, powers[end]*powers[end])
    end
    out = one(p)
    for i in 0:max_power
        if n & (1 << i) != 0 #if n has a 1 in the ith bit
            out*=powers[i+1]
        end
    end
    return out
end



l_pow(s1modp, 2)
pow_mod_new(s1modp,2)

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

