using Pkg
Pkg.activate(".")

#COMMENT OUT ALL ABOVE

include("poly_factorization_project.jl")

##### TESTING #####
s1 = PolynomialSparse([Term(3,1), Term(4,0)])
s2 = PolynomialSparse([Term(6,1), Term(5,0)])
s3  = PolynomialSparse([Term(3,2), Term(9,1), Term(3,0)])
d1 = PolynomialDense([Term(4,3), Term(2,9), Term(3,4)])
s1modp = PolynomialModP(s1,5)
s2modp = PolynomialModP(s2, 5)
s3modp = PolynomialModP(s3,11)

function CRT_B(a::PolynomialSparse, b::PolynomialSparse, n::Integer, m::Integer)::PolynomialSparse
    x = x_polysparse()
    c = PolynomialSparse()
    while !isempty(a.terms) || !isempty(b.terms)
        k = max(degree(a), degree(b))
        if k > degree(a)
            ak = 0
        else
            ak = leading(a).coeff
            delete_element!(a.terms, a.dict, leading(a).degree)
        end
        if k > degree(b)
            bk = 0
        else
            bk = leading(b).coeff
            delete_element!(b.terms, b.dict, leading(b).degree)
        end
        @show ck = CRT_int([Int(ak), Int(bk)], [n, m])
        c = c + ck*x^k
    end
    return c
end

CRT_L(s1,s2,5,7)
CRT_B(s1, s2, 5, 7)

a = 1
Bool(0)
Int(trunc(log2(129)))

maxpower = Int(trunc(log2(64)))
binary = [2^i for i in 0:maxpower]
binary_string = digits(2, base=2, pad=length(binary))


function pow_mod(p::PolynomialModP, n::Int)
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

