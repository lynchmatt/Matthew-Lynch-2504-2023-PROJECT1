"""
Factors a polynomialsparse128 over the field Z_p, using CRT and repeated squares functionality.

Returns a vector of tuples of (irreducible polynomials (mod p), multiplicity) such that their product of the list (mod p) is f. Irreducibles are fixed points on the function factor.
"""
function CRT_factor(f::PolynomialSparse128, prime::Integer)::Vector{Tuple{PolynomialSparse128,Integer}}
    #Cantor Zassenhaus factorization

    f_modp = mod(f, prime)
    degree(f_modp) ≤ 1 && return [(f_modp,1)]

    # make f primitive
    ff = prim_part(f_modp)(prime)      
    #@show "after prim:", ff

    # make f square-free
    squares_poly = gcd(f, derivative(ff), prime) 
    ff = (÷(ff, squares_poly))(prime) 
    # @show "after square free:", ff

    # make f monic
    old_coeff = leading(ff).coeff
    ff = (÷(ff,old_coeff))(prime)
    # @show "after monic:", ff

    dds = CRT_dd_factor(ff, prime)

    ret_val = Tuple{PolynomialSparse128,Integer}[]

    for (k,dd) in enumerate(dds)
        sp = CRT_dd_split(dd, k, prime)
        sp = map((p)->(p ÷ leading(p).coeff)(prime),sp) #makes the polynomials inside the list sp, monic
        for mp in sp
            push!(ret_val, (mp, CRT_multiplicity(f_modp,mp,prime)) )
        end
    end

    #Append the leading coefficient as well
    push!(ret_val, (multiplication(PolynomialSparse128(Term128(leading(f_modp).coeff, 0)), one(PolynomialSparse128)), 1) )

    return ret_val
end

"""
first(factorization[1])^last(factorization[1])
Expand a factorization for polynomialsparse128 using CRT
"""
function CRT_expand_factorization(factorization::Vector{Tuple{PolynomialSparse128,Integer}})::PolynomialSparse128
    length(factorization) == 1 && return ^(first(factorization[1]), last(factorization[1]))
    product = 1
    for tt in factorization
        pow = repsq_power(first(tt), last(tt))
        product = product*pow
    end
    return product
end

"""
Compute the number of times g divides f, for polynomialsparse128 using CRT
"""
function CRT_multiplicity(f::PolynomialSparse128, g::PolynomialSparse128, prime::Integer)::Integer
    degree(gcd(f, g, prime)) == 0 && return 0
    return 1 + CRT_multiplicity((÷(f,g)(prime)), g, prime)
end

"""
Distinct degree factorization for polynomialsparse128

Given a square free polynomial `f` returns a list, `g` such that `g[k]` is a product of irreducible polynomials of degree `k` for `k` in 1,...,degree(f) ÷ 2, such that the product of the list (mod `prime`) is equal to `f` (mod `prime`).
"""

function CRT_dd_factor(f::PolynomialSparse128, prime::Integer)::Array{PolynomialSparse128}
    x = x_polysparse128()
    w = deepcopy(x)
    g = Array{PolynomialSparse128}(undef,degree(f)) #Array of polynomials indexed by degree

    #Looping over degrees
    for k in 1:degree(f)
        w = rem(repsq_pow_mod(w,prime,prime), f)(prime) # remainder of x to the power of prime mod prime divided by f, mod prime
        g[k] = gcd(w - x, f, prime) # kth item in the polynomial array is the gcd of w minus unit polynomial with f, mod prime
        f = (÷(f,g[k]))(prime)
    end

    #edge case for final factor
    f != one(PolynomialSparse128) && push!(g,f)
    
    return g
end

"""
Distinct degree split for polynomialsparse128

Returns a list of irreducible polynomials of degree `d` so that the product of that list (mod prime) is the polynomial `f`.
"""
function CRT_dd_split(f::PolynomialSparse128, d::Integer, prime::Integer)::Vector{PolynomialSparse128}
    f = mod(f,prime)
    degree(f) == d && return [f]
    degree(f) == 0 && return []
    w = rand(PolynomialSparse128, degree = d, monic = true)
    w = mod(w,prime)
    n_power = (prime^d-1) ÷ 2
    g = gcd(repsq_pow_mod(w,n_power,prime) - one(PolynomialSparse128), f, prime)
    ḡ = (÷(f,g))(prime) # g\bar + [TAB]
    return vcat(dd_split(g, d, prime), dd_split(ḡ, d, prime) )
end
