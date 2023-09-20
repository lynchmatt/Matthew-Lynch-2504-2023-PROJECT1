#################################################################################################
#################################################################################################
#
# This file implements the CRT for integers, polynomials and polynomial_sparse_128 multiplication
#                                                                               
#################################################################################################
#################################################################################################


###########################
# SYMMETRIC MOD FUNCTIONS #
###########################

"""
Create a method for the symmetric mod for integers
"""
function smod(a::Integer, m::Integer)
    if mod(a,m) <= m//2
        return mod(a,m)
    else
        return mod(a,m) - m
    end
end

"""
Compute the symmetric mod of a Term128 with an integer.
"""
smod(t::Term128, p::Integer) = Term128(Int128(smod(t.coeff,p)), Int128(t.degree))

"""
Create a method for the symmetric mod for a polynomialsparse128 and an integer
"""
function smod(f::PolynomialSparse128, p::Integer)::PolynomialSparse128
    smod_vector = Vector{Term128}(undef, length(f.terms))
    for (i,t) in enumerate(f.terms)
        smod_vector[i] = smod(t, p)
    end
    return PolynomialSparse128(smod_vector)
end

#################
# CRT FUNCTIONS #
#################

"""
Implments the CRT using two 2-element vectors of integers.
"""
function CRT_int(ui::Vector{T}, m::Vector{T}) where T <: Union{Int64, Int128}
    ui, m = Integer.(ui), Integer.(m)
    @assert length(ui) == length(m)
    v = Vector{Integer}(undef, length(ui))
    v[1] = ui[1]
    v[2] = (ui[2] -v[1])*int_inverse_mod(m[1], m[2]) % m[2]
    u = v[1] + v[2]*m[1]
    return u
end


"""
Implements the CRT for two polynomialsparse128s and two primes, using the CRT_int function.
"""
function CRT_poly(p1::PolynomialSparse128, p2::PolynomialSparse128, n::Integer, m::Integer)::PolynomialSparse128
    n, m = Int128(n), Int128.(m)
    a = deepcopy(p1)
    b = deepcopy(p2)
    x = x_polysparse128()
    c = PolynomialSparse128([Term128(0,0)])
    delete_element!(c.terms, c.dict, 0)
    while !iszero(leading(a)) || !iszero(leading(b))
        k = max(degree(a), degree(b))
        if k > degree(a)
            ak = 0
        else
            ak = leading(a).coeff
            if iszero(leading(a))
                nothing
            else
                delete_element!(a.terms, a.dict, leading(a).degree)
            end
        end
        if k > degree(b)
            bk = 0
        else
            bk = leading(b).coeff
            if iszero(leading(b))
                nothing
            else
                delete_element!(b.terms, b.dict, leading(b).degree)
            end
        end
        ck = CRT_int([Int128(ak), Int128(bk)], [Int128(n), Int128(m)])
        c = c + ck*x^k
    end
    return c
end

"""
Implements the CRT for two polynomialmodp128s (with their primes inbuilt), using the CRT_int function.
"""
function CRT_poly(p1::PolynomialModP128, p2::PolynomialModP128)::PolynomialModP128
    a = deepcopy(p1).polynomial
    b = deepcopy(p2).polynomial
    n = p1.prime
    m = p2.prime
    x = x_polysparse128()
    c = PolynomialSparse128([Term128(0,0)])
    delete_element!(c.terms, c.dict, 0)
    while !iszero(leading(a)) || !iszero(leading(b))
        k = max(degree(a), degree(b))
        if k > degree(a)
            ak = 0
        else
            ak = leading(a).coeff
            if iszero(leading(a))
                nothing
            else
                leading(a).degree
                delete_element!(a.terms, a.dict, leading(a).degree)
            end
        end
        if k > degree(b)
            bk = 0
        else
            bk = leading(b).coeff
            if iszero(leading(b))
                nothing
            else
                delete_element!(b.terms, b.dict, leading(b).degree)
            end
        end
        ck = CRT_int([Integer(ak), Integer(bk)], [n, m])
        c = c + ck*x^k
    end
    return PolynomialModP(c, n*m)
end

######################
# CRT MULTIPLICATION #
######################

"""
Polynomial128 Multiplication using the CRT
"""
function multiplication(a::PolynomialSparse128, b::PolynomialSparse128)
    height_a = maximum(abs.(coeffs(a)))
    height_b = maximum(abs.(coeffs(b)))
    B = Int128(2*height_a*height_b*min(degree(a)+1, degree(b)+1))
    p = Int128(3)
    M = Int128(p)
    moda = PolynomialModP128(a,M)
    modb = PolynomialModP128(b,M)
    c = (moda * modb).polynomial
    while M < B
        p = nextprime(p,2)
        moda = PolynomialModP128(a,p)
        modb = PolynomialModP128(b,p)
        c_prime = moda*modb
        c = CRT_poly(c, c_prime.polynomial, Int128(M), Int128(p))
        M *= p
    end
    return smod(c,M)
end