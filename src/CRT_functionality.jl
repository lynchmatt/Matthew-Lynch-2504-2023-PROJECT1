#################################################################################################
#################################################################################################
#
# This file implements the CRT for integers, polynomials and polynomial_sparse_128 multiplication
#                                                                               
#################################################################################################
#################################################################################################

"""
Implments the CRT using two 2-element vectors of integers.
"""
function CRT_int(ui::Vector{Int}, m::Vector{Int})
    @assert length(ui) == length(m)
    v = Vector{Int}(undef, length(ui))
    v[1] = ui[1]
    v[2] = (ui[2] -v[1])*invmod(m[1], m[2]) % m[2]
    u = v[1] + v[2]*m[1]
    return u
end


"""
Implements the CRT for two polynomialmodps, using the CRT_int function. Have to have the same prime?
"""
function CRT_poly_L(p1::PolynomialSparse, p2::PolynomialSparse, m::Int, n::Int)
    v1 = p1
    v2 = p2 - v1
    v2 = v2*int_inverse_mod(m,n)
    v2 = mod(v2, n)
    p = v1 + v2*m
    return p
end




# 3, 4 and 0 - coefficients of the leading terms using CRT on each polynomial
# two polys, input ab and nm. 3mod5 and 4mod7. dont need to 
# two polys, a and b. a = x mod m and b= x mod n for some polynomial x, for some integers m and n
# consecutive pairs to get past the threshold. will produce c, where c = x mod mn







# for polymodps a and b, pass in a.polynomial, b.polynomial, a.prime, b.prime?

