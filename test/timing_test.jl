################################################################################
################################################################################
#
# A script that compares the runtime of PolynomialSparse and PolynomialSparse128
#                                                                               
################################################################################
################################################################################

# create polysparse
a64 = PolynomialSparse([Term(9223372036854775807, 2), Term(9223372036854775807, 3)])
b64 = PolynomialSparse([Term(9223372036854775807, 2), Term(9223372036854775807, 1)])
println("Time the multiplication of the Int64 PolynomialSparses $a64 and $b64")
@time a64*b64

# create polysparse128
a128 = PolynomialSparse128([Term128(9223372036854775807, 2), Term128(9223372036854775807, 3)])
b128 = PolynomialSparse128([Term128(9223372036854775807, 2), Term128(9223372036854775807, 1)])
println("Time the multiplication of the Int128 PolynomialSparses $a128 and $b128")
c = @time a128*b128