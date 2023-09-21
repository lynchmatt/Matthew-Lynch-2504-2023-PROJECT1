#############################################################################
#############################################################################
#
# A script that runs all unit tests in the project.
#                                                                               
#############################################################################
#############################################################################

using Pkg
Pkg.activate(".")

include("../poly_factorization_project.jl")

####
# Execute unit tests for integers
###
include("integers_test.jl")
test_euclid_ints()
test_ext_euclid_ints()

####
# Execute unit tests for polynomialdense
####
include("polynomialdense_test.jl")


####
# Execute unit tests for polynomialsparse
####
include("polynomialsparse_test.jl")

####
# Execute unit tests for polynomialsparse128
####
include("polynomialsparse128_test.jl")

####
# Execute unit tests for polynomialmodp
####
include("polynomialmodp_test.jl")

####
# Execute unit tests for all types of polynomial factorization
####
include("factorization_test.jl")
factor_test_polydense()
factor_test_polysparse()
factor_test_polysparse128()
factor_test_polymodp()
CRT_factor_test()


####
# Execute unit tests for CRT multiplication
###
include("CRT_test.jl")

