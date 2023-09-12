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
prod_test_polydense()
prod_derivative_test_polydense()
ext_euclid_test_polydense()
division_test_polydense()

####
# Execute unit tests for polynomialsparse
####
include("polynomialsparse_test.jl")


####
# Execute unit tests for polynomial factorization
####
include("factorization_test.jl")
factor_test_poly()