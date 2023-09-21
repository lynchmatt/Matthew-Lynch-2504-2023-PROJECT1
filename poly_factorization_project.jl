#############################################################################
#############################################################################
#
# This is the main project file for polynomial factorization
#                                                                               
#############################################################################
#############################################################################

using Distributions, StatsBase, Random, DataStructures, Primes, BenchmarkTools

import Base: %
import Base: push!, pop!, iszero, show, isless, map, map!, iterate, length, last
import Base: +, -, *, mod, %, ÷, ==, ^, rand, rem, zero, one
import Primes: factor

global lowest_to_highest = false

include("src/general_alg.jl")
include("src/term.jl")
include("src/term128.jl")
include("src/sorted_linked_list.jl")
include("src/polynomial_dense.jl")
include("src/polynomial_sparse.jl")
include("src/polynomial_sparse_128.jl")
include("src/polynomial_mod_p.jl")
include("src/polynomial_mod_p_128.jl")
include("src/CRT/CRT_functionality.jl")
include("src/CRT/CRT_factor.jl")
    include("src/basic_polynomial_operations/polynomial_addition.jl")
    include("src/basic_polynomial_operations/polynomial_multiplication.jl")
    include("src/basic_polynomial_operations/polynomial_division.jl")
    include("src/basic_polynomial_operations/polynomial_gcd.jl")
include("src/polynomial_factorization/factor.jl")

nothing