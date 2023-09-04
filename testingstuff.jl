using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly()

show(Term(2,3)) # will display constants correctly when the coefficient is non zero, ie without x^0

y1 = -1x^1 + 2x^3 -3x^0 -4x^2  #will display without +/-1, without powers of zero, and without x^1 when put into polynomial, but not when using show?
y2 = 4x^0 + -3x^2
show(y1)

a = (Term(2,1)) # will only display correctly if i do show(), is this fine?

show(y1)
show(y2)

y3 = -x^3 + 5x + (-5)

show(y3)

# bring up the base.show issue?
# also bring up the 'works' issue?