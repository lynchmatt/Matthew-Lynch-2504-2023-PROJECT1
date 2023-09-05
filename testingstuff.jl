using Pkg
Pkg.activate(".")

include("poly_factorization_project.jl")

x = x_poly()

show(Term(2,3)) # will display constants correctly when the coefficient is non zero, ie without x^0


y1 = -1x^1 + 2x^3 -3x^0 -4x^2  #will display without +/-1, without powers of zero, and without x^1 when put into polynomial, but not when using show?
y2 = 1x^0 + -3x^2 # need to print 1 as constant, currently won't
show(y1)

show(Term(-1,1))

a = (Term(1,0)) # will only display correctly if i do show(), is this fine?

show(y1)
show(y2)

y3 =  5x -x^3 + (-5)

show(y3)

import Pkg; Pkg.add("Subscripts")
using Subscripts
super("1")

println(÷)

# bring up the base.show issue?
# ask if it's valid to use an additional package to easily convert to superscripts?


# function globaltesting()
#     if a == false || a === nothing
#         println("a is false or has not been set")
#     else
#         println("a is true")
#     end
# end

# globaltesting()

# introduce another local variable to keep == out of if-else? works
# sort out the global part to make sure works in a function


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

alphabetprinter()
