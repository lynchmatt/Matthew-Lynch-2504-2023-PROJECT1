using DataStructures, Random

# This code here is for meant to help with Project 1.
# It is not code that should be directly used in the project.
# Students can rather copy the relevant bits and incorporate it in their own code.

# It is basically a form of "wrapping" of the MutableLinkedList type from DataStructures.jl, which supports staying sorted, iteration, and dictionary access.

#############
# iteration #
#############

## copy from datastructures
# Base.iterate(l::MutableLinkedList) = l.len == 0 ? nothing : (l.node.next.data, l.node.next.next)
# Base.iterate(l::MutableLinkedList, n::ListNode) = n === l.node ? nothing : (n.data, n.next)

# """
# Makes the linked list iterable. 
# """
# function iterate(lst::MutableLinkedList{T}, state=lst.node.next) where T
#     state == lst.node.prev ? nothing : (lst.node.data, lst.node.next)
# end


# """
# Makes the linked list iterable. 
# """
# function iterate(lst::MutableLinkedList{T}, state=lst.node.next) where T <: Term
#     state == lst.node.prev ? nothing : (lst.node.data, lst.node.next)
# end



#########################
# Demonstrate iteration #
#########################


# for d in l
#     print(d, "\t")
# end
# println()

##################################################
# Enforcing sorting and connecting to dictionary #
##################################################

"""
Assumes that the element `V` has an order (i.e. you can use <).

Assumes that the linked list `lst` is already sorted.

Assumes the dictionary `dict` does not have the key `key`.

Inserts the new element `value` to the list in its sorted position and updates the dictionary with `key` to point at the node of the linked list.
"""
function insert_sorted!(    lst::MutableLinkedList{Term}, 
                            dict::Dict{Int, DataStructures.ListNode{Term}},
                            key::Int,
                            value::Term)::Nothing

    #Note that MutableLinkedList is implemented as a doubly pointed linked list 
    #The element lst.node is a root which points at the head and the tail.
    # lst.node.prev is the last element of the list
    # lst.node.next is the first element of the list

    haskey(dict, key) && error("Key is already in dict")
    
    #If list is empty or the value is greater than end of list, push at end
    if isempty(lst) || last(lst) <= value
        push!(lst, value)
        dict[key] = lst.node.prev #point to last since value just added to last
        return nothing
    end
    
    #if here then lst is not empty
    current_node = lst.node.next #point at first node
    i = 1
    #iterate to find 
    while current_node.data <= value
        if current_node == lst.node.prev #if got to last node
            push!(lst, value) #just push at end
            dict[key] = lst.node.prev #point to last
            return nothing
        end
        current_node = current_node.next #move to next node
    end

    #if here then current_node points at right place
    new_node = DataStructures.ListNode{Term}(value) #create a new node

    #tie new_node between current_node.prev and current_node
    new_node.prev = current_node.prev 
    new_node.next = current_node

    #tie prev to new_node
    new_node.prev.next = new_node

    #tie next to new_node
    new_node.next.prev = new_node

    lst.len += 1
    dict[key] = new_node
    return nothing
 end


"""
Assumes the dictionary `dict` has the key `key`. and it is pointing at the linked list.
"""
function delete_element!(   lst::MutableLinkedList{Term}, 
                            dict::Dict{Int, DataStructures.ListNode{Term}},
                            key::Int)::Nothing

    haskey(dict, key) || error("Key is not in dict")

    node = dict[key]
    delete!(dict, key)

    node.prev.next = node.next
    node.next.prev = node.prev
    lst.len -= 1

    return nothing
end

"""
Returns the value associated with the key `key` or `nothing` if not in list.
"""
get_element(    lst::MutableLinkedList{Term}, 
                dict::Dict{Int, DataStructures.ListNode{Term}},
                key::Int) = haskey(dict, key) ? dict[key].data : nothing

"""
Returns the element following the element associated with the key `key` or `nothing` if there is no next element.
Throws an error if `key` is not available.
"""
function get_next_element(  lst::MutableLinkedList{Term}, 
                            dict::Dict{Int, DataStructures.ListNode{Term}},
                            key::Int) 
    
    haskey(dict, key) || error("Key is not in dict")
    node = dict[key]
    node.next == lst.node && return nothing
    return node.next.data
end

####################################################################
# Demonstrating enforcing sorting and connecting to dictionary     #
# In this example the values are strings and the keys are integers #
####################################################################

# l = MutableLinkedList{String}()
# d = Dict{Int, DataStructures.ListNode{String}}()
# Random.seed!(0)
# rand_string() = String(rand('a':'z', 5))

# for i in 10:10:100
#     value = rand_string()
#     @show (i, value)
#     insert_sorted!(l, d, i, value)
# end

# @show l
# @show d[50]
# @show d[90]
# @show get_element(l, d, 50)
# delete_element!(l, d, 50)
# @show l
# @show length(l)
# @show get_element(l, d, 50)
# @show get_next_element(l, d, 90)
# @show get_element(l, d, 40)
# @show get_next_element(l, d, 40)
# nothing;