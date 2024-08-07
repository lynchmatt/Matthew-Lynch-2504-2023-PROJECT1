using DataStructures, Random

#############
# iteration #
#############

"""
Makes the linked list iterable. 
"""
function iterate(lst::MutableLinkedList{T}, state=lst.node.next) where T 
    state == lst.node.prev ? nothing : (lst.node.data, lst.node.next)
end


##################################################
# Enforcing sorting and connecting to dictionary #
##################################################

"""
Assumes that the element `V` has an order (i.e. you can use <).

Assumes that the linked list `lst` is already sorted.

Assumes the dictionary `dict` does not have the key `key`.

Inserts the new element `value` to the list in its sorted position and updates the dictionary with `key` to point at the node of the linked list.
"""
function insert_sorted!(    lst::MutableLinkedList{V}, 
                            dict::Dict{K, DataStructures.ListNode{V}},
                            key::K,
                            value::V)::Nothing where {K,V}


    haskey(dict, key) && error("Key is already in dict")
    
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
        # if current_node == lst.node.prev #if got to last node
        #     push!(lst, value) #just push at end
        #     dict[key] = lst.node.prev #point to last
        #     return nothing
        # end
        current_node = current_node.next #move to next node
    end

    #if here then current_node points at right place
    new_node = DataStructures.ListNode{V}(value) #create a new node

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
function delete_element!(   lst::MutableLinkedList{V}, 
                            dict::Dict{K, DataStructures.ListNode{V}},
                            key::K)::Nothing where {K,V}

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
get_element(    lst::MutableLinkedList{V}, 
                dict::Dict{K, DataStructures.ListNode{V}},
                key::K) where {K,V} = haskey(dict, key) ? dict[key].data : nothing


