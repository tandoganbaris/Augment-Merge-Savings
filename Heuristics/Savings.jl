#=
function SavingsAlgo
   tourlist= Initilize all tours 0-i-0
    calculate savingslist (or dict)
   while savingslist is not empty
    for each saving in savingsmap (highest savings first)
      tour1 = searchtour that has node i right after or right before depot(input tourlist and saving)
      tour2 = searchtour that has node j right after or right before depot(input tourlist and saving)
      mergedtour = empty tour
      bool merge(tour1, tour2, out mergedtour) check if can be merged and if so return true , if not false
      delete saving
      if merge false continue
      else if merge true
         delete tour1 and tour2 from tourlist, add mergedtour
         update savings
    
return tourlist
=#

function savingsalgo(instance)
    tours = initializetours(instance)
    savingsdict = calculatesavings(tours, instance)
    while !isempty(savingsdict)
        key, value = popfirst!(savingsdict) #remove best saving to used
        found1, tour1 = findtour(key[1], tours)
        found2, tour2 = findtour(key[2], tours)
        capacityfits =(tour1.demand + tour2.demand <= instance.capacity)
        if((found1 && found2) &&tour1!=tour2 &&capacityfits)
            newtour = mergetours(instance, key, tour1, tour2)
            savingsdict = updatesavings(key, newtour, savingsdict)
            filter!(tour -> !(tour in [tour1, tour2]), tours)
            push!(tours, newtour)
        end
    end
  
    return tours 
end
function updatesavings(key, newtour, savingsdict)
    remove1=false
    remove2=false
    filtered_dict = Dict{Tuple{Int, Int}, Float64}()
    if !(key[1] == newtour.nodes[1].id || key[1] == newtour.nodes[end].id) 
        remove1=true        
    end 
    if  !(key[2] == newtour.nodes[1].id || key[2] == newtour.nodes[end].id)
        remove2=true
    end
    if remove1 && remove2==false
        for (k, v) in savingsdict
            if key[1] ∉ k
                filtered_dict[k] = v
            end
        end
    elseif remove1==false && remove2
        for (k, v) in savingsdict
            if key[2] ∉ k
                filtered_dict[k] = v
            end
        end

    elseif remove1 && remove2
        for (k, v) in savingsdict
            if key[1] ∉ k && key[2] ∉ k
                filtered_dict[k] = v
            end
        end 
    else 
        return savingsdict
    end
    sorted_pairs = sort([(k, v) for (k, v) in filtered_dict], by=x -> (x[2], -sum(x[1])), rev=true) #
    filtered_dict = DataStructures.OrderedDict(sorted_pairs)
    return filtered_dict

end

function savingsalgo_output(tours)::String
    total_distance = sum(tour.distance for tour in tours)
    output = "\nTotal distance : $total_distance\n"
    for tour in tours
        node_str = join([string(node.id) for node in tour.nodes], ", ")
        thistourstring = "Tour: [Nodes: $node_str, Tour Distance: $(tour.distance), Total Demand: $(tour.demand)] \n"
        output = string(output,thistourstring)
    end

    return output

end

function mergetours(instance, key, tour1::VRPTour, tour2::VRPTour)

newtour = createVRPTour(0,0,0.0,Vector{Node}())

if tour1.nodes[end].id == key[1] && tour2.nodes[1].id == key[2] # ...i <-> j....
    nodes_combined = vcat(deepcopy(tour1.nodes), deepcopy(tour2.nodes))
    total_demand = tour1.demand + tour2.demand
    newtour = createVRPTour(instance.capacity,total_demand,0.0,nodes_combined)
    newtour.distance = calculate_tour_distance(instance, newtour)


elseif tour1.nodes[end].id == key[1] && tour2.nodes[end].id == key[2] # ...i <-> ....j
    nodes_combined = vcat(deepcopy(tour1.nodes), reverse(deepcopy(tour2.nodes)))
    total_demand = tour1.demand + tour2.demand
    newtour = createVRPTour(instance.capacity,total_demand,0.0,nodes_combined)
    newtour.distance = calculate_tour_distance(instance, newtour)

elseif tour1.nodes[1].id == key[1] && tour2.nodes[1].id == key[2] # i... <-> j....
    nodes_combined = vcat(reverse(deepcopy(tour1.nodes)), deepcopy(tour2.nodes))
    total_demand = tour1.demand + tour2.demand
    newtour = createVRPTour(instance.capacity,total_demand,0.0,nodes_combined)
    newtour.distance = calculate_tour_distance(instance, newtour)

elseif tour1.nodes[1].id == key[1] && tour2.nodes[end].id == key[2] # i... <-> ....j
    nodes_combined = vcat(deepcopy(tour2.nodes), deepcopy(tour1.nodes))
    total_demand = tour1.demand + tour2.demand
    newtour = createVRPTour(instance.capacity,total_demand,0.0,nodes_combined)
    newtour.distance = calculate_tour_distance(instance, newtour)
end
#if(newtour.distance> tour1.distance + tour2.distance)
#    here = 0
#end

return newtour
end

function findtour(node::Int, tours::Vector{VRPTour})
    emptytour = createVRPTour(0,0,0.0,Vector{Node}())
    for tour in tours
        if tour.nodes[1].id == node || tour.nodes[end].id == node
            return true, tour
        end
    end
    return false, emptytour
end

function calculatesavings(tours::Vector{VRPTour}, instance) 
    savingsdict = Dict{Tuple{Int, Int}, Float64}()
    distances = instance.distancematrix
    num_nodes = instance.numNodes

    for i in 1:num_nodes
        if(i==instance.destDepot)
            continue
        else
            for j in (i+1):num_nodes
                # Calculate savings
                saving = distances[instance.destDepot , i] + distances[instance.destDepot, j] - distances[i, j]

                # Create the key as a tuple of two integers
                key = (i, j)

                # Add the key-value pair to the dictionary
                savingsdict[key] = saving
            end
        end
    end
    sorted_pairs = sort([(k, v) for (k, v) in savingsdict], by=x -> (x[2], -sum(x[1])), rev=true) #
    savingsdict = DataStructures.OrderedDict(sorted_pairs)



    return savingsdict
end

function initializetours(instance)
    tours = Vector{VRPTour}()
    for node in instance.nodes
        if node.demand>0
            
            newtour = createVRPTour(instance.capacity, node.demand,0.0,[node])
            newtour.distance = 2*instance.distancematrix[instance.destDepot, node.id]
            push!(tours, newtour)
        end
    end



    return tours
    
end

function calculate_tour_distance(instance, tour::VRPTour)
    distance_matrix = instance.distancematrix
    total_distance = 0.0
    num_nodes = length(tour.nodes)
    depotid = instance.destDepot
    
    # Iterate over the tour nodes

    for i in 1:num_nodes -1
        node1 = tour.nodes[i]
        node2 = tour.nodes[i+1]
        total_distance += distance_matrix[node1.id, node2.id]
    end

    #add distances to/from depot
    total_distance += distance_matrix[depotid, tour.nodes[1].id]
    total_distance += distance_matrix[tour.nodes[end].id, depotid]
    
    return total_distance
end