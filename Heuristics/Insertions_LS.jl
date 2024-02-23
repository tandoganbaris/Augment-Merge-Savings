function compute_profit_relocate(cur_sol::CVRPSolution, tournr1::Int, tournr2::Int, pos1::Int, pos2::Int, instance::CVRPInstance)
    # relocating one node from tour 1 into tour2 
    #current state: (i-1)=>i=>(i+1)      (j-1)=>j=>(j+1)
    #wished state: (i-1)=>(i+1)         (j-1)=>j=>i=>(j+1)
    #THIS FUNCTION: checks profit of this move
	

    # Take the vectors of nodes from tours for calculating. (depot excluded)
    cur_tour1 = cur_sol.tour_plan[tournr1].nodes
    cur_tour2 = cur_sol.tour_plan[tournr2].nodes

    return instance.distancematrix(cur_tour1[pos1 - 1].id, cur_tour1[pos1].id) 
    + instance.distancematrix(cur_tour1[pos1].id, cur_tour1[pos1 + 1].id) 
    + instance.distancematrix(cur_tour2[pos2].id, cur_tour2[pos2 + 1].id)
    - instance.distancematrix(cur_tour2[pos2].id, cur_tour1[pos1].id) 
    - instance.distancematrix(cur_tour1[pos1].id, cur_tour2[pos2 + 1].id) 
    - instance.distancematrix(cur_tour1[pos1 - 1].id, cur_tour1[pos1 + 1].id)
end



function check_feasibility_relocate(cur_sol::CVRPSolution,tournr1::Int, tournr2::Int, pos1::Int, pos2::Int, instance::CVRPInstance)::Bool
    # relocating one node from tour 1 into tour2 
    #current state: (i-1)=>i=>(i+1)      (j-1)=>j=>(j+1)
    #wished state: (i-1)=>(i+1)         (j-1)=>j=>i=>(j+1)
    #THIS FUNCTION: checks demand feasibility of this move
    cur_tour1 = cur_sol.tour_plan[tournr1].nodes

    #find the corresponding node from tour1 
    nodetoinsert = findfirst(n -> n.id == cur_tour1[pos1].id, instance.nodes)

    # Get demand of tour tournr2
    demand = cur_sol.tour_plan[tournr2].demand + nodetoinsert.demand

    # Check if demand exceeds capacity
    if demand > cur_sol.tour_plan[tournr2].capacity
        return false
    else
        return true
    end
end
function PerformRelocate(cur_sol::CVRPSolution, besttour1::Int, besttour2::Int, bestpos1::Int, bestpos2::Int,instance::CVRPInstance)
    # Change the tours by relocating the vertex at bestpos1 in besttour1
    # after the vertex at bestpos2 in besttour2
    # The order of inserting and erasing if besttour1=besttour2 is important
    
    cur_tour1 = cur_sol.tour_plan[tournr1].nodes
    nodeToRelocate = deepcopy(cur_tour1[bestpos1])
    
    if besttour1 === besttour2 && bestpos1 < bestpos2
        # Warning: you need to do the action later 
        # in the tour first otherwise positions are wrong
        insert!(cur_sol.tour_plan[tournr2].nodes, bestpos2+1, nodeToRelocate)
        filter!(node -> node.id != nodeToRelocate.id,  cur_sol.tour_plan[tournr1].nodes)

    else
        filter!(node -> node.id != nodeToRelocate.id,  cur_sol.tour_plan[tournr1].nodes)
        insert!(cur_sol.tour_plan[tournr2].nodes, bestpos2+1, nodeToRelocate)
    end

    recomputeCostAndDemand!(cur_sol, instance)

    return cur_sol.feasibility
end
function ConstructSolutionByInsertion(best_sol::CVRPSolution, instance::CVRPInstance, mu::Float64, lambda::Float64)
    # Implement an insertion heuristic,
    # the first customer in each tour is the one with the largest demand
    
    cur_tour = createVRPTour(instance.capacity,0,0.0, Vector{Node}())

    cust_inserted = falses(instance.numNodes)
    cust_inserted[instance.destDepot] = true
    cust_inserted[instance.origDepot] = true

    solution = Vector{VRPTour}()

    while any(x -> !x, cust_inserted)
        if length(cur_tour.nodes) == 0
            largest_demand = 0
            best_cust = -1
            for j in 1:instance.numNodes
                if cust_inserted[j] 
                    continue
                end
                if instance.nodes[j].demand > largest_demand #find the highest demand node
                    largest_demand = instance.nodes[j].demand
                    best_cust = j
                end
            end
            insert!(cur_tour.nodes, 1, instance.nodes[best_cust])
            cust_inserted[best_cust] = true
            cur_tour.demand += instance.nodes[best_cust].demand
        else
            best_cust = -1
            best_cost = Inf
            best_pos = -1
            for j in 1:instance.numNodes
                if cust_inserted[j] 
                    continue
                end
                # Check feasibility regarding capacity
                if cur_tour.demand + instance.nodes[j].demand <= instance.capacity && (instance.nodes[j].id != instance.destDepot) 
                    #check cost of inserting the node j to before the node i in the current tour
                    for i in eachindex(cur_tour.nodes[1:end]) #changed from 2:end
                        index_before = (i-1)==0 ? instance.origDepot : i-1 # if before is 0 that means its from depot to i 
                        idbefore = cur_tour.nodes[index_before].id
                        cur_id = cur_tour.nodes[i].id
                        ins_cost = (-instance.distancematrix[idbefore, cur_id] * mu) + # mu * (i-1)=> i
                                   instance.distancematrix[idbefore, j] + instance.distancematrix[j, cur_id] - # (i-1)=>j=>i
                                   (lambda * (instance.distancematrix[instance.origDepot, j] + # lambda* depot=>j=>depot
                                             instance.distancematrix[j, instance.destDepot]))
                        if ins_cost < best_cost
                            best_cost = ins_cost
                            best_cust = j
                            best_pos = i
                        end
                    end
                end
            end
            if best_cust > -1
                insert!(cur_tour.nodes, best_pos, instance.nodes[best_cust])
                cust_inserted[best_cust] = true
                cur_tour.demand += instance.nodes[best_cust].demand
            else
                recomputeCostAndDemand!(cur_tour, instance)
                push!(solution, deepcopy(cur_tour))
                cur_tour = createVRPTour(instance.capacity,0,0.0, Vector{Node}())
            end
        end
    end
    recomputeCostAndDemand!(cur_tour, instance)
    push!(solution, deepcopy(cur_tour))
    
    cur_sol = CVRPSolution(instance, solution)

    # If it is the best solution found so far, update best_sol
    if cur_sol.distance < best_sol.distance || length(best_sol.tour_plan) == 0
        best_sol = deepcopy(cur_sol)
    end

    return best_sol
end
function LocalSearchWithRelocate(cur_sol::CVRPSolution, best_sol::CVRPSolution, instance::CVRPInstance)
    counter = 0
    Bestprofit = 0
    bestpos1 = -1
    bestpos2 = -1
    besttour1 = -1
    besttour2 = -1

    while Bestprofit > 0
        # Reset values before searching for the next swap
        Bestprofit = 0
        bestpos1 = -1
        bestpos2 = -1
        besttour1 = -1
        besttour2 = -1

        for tournr1 in 1:length(cur_sol.tour_plan)
            cur_tour1 = cur_sol.tour_plan[tournr1].nodes
            for pos1 in 2:length(cur_tour1) - 1
                for tournr2 in 1:length(cur_sol.tour_plan)
                    cur_tour2 = cur_sol.tour_plan[tournr2].nodes
                    for pos2 in 1:length(cur_tour2) - 1
                        if tournr1 == tournr2 && (pos2 == pos1 || pos2 == pos1 - 1)
                            continue
                        end
                        
                        profit = compute_profit_relocate(cur_sol, tournr1, tournr2, pos1, pos2, instance)
                        if profit > Bestprofit
                            isFeasible = true
                            if tournr1 != tournr2
                                isFeasible =  check_feasibility_relocate(cur_sol, tournr1, tournr2, pos1, pos2, instance)
                            end
                            if isFeasible
                                Bestprofit = profit
                                bestpos1 = pos1
                                bestpos2 = pos2
                                besttour1 = tournr1
                                besttour2 = tournr2
                            end
                        end
                    end
                end
            end
        end

        # Do best move
        if Bestprofit > 0
            counter += 1
            println("Relocate customer ", cur_sol.tour_plan[besttour1].nodes[bestpos1].id, 
                    " after customer ", cur_sol.tour_plan[besttour2].nodes[bestpos2].id)
            PerformRelocate(cur_sol, besttour1, besttour2, bestpos1, bestpos2, instance)

            if !cur_sol.feasibility
                println("Error: tour is infeasible after relocate")
            end
        end
    end

    # Set BestSolution if this one is better
    if cur_sol.distance < best_sol.distance
        best_sol = deepcopy(cur_sol)
    end

    return counter
end

function ComputeProfitSwap(cur_sol::CVRPSolution, tournr1::Int, tournr2::Int, pos1::Int, pos2::Int, instance::CVRPInstance)
    cur_tour1 = cur_sol.tour_plan[tournr1].nodes
    cur_tour2 = cur_sol.tour_plan[tournr2].nodes

    return instance.distancematrix[cur_tour1[pos1 - 1].id, cur_tour1[pos1].id] +
           instance.distancematrix[cur_tour1[pos1].id, cur_tour1[pos1 + 1].id] +
           instance.distancematrix[cur_tour2[pos2 - 1].id, cur_tour2[pos2].id] +
           instance.distancematrix[cur_tour2[pos2].id, cur_tour2[pos2 + 1].id] -
           instance.distancematrix[cur_tour1[pos1 - 1].id, cur_tour2[pos2].id] -
           instance.distancematrix[cur_tour2[pos2].id, cur_tour1[pos1 + 1].id] -
           instance.distancematrix[cur_tour2[pos2 - 1].id, cur_tour1[pos1].id] -
           instance.distancematrix[cur_tour1[pos1].id, cur_tour2[pos2 + 1].id]
end

function CheckFeasibilitySwap(cur_sol::CVRPSolution, routenr1::Int, routenr2::Int, pos1::Int, pos2::Int,instance::CVRPInstance)
    cur_tour1 = cur_sol.tour_plan[routenr1].nodes
    cur_tour2 = cur_sol.tour_plan[routenr2].nodes
    node1demand = cur_tour1[pos1].demand
    node2demand = cur_tour2[pos2].demand

    if routenr1 != routenr2
        demand1 = cur_sol.tour_plan[routenr1].demand - (node1demand) + (node2demand)
        demand2 = cur_sol.tour_plan[routenr2].demand + (node1demand) - (node2demand)

        if demand1 > instance.capacity || demand2 > instance.capacity
            return false
        end
    end

    return true
end
function PerformSwap(cur_sol::CVRPSolution, besttour1::Int, besttour2::Int, bestpos1::Int, bestpos2::Int,instance::CVRPInstance)
    # Change the tours by relocating the vertex at bestpos1 in besttour1
    # after the vertex at bestpos2 in besttour2
    # The order of inserting and erasing if besttour1=besttour2 is important
    
    #for debugging 
    #tour1=cur_sol.tour_plan[besttour1]
    #tour2=cur_sol.tour_plan[besttour2]

    nodeToSwap1 = deepcopy(cur_sol.tour_plan[besttour1].nodes[bestpos1])
    nodeToSwap2 = deepcopy(cur_sol.tour_plan[besttour2].nodes[bestpos2])

    insert!(cur_sol.tour_plan[besttour1].nodes, bestpos1+1, nodeToSwap2)
    deleteat!(cur_sol.tour_plan[besttour1].nodes, bestpos1)
    insert!(cur_sol.tour_plan[besttour2].nodes, bestpos2+1, nodeToSwap1)
    deleteat!(cur_sol.tour_plan[besttour2].nodes, bestpos2)


    recomputeCostAndDemand!(cur_sol, instance)

    return cur_sol.feasible
end

function LocalSearchWithSwap(cur_sol::CVRPSolution, best_sol::CVRPSolution, instance::CVRPInstance)
    if !cur_sol.feasible
        println("Current Solution is infeasible!")
    end
    
    
    
    counter = 0
    Bestprofit = 1
    bestpos1 = -1
    bestpos2 = -1
    besttour1 = -1
    besttour2 = -1

    while Bestprofit > 0
        # Reset values before searching for the next swap
        Bestprofit = 0
        bestpos1 = -1
        bestpos2 = -1
        besttour1 = -1
        besttour2 = -1

        for tournr1 in 1:length(cur_sol.tour_plan)
            cur_tour1 = cur_sol.tour_plan[tournr1].nodes
            for pos1 in 2:length(cur_tour1) - 1
                if pos1 + 2 < length(cur_tour1) #we check swapping position of node inside a tour (?)
                    profit = instance.distancematrix[cur_tour1[pos1 - 1].id, cur_tour1[pos1].id]  +
                    instance.distancematrix[cur_tour1[pos1].id, cur_tour1[pos1 + 1].id]  +
                    instance.distancematrix[cur_tour1[pos1+1].id, cur_tour1[pos1 + 2].id] -
                             instance.distancematrix[cur_tour1[pos1 - 1].id, cur_tour1[pos1 + 1].id] -
                             instance.distancematrix[cur_tour1[pos1 + 1].id, cur_tour1[pos1].id] -
                             instance.distancematrix[cur_tour1[pos1].id, cur_tour1[pos1 + 2].id]

                    if profit > Bestprofit
                        Bestprofit = profit
                        bestpos1 = pos1
                        bestpos2 = pos1 + 1
                        besttour1 = tournr1
                        besttour2 = tournr1
                    elseif profit<-100000
                        println("something is wrong in swap")
                    end
                end
                startpos = pos1 + 2
                for tournr2 in tournr1:(length(cur_sol.tour_plan))
                    cur_tour2 = cur_sol.tour_plan[tournr2].nodes
                    for pos2 in startpos:(length(cur_tour2) - 1)
                        profit = ComputeProfitSwap(cur_sol, tournr1, tournr2, pos1, pos2, instance)
                        if profit > Bestprofit
                            isFeasible = CheckFeasibilitySwap(cur_sol, tournr1, tournr2, pos1, pos2,instance)
                            if isFeasible
                                Bestprofit = profit
                                bestpos1 = pos1
                                bestpos2 = pos2
                                besttour1 = tournr1
                                besttour2 = tournr2
                            elseif profit<-100000
                                println("something is wrong in swap")
                            end
                        end
                    end
                    startpos = 2
                end
            end
        end

        if Bestprofit > 0
            counter += 1
            println("Swap customer $(cur_sol.tour_plan[besttour1].nodes[bestpos1].id) and customer $(cur_sol.tour_plan[besttour2].nodes[bestpos2].id)")

            PerformSwap(cur_sol, besttour1, besttour2, bestpos1, bestpos2, instance)

            if !cur_sol.feasible
                println("Error: tour is infeasible after swap")
            end
        else
            break
        end
    end

    if cur_sol.distance < best_sol.distance
        best_sol = deepcopy(cur_sol)
    end

    return counter , best_sol
end

function simpleLS(instance::CVRPInstance)
    # Initialize counter and random values for mu and lambda
    counter = -1
    lambda = rand(0.0:0.01:2.0)  # Random value between 0.0 and 2.0, adjust step as needed
    mu = rand(0.0:0.01:2.0)      # Same here

    # Start timing
    start_time = time()
    cur_sol = CVRPSolution(instance, Vector{VRPTour}())
    # Insertion Heuristic
    cur_sol= ConstructSolutionByInsertion(cur_sol,instance,mu, lambda)
    best_sol = deepcopy(cur_sol)
    println("Starting solution:")
    screenOutput(best_sol)  
    println("_________________________________________________________________________\n")

    
    # Local Search with Swap
    counter, best_sol = LocalSearchWithSwap(cur_sol, best_sol,instance)
    println("$counter Improving moves were performed with Swap. New solution:")
    screenOutput(best_sol)
    println("_________________________________________________________________________\n")
    #=
    # Local Search with OrOpt
    counter = localSearchWithOrOpt(3)  # Assuming 3 is the parameter for the Or-opt function
    println("$counter Improving moves were performed with Or-opt(k=3). New solution:")
    screenOutput(getBestSolution())
    println("_________________________________________________________________________\n")
    =#
    # Stop timing
    end_time = time()
    elapsed_time = end_time - start_time
    println("\nOptimization took $elapsed_time seconds")
end
