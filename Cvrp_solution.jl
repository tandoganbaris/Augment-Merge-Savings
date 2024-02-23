mutable struct CVRPSolution
    tour_plan::Vector{VRPTour}
    distance::Float64
    feasible::Bool
    instance::CVRPInstance
end

function CVRPSolution(instance::CVRPInstance, tour_plan::Vector{VRPTour})
    solution = CVRPSolution(tour_plan, 0, true, instance)
    recomputeCostAndDemand!(solution, instance)
    return solution
end

function recomputeCostAndDemand!(solution::CVRPSolution, instance)
    if(length(solution.tour_plan)>0)

        for tour in solution.tour_plan
            recomputeCostAndDemand!(tour, instance)
        end
        solution.distance = sum(tour.distance for tour in solution.tour_plan)

        solution.feasible = all(tour.feasible for tour in solution.tour_plan)&& checkAllCustomersVisited(solution, instance)
    end
end
function checkAllCustomersVisited(solution::CVRPSolution, instance)::Bool
    all_nodes_with_demand = [node.id for node in instance.nodes if node.demand > 0]
    seen_ids = Vector{Int}()

    for node_id in all_nodes_with_demand
        found = false
        for tour in solution.tour_plan
            if any(node -> node.id == node_id, tour.nodes)
                if node_id in seen_ids
                    #println("Node $node_id appears more than once in the tours.")
                    return false
                end
                found = true
                push!(seen_ids, node_id)
                break
            end
        end
        if !found
            #println("Node $node_id with positive demand is not included in any tour.")
            return false
        end
    end
    #println("All nodes with positive demand are included in the tours and appear only once.")
    return true


end

function screenOutput(solution::CVRPSolution)
    println(solution.feasible ? "The tour plan is feasible " : "The tour plan is NOT feasible ", 
    "with overall cost: ", solution.distance)
    for tour in solution.tour_plan
        show(tour)
    end
  
end
#=
function Base.show(io::IO, solution::CVRPSolution)
    println(solution.feasible ? "The tour plan is feasible " : "The tour plan is NOT feasible ", 
            "with overall cost: ", solution.distance)
    for tour in solution.tour_plan
        show(tour)
    end
end
=#

