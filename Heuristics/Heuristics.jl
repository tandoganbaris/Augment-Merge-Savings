module Heuristics
using DataStructures

import .Main.CVRP_instance: CVRPInstance, VRPTour, Node, createVRPTour, recomputeCostAndDemand!

include("../Cvrp_solution.jl")
include("Savings.jl")
include("Insertions_LS.jl")

export savingsalgo, savingsalgo_output, ConstructSolutionByInsertion, LocalSearchWithRelocate, simpleLS


end  # End of module Heuristics