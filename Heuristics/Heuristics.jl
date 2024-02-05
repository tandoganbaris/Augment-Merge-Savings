module Heuristics
using DataStructures


#include("../Cvrp_instance.jl")
using .Main.CVRP_instance
#using .Main.Cvrp_instance
#import .Main.Cvrp_instance::CVRPInstance, VRPTour, Node
include("Savings.jl")

export savingsalgo, savingsalgo_output


end