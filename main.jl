
using TimerOutputs
using Statistics

include("Cvrp_instance.jl")
import .CVRP_instance as cv
include("./Heuristics/Heuristics.jl")
using .Heuristics 
using Revise


filepath = raw"C:\codestuff\TilkVRP\CVRP_heu framework\instances\Vrp-Set-X\X\X-n153-k22.vrp"
currentinstance = cv.CVRPInstance(filepath)



##Savings Test Part

#tours = savingsalgo(currentinstance)
#print(savingsalgo_output(tours))
#=
const to = TimerOutput()

outputs = Dict()
for i in 1:5
    if(i<2)
        tours = savingsalgo(currentinstance)
       # verified = verify_augmentmerge(tours, currentinstance) #verify that all demand is covered
       # println("Verified that all edges are covered: ", verified)
    else
       @timeit to "savingsalgo" tours=savingsalgo(currentinstance)
    
    end
    outputs[i] = savingsalgo_output(tours)
    #println( outputs[i] )
end
print(outputs[1])
show(to)
=#
 simpleLS(currentinstance)