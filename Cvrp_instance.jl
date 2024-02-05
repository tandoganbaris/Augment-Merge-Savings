module CVRP_instance
using LinearAlgebra
import Base

export Node, CVRPinstance, VRPTour 
mutable struct Node
    X::Int
    Y::Int
    demand::Int
    id::Int
    TimeWindow::Tuple{Int, Int}
end

function Node(X::Int, Y::Int, demand::Int, id::Int)::Node #constructor without time window
    time_window = (0,0)
    return Node(X, Y, demand, id, time_window)
end

mutable struct CVRPInstance
    capacity::Int64	
    destDepot::Int64
    numCustomers::Int64
    numNodes::Int64
    origDepot::Int64
    instanceName::String
    nodes::Vector{Node}
    distancematrix::Matrix{Float64}
end

mutable struct VRPTour
    capacity::Int64
    demand::Int64
    distance::Float64
    nodes::Vector{Node} #careful about including or excluding depot
end
function VRPTour(capacity::Int64, demand::Int64, distance::Float64, nodes::Vector{Node})
    return VRPTour(capacity, demand, distance, nodes)
end

function Base.show(io::IO, tour::VRPTour)
    node_str = join([string(node.id) for node in tour.nodes], ", ")
    println(io, "Tour: [Nodes: $node_str, Total Distance: $(tour.distance), Total Demand: $(tour.demand)]")
end

function read_section(file, keyword::String)
    line = readline(file)
    while !startswith(line, keyword) && !eof(file)
        line = readline(file)
    end
    return line
end



function CVRPInstance(filename::String)::CVRPInstance
    capacity = -1
    numNodes = -1
    nodes = Vector{Node}()
    origDepot = -1
    destDepot = -1

    file = open(filename, "r") 
    
    try
        #line = readline(file)          
        while !eof(file)
    
            linewithdata = read_section(file, "DIMENSION :")
            lineinfo =split(linewithdata)
            numNodes = parse(Int, lineinfo[3])

            linewithdata = read_section(file, "CAPACITY :")
            lineinfo =split(linewithdata)
            capacity = parse(Int, lineinfo[3])

            # Read NODE_COORD_SECTION
           
            try 
                read_section(file, "NODE_COORD_SECTION")
                for _ in 1:numNodes
                    lineinfo = split(readline(file))
                    if !isempty(lineinfo)  
                        id = parse(Int, lineinfo[1])
                        xCoord = parse(Int, lineinfo[2])
                        yCoord = parse(Int, lineinfo[3])
                        newnode = Node(xCoord, yCoord, 0, id)  # Assuming 0 demand initially
                        push!(nodes, newnode)
                    end
                end
                read_section(file, "DEMAND_SECTION")
                for _ in 1:numNodes
                    dem = split(readline(file))
                    if !isempty(dem)  # Check if line contains numeric data
                        id = parse(Int, dem[1])
                        demand = parse(Int, dem[2])
                        for node in nodes
                            if node.id == id
                                node.demand = demand
                                break
                            end
                        end
                    end
                end
        
            catch
                #println( "error while reading nodes")
                continue
            end
            read_section(file, "DEPOT_SECTION")
            origDepot = parse(Int, split(readline(file))[1])
            destDepot = origDepot
        end

    catch
        # Handle the case where readline fails
       # println("Failed to read line. Closing file...")
       
    end

    close(file)
    numCustomers = numNodes - 1
    # Compute costs
    distancematrix = zeros(Float64, numNodes, numNodes)
    for i in 1:numNodes
        for j in 1:numNodes         
            
            if i == j
                distancematrix[i, j] = typemax(Int)
            else
                node1= nodes[i]
                node2= nodes[j]
                distancematrix[i, j] = sqrt((node1.X - node2.X)^2 + (node1.Y - node2.Y)^2)
            end
        end
    end

    #println("Instance $filename was constructed successfully.")
    return CVRPInstance(capacity, destDepot, numCustomers, numNodes, origDepot, filename, nodes, distancematrix)
end
end