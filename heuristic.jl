include("data.jl") # informations about the instance

using ArgParse
import Unicode
using Hygese
using CPUTime
using Random
function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage = "##### m-CTP #####\n\n"*
        "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help = false)
    @add_arg_table s begin
       "instance"
          help = "Instance file path"
       "--iter","-i"
          help = "number of iterations"
          arg_type = Int64
          default = 100000000
        "--time","-t"
        help = "time limit in seconds"
        arg_type = Int64
        default = 180
        "--duration","-q"
        help = "problem model"
        arg_type = Int64
        default = 250

    end
    return parse_args(args_array, s)
end

function heuristic(data, iter::Int64, time_lim::Int64, q::Int64)
    C, O, M = data.C, data.O, data.M
    bs_O = BitSet(O)
    p_ = data.p
    time = 0
    all_covered = [BitSet(cover_by_optional(data, o)) for o in O]
    get_covered_by_optional(o) = all_covered[o - length(M)]
    all_covered2 = [cover_by_optional(data, o) for o in O]
    get_covered_by_optional2(o) = all_covered2[o - length(M)]
    cover_aux = zeros(length(C) + length(O) + length(M))

    maxL = 0.0
    for i in vcat(O,M)
        e = (0,i)
        if c(data,e) > maxL
            maxL = c(data,e)
        end
    end
    function CPUtok()
        time += CPUtoq()
        CPUtic()
        return time
    end

    function build_sol(S::Vector{Int})
        cov_S = BitSet()
        for i = 1:length(S)
            cov_S = cov_S ∪ get_covered_by_optional(S[i])
        end
        cont = true
        aux = Set()
        for o in O
            if !isempty(setdiff(get_covered_by_optional(o), cob_S))
                push!(aux, o)
            end
        end
        while cont
            j = rand(aux) # select a random facility to enter S
            c_j = get_covered_by_optional(j)

            cov_S_j = cov_S ∪ c_j
            if length(cov_S_j) > length(cov_S)
                push!(S, j)
                delete!(aux, j)
                for o in aux
                    if isempty(setdiff(get_covered_by_optional(o), cov_S))
                        delete!(aux, o)
                    end
                end
                cov_S = cov_S_j
            end
            if length(cov_S) >= length(C) # if all customers are covered
                cont = false
                
                cover_aux = fill!(cover_aux, 0)
                for i in S
                    for c in get_covered_by_optional(i)
                        cover_aux[c] += 1
                    end
                end

                i = 1
                k = length(S)
                # removing redundant facilities
                while i < k
                    can_delete = true
                    for c in get_covered_by_optional2(S[i])
                        if cover_aux[c] == 1
                            can_delete = false
                        end
                    end
                    if can_delete
                        for c in get_covered_by_optional2(S[i])
                            cover_aux[c] -= 1
                        end
                        deleteat!(S, i)
                        k -= 1
                    end
                    i += 1
                end                
            end
        end
        return S
    end
    
    function perturbate(S::Vector{Int}, k::Int64, t::Float64, nm::Int, costS::Float64, time_lim::Int)
        n_m = 0
        while n_m < nm
            CPUtic()
            dist_matS = zeros(length(S))
            j = rand(1:length(S))
            for i = 1:length(S)
                e = (S[i], S[j])
                dist_matS[i] = c(data, e)
            end
            
            kclose = sortperm(dist_matS, alg=MergeSort)[1:k]
    
            aux = setdiff(S, S[kclose])
            shuffle!(aux)
            S′ = build_sol(aux)
            
            costS′ = call_cvrp(S′, t)
            
            if costS′ < costS - 0.001
                costS = costS′
                S = S′
                @show costS, CPUtok(), n_m, length(S)
                n_m = 0
                
            else
                n_m += 1
            end
            if CPUtok()>time_lim
                break
            end
            
        end
        return S, costS
    end
    
    
    
    function call_cvrp(S::Vector{Int}, t::Float64)
        V = vcat([0], M, S)
        dist_mat = zeros(length(V), length(V))
    
        for i = 1:length(V)
            for k = 1:length(V)
                e = (V[i], V[k])
                dist_mat[i,k]  =  c(data, e)
            end
        end
    
        demand = Float64[]
        push!(demand, 0.0)
        for i = 2:length(V)
            push!(demand, 1.0)
        end
        cap = p_
        veh = length(V)
        ap = AlgorithmParameters(timeLimit = t, nbIter = 1, mu = 3, lambda = 10) # seconds
        result = solve_cvrp(dist_mat, demand, cap, ap, verbose = false)
        
        return result.cost
    end
    
    
    opt = 100000000.0
    i_ = 0
    CPUtic()
    while i_<iter && CPUtok()<time_lim
        
        S = build_sol(Int[])
        costS = call_cvrp(S, 1.0)
        

        @show i_, costS, CPUtok()
        if costS < opt - 0.001
            opt = costS
        end
        
        if CPUtok()<time_lim                  #
            S′, costS′ = perturbate(S, 3, 1.0, 200, costS, time_lim)
            if costS′ < opt - 0.001
                opt = costS′
                S = S′
            end
        end
        
        i_ += 1
    end
    @show i_, opt, CPUtok()
    return opt, CPUtok(), i_
end


function main()
    appfolder = dirname(@__FILE__)
    app = parse_commandline(ARGS, appfolder)
    instance_name = split(basename(app["instance"]), ".")[1]
    data = read_data(app)
    opt, time, i_ = heuristic(data, app["iter"], app["time"], app["duration"])
    @show instance_name, i_, opt, time
end
main()
