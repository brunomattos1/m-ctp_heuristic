include("data.jl")

using ArgParse
import Unicode
using Hygese
using CPUTime
using Random
function parse_commandline(args_array::Array{String,1}, appfolder::String)
    s = ArgParseSettings(usage = "##### VRPSolver #####\n\n"*
        "  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help = false)
    @add_arg_table s begin
       "instance"
          help = "Instance file path"
          default = "data/CAB25.txt"
       "--iter","-i"
          help = "number of iterations"
          arg_type = Int64
          default = 100000000
        "--tempo","-t"
        help = "time limit in seconds"
        arg_type = Int64
        default = 180
        "--model","-m"
        help = "problem model"
        arg_type = Int64
        default = 1
        "--duration","-q"
        help = "problem model"
        arg_type = Int64
        default = 250

    end
    return parse_args(args_array, s)
end

function heuristic(data, iter::Int64, tempo_lim::Int64, q::Int64)
    C, O, M = data.C, data.O, data.M
    bs_O = BitSet(O)
    p_ = data.p
    tempo = 0
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
        tempo += CPUtoq()
        CPUtic()
        return tempo
    end

    function monta_sol(S::Vector{Int})
        # @show S
        cob_S = BitSet()
        for i = 1:length(S)
            cob_S = cob_S ∪ get_covered_by_optional(S[i])
        end
        cont = true
        aux = Set()
        for o in O
            if !isempty(setdiff(get_covered_by_optional(o), cob_S))
                push!(aux, o)
            end
        end
        while cont
            j = rand(aux) # sortear quem cobre, fora do S
            #@show j
            c_j = get_covered_by_optional(j)

            cob_S_j = cob_S ∪ c_j
            if length(cob_S_j) > length(cob_S)
                push!(S, j)
                delete!(aux, j)
                for o in aux
                    if isempty(setdiff(get_covered_by_optional(o), cob_S))
                        delete!(aux, o)
                    end
                end
                cob_S = cob_S_j
            end
            #@show length(cob_S)
            if length(cob_S) >= length(C)
                cont = false
                
                cover_aux = fill!(cover_aux, 0)
                for i in S
                    for c in get_covered_by_optional(i)
                        cover_aux[c] += 1
                    end
                end

                i = 1
                k = length(S)
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
        # @show S
        return S
    end
    
    function perturba(S::Vector{Int}, k::Int64, t::Float64, nm::Int, custoS::Float64, tempo_lim::Int)
        n_m = 0
        while n_m < nm
            CPUtic()
            dist_matS = zeros(length(S))
            j = rand(1:length(S))
            for i = 1:length(S)
                e = (S[i], S[j])
                dist_matS[i] = c(data, e)
            end
            
            kpertos = sortperm(dist_matS, alg=MergeSort)[1:k]
    
            aux = setdiff(S, S[kpertos])
            shuffle!(aux)
            #@show aux
            S′ = monta_sol(aux)
            
            custoS′ = call_cvrp(S′, t)
            
            if custoS′ < custoS - 0.001
                custoS = custoS′
                S = S′
                @show custoS, CPUtok(), n_m, length(S)
                n_m = 0
                
            else
                n_m += 1
            end
            if CPUtok()>tempo_lim
                break
            end
            
        end
        return S, custoS
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
       # dur_limit = 2*maxL + q
        ap = AlgorithmParameters(timeLimit = t, nbIter = 1, mu = 3, lambda = 10) # seconds
        result = solve_cvrp(dist_mat, demand, cap, ap, verbose = false)
        
        return result.cost
    end
    
    
    opt = 100000000.0
    i_ = 0
    CPUtic()
    while i_<iter && CPUtok()<tempo_lim
        
        S = monta_sol(Int[])
        custoS = call_cvrp(S, 1.0)
        

        @show i_, custoS, CPUtok()
        if custoS < opt - 0.001
            opt = custoS
        end
        
        if CPUtok()<tempo_lim                  #
            S′, custoS′ = perturba(S, 3, 1.0, 200, custoS, tempo_lim)
            if custoS′ < opt - 0.001
                opt = custoS′
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
    "$appfolder/540-1/X11119A"

    appfolder = dirname(@__FILE__)
    @show appfolder
    app = parse_commandline(ARGS, appfolder)
    instance_name = split(basename(app["instance"]), ".")[1]
    data = read_data(app)
    opt, tempo, i_ = heuristic(data, app["iter"], app["tempo"], app["duration"])
    @show instance_name, i_, opt, tempo
end
main()
