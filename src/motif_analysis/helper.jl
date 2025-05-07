# create SPF
# seq => pos => filter (nested Dict)
function create_spf(nz_components::Dict{Int,Vector{CartesianIndex{2}}})
    spf = Dict{Int, Dict{Int,Set{Int}}}()
    for seq in keys(nz_components)
        for c_ind in nz_components[seq]
            pos, fil = c_ind[1], c_ind[2]
            spf[seq] = get!(spf, seq, Dict{Int, Set{Int}}())
            seq_seq_pos = get!(spf[seq], pos, Set{Int}())
            push!(seq_seq_pos, fil)
        end
    end
    return spf
end

# _infer_size_(xplet_configs) = 
#     length(first(xplet_configs)) ÷ 2 + 1

function config2namedtup_fil(config)
    vals = Vector{String}(undef, length(config)÷2 + 1)
    for i = 1:(length(config)÷2 + 1)
        fi = config[2*(i-1)+1]
        vals[i] = "$fi"
    end
    names = ["m$i" for i in eachindex(vals)]
    return create_named_tuple(names, vals)
end

#=
check if a configuration exists in the SPF
    if it does, then it is part of a motif
    otherwise, it is not
seq: sequence
pos: position in the sequence
config: configuration
ds: number of distances in the configuration
=#
function check_existance(spf, seq, pos, config, ds; i = 1, l = 7)
    @inbounds while i ≤ ds
        next_pos = pos + l + config[2*i]
        next_fil = config[2*(i-1)+3]
        if haskey(spf[seq], next_pos) && next_fil ∈ spf[seq][next_pos]
            (i == ds) && return true
            pos = next_pos
            i += 1
        else
            break
        end
    end
    return false
end

obtain_onehot_subarray(onehotarr, start_pos, end_pos, seq) =
    @view onehotarr[1:4, start_pos:end_pos, 1, seq]

function insert_!(x_motifs, x_tuple_here, x_distance_here, onehot_subarray)
    distances = get!(x_motifs, x_tuple_here, Dict{typeof(x_distance_here), Matrix{float_type}}())
    distances[x_distance_here] = 
        get!(distances, x_distance_here, zeros(float_type, size(onehot_subarray))) .+ onehot_subarray
    return nothing
end

function insert_power_index!(power_info, x_tuple_here, x_distance_here, power_index)
    distances = get!(power_info, x_tuple_here, Dict{typeof(x_distance_here), Vector{float_type}}())
    distances[x_distance_here] = get!(distances, x_distance_here, Vector{float_type}())
    push!(distances[x_distance_here], power_index)
end

function push_placements!(vec_here, seq, placements, is_reversed_complemented::Int)
    for r in eachrow(placements)
        # this is because of 0-based indexing in javascript
        push!(vec_here, (seq - 1, r[1] - 1, r[2], is_reversed_complemented))
    end
end

function add_pos_info!(pos_info, x_tuple_here, x_distance_here, seq::Int, placements)
    if !isnothing(pos_info)
        distances = get!(pos_info, x_tuple_here, 
            Dict{typeof(x_distance_here), Dict{Int, Vector{NTuple{4,Int}}}}())
        distances_x_distance_here = 
            get!(distances, x_distance_here, Vector{NTuple{4,Int}}())
        push_placements!(distances_x_distance_here, seq, placements, 0)
        return nothing
    end
end

function filter_low_count_configurations!(
    x_motifs;
    count_threshold = 0,
)
    for k in keys(x_motifs)
        for d in keys(x_motifs[k])
            count_here = sum(@view x_motifs[k][d][:, 1])
            count_here ≤ count_threshold && delete!(x_motifs[k], d)
            isempty(x_motifs[k]) && delete!(x_motifs, k)
        end
    end
    return x_motifs
end

function update_avg_contrib_and_is_pos_contrib!(_m_, i, contribs_i, x_tuple_here, x_distance_here)
    avg_contrib_dict = get!(_m_.avg_contribs[i], x_tuple_here, Dict{typeof(x_distance_here), float_type}())
    is_pos_contrib_dict = get!(_m_.is_pos_contrib[i], x_tuple_here, Dict{typeof(x_distance_here), Bool}())
    avg_contrib_dict[x_distance_here] = StatsBase.mean(contribs_i[x_tuple_here][x_distance_here]) |> round2
    is_pos_contrib_dict[x_distance_here] = avg_contrib_dict[x_distance_here] ≥ 0 ? true : false
end

function delete_k_d!(_m_, i, config, d; pos_contrib_delete = true)    
    delete!(_m_.pfms[i][config], d) # noise
    delete!(_m_.pos[i][config], d) 
    delete!(_m_.avg_contribs[i][config], d)
    if pos_contrib_delete
        delete!(_m_.pos_contribs[i][config], d)
    else
        delete!(_m_.neg_contribs[i][config], d)
    end
    # delete empty configurations
    if isempty(_m_.pfms[i][config])
        delete!(_m_.pfms[i], config)
        delete!(_m_.pos[i], config)
        delete!(_m_.avg_contribs[i], config)
        if pos_contrib_delete
            delete!(_m_.pos_contribs[i], config)
        else
            delete!(_m_.neg_contribs[i], config)
        end
    end
end

#= 
update the motif struct to get the mean contribution and 
    whether the contribution is positive or negative
=#
function update_avg_contrib_and_is_pos_contrib_all!(_m_)
    # for singletons
    for fil in keys(_m_.pos_contribs[1])
        _m_.avg_contribs[1][fil] = StatsBase.mean(_m_.pos_contribs[1][fil]) |> round2
        _m_.is_pos_contrib[1][fil] = _m_.avg_contribs[1][fil] ≥ 0 ? true : false
    end
    for i = 2:_m_.upto
        for config in keys(_m_.pos_contribs[i])
            for d in keys(_m_.pos_contribs[i][config])
                update_avg_contrib_and_is_pos_contrib!(_m_, i, _m_.pos_contribs[i], config, d)
                !_m_.is_pos_contrib[i][config][d] && delete_k_d!(_m_, i, config, d; pos_contrib_delete=true) # noise
            end
        end
        for config in keys(_m_.neg_contribs[i])
            for d in keys(_m_.neg_contribs[i][config])
                update_avg_contrib_and_is_pos_contrib!(_m_, i, _m_.neg_contribs[i], config, d)
                _m_.is_pos_contrib[i][config][d] && delete_k_d!(_m_, i, config, d; pos_contrib_delete=false)
            end
        end
    end
end


function print_mean_contribution_singletons(_m_; operator=StatsBase.mean)
    for k in keys(_m_.pos_contribs[1])
        mean_here = _m_.pos_contribs[1][k] |> operator
        @info "k: $k, median: $(round2(mean_here))"
    end
end

#=
Rid of the (ambiguous) configurations that are both in positive and negative configs
=#
function rid_of_ambiguous_configs!(pos_configs, neg_configs)
    pos_configs, neg_configs = Set(pos_configs), Set(neg_configs)
    ambiguous_configs = intersect(pos_configs, neg_configs)
    for config in ambiguous_configs
        delete!(pos_configs, config)
        delete!(neg_configs, config)
    end
    return pos_configs, neg_configs
end
