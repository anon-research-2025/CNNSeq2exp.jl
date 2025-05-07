
# define the types to contain the motifs
const pair_motifs_t = NamedTuple{(:m1, :m2),NTuple{2, String}}
const triplet_motifs_t = NamedTuple{(:m1, :m2, :m3),NTuple{3, String}}
const quadruplet_motifs_t = NamedTuple{(:m1, :m2, :m3, :m4),NTuple{4, String}}
const quintuplet_motifs_t = NamedTuple{(:m1, :m2, :m3, :m4, :m5),NTuple{5, String}}
const sextuplet_motifs_t = NamedTuple{(:m1, :m2, :m3, :m4, :m5, :m6),NTuple{6, String}}

const triplet_d_t = NamedTuple{(:d12, :d23),NTuple{2,Int}}
const quadruplet_d_t = NamedTuple{(:d12, :d23, :d34),NTuple{3,Int}}
const quintuplet_d_t = NamedTuple{(:d12, :d23, :d34, :d45),NTuple{4,Int}}
const sextuplet_d_t = NamedTuple{(:d12, :d23, :d34, :d45, :d56),NTuple{5,Int}}

const keys_t = [
    String,
    pair_motifs_t,
    triplet_motifs_t,
    quadruplet_motifs_t,
    quintuplet_motifs_t,
    sextuplet_motifs_t,
]

const values_t = [
    Matrix{float_type},
    Dict{Int,Matrix{float_type}},
    Dict{triplet_d_t,Matrix{float_type}},
    Dict{quadruplet_d_t,Matrix{float_type}},
    Dict{quintuplet_d_t,Matrix{float_type}},
    Dict{sextuplet_d_t,Matrix{float_type}},
]

#=
NTuple{4, Int} refers to (seq, pos_start, pos_end, is_reverse_complemented)
 - This is used for the highlighted positions in the sequence
 - is_reverse_complemented: 0 no, 1 yes
=#
const values_pos_t = [
    Vector{NTuple{4, Int}},
    Dict{Int,Vector{NTuple{4, Int}}},
    Dict{triplet_d_t,Vector{NTuple{4, Int}}},
    Dict{quadruplet_d_t,Vector{NTuple{4, Int}}},
    Dict{quintuplet_d_t,Vector{NTuple{4, Int}}},
    Dict{sextuplet_d_t,Vector{NTuple{4, Int}}},
]

const _contrib_t_ = [
    Vector{float_type},
    Dict{Int,Vector{float_type}},
    Dict{triplet_d_t,Vector{float_type}},
    Dict{quadruplet_d_t,Vector{float_type}},
    Dict{quintuplet_d_t,Vector{float_type}},
    Dict{sextuplet_d_t,Vector{float_type}},
]

# to store the average contribution of a configuration/coalition
const _avg_contrib_t_ = [
    float_type,
    Dict{Int,float_type},
    Dict{triplet_d_t, float_type},
    Dict{quadruplet_d_t, float_type},
    Dict{quintuplet_d_t, float_type},
    Dict{sextuplet_d_t, float_type},
]

# to check if a configuration/coalition contributes positively
const _is_pos_contrib_t_ = [
    Bool,
    Dict{Int, Bool},
    Dict{triplet_d_t, Bool},
    Dict{quadruplet_d_t, Bool},
    Dict{quintuplet_d_t, Bool},
    Dict{sextuplet_d_t, Bool},
]

const _contrib_t_w_exclusions_ = [
    @NamedTuple{excluded::Vector{Int}, contrib::Vector{float_type}},
    @NamedTuple{excluded::Vector{Int}, contrib::Dict{Int,Vector{float_type}}},
    @NamedTuple{excluded::Vector{Int}, contrib::Dict{triplet_d_t,Vector{float_type}}},
    @NamedTuple{excluded::Vector{Int}, contrib::Dict{quadruplet_d_t,Vector{float_type}}},
    @NamedTuple{excluded::Vector{Int}, contrib::Dict{quintuplet_d_t,Vector{float_type}}},
    @NamedTuple{excluded::Vector{Int}, contrib::Dict{sextuplet_d_t,Vector{float_type}}}
]

const freq_matrix_t = Union{[Dict{k,v} for (k, v) in zip(keys_t, values_t)]...}
const freq_pos_t = Union{[Dict{k,v} for (k, v) in zip(keys_t, values_pos_t)]...}
const contrib_t = Union{[Dict{k,v} for (k, v) in zip(keys_t, _contrib_t_)]...}
const contrib_t_w_exclusions = Union{[Dict{k,v} for (k, v) in zip(keys_t, _contrib_t_w_exclusions_)]...}
const avg_contrib_t = Union{[Dict{k,v} for (k, v) in zip(keys_t, _avg_contrib_t_)]...}
const is_pos_contrib_t = Union{[Dict{k,v} for (k, v) in zip(keys_t, _is_pos_contrib_t_)]...}
const weights_t = Union{[Dict{k,Vector{Float64}} for k in keys_t]...}
const distances_t = Union{[Dict{k,Matrix{Float64}} for k in keys_t]...}
const pvals_t = Union{[Dict{k,Float64} for k in keys_t]...}
infer_size(xplet_configs) = length(first(xplet_configs)) ÷ 2 + 1

create_named_tuple(names::Vector{String}, vals::Union{Vector{String}, Vector{Int}}) =
    NamedTuple(zip(Symbol.(names), vals))

function config2namedtup_fil(config, fil2ind, hp)
    vals = Vector{String}(undef, length(config)÷2 + 1)
    for i = 1:(length(config)÷2 + 1)
        fi = config[2*(i-1)+1]
        if fi > hp.M # so filters for RNA problems will never have this modification
            vals[i] = "$(fil2ind[fi-hp.M])r"
        else
            vals[i] = "$(fil2ind[fi])"
        end
    end
    names = ["m$i" for i in eachindex(vals)]
    return create_named_tuple(names, vals)
end

parse_fil_inds(tup) = parse.(Int, (values(tup)))
get_namedtup_ints(tup) = parse_fil_inds(tup) |> unique

#=
c: @NamedTuple of filter indices
d: distance tuple
    Turn a nametuple of filter indices and distances into a config
=#
function to_config(c, d)
    fil_inds = parse_fil_inds(c)
    config = Vector{Int}(undef, length(c) + length(d))
    i, j = 1, 1
    d_scalar = length(d)==1  # Check if d is a scalar    
    for k in 1:length(config)
        if k % 2 == 0
            config[k] = d_scalar ? d : d[j]
            j += 1
        else
            config[k] = fil_inds[i]
            i += 1
        end
    end
    return Tuple(config)
end

function config2namedtup_d(config)
    length(config) == 3 && return config[2]
    vals = [Int(config[2*i]) for i = 1:(length(config)÷2)]
    names = ["d$(i)$(i+1)" for i in eachindex(vals)]
    return create_named_tuple(names, vals)
end

