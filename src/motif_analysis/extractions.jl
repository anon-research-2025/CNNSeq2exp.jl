#=
    Extract the necessary data structures to get the powerful coalitions (motifs)
=#

const B_flatten_t =  @NamedTuple{seq::Int, fil::Int, pos::Int, contribution::float_type, code_mag::float_type}
const nz_dict_t = Dict{Int, Vector{CartesianIndex{2}}} # for CMS

function get_flattened_B(B; lvl=0)
    B_flattened = Vector{B_flatten_t}()
    last_sequence_index = 1
    @inbounds for fil in keys(B[lvl])
        for seq in keys(B[lvl][fil])
            seq > last_sequence_index && (last_sequence_index = seq)
            for pos in keys(B[lvl][fil][seq])
                contribution_here = contribution.(B[lvl][fil][seq][pos]) |> sum # could be negative
                code_mag = B[lvl][fil][seq][pos][1].code_mag
                push!(B_flattened, (seq=seq, fil=fil, pos=pos, contribution=contribution_here, code_mag=code_mag))
            end
        end
    end
    return B_flattened, last_sequence_index
end

get_contribution(b::B_flatten_t) = b.contribution

get_contributions(B_flattened::Vector{B_flatten_t}) = get_contribution.(B_flattened)

get_abs_contributions(B_flattened::Vector{B_flatten_t}) = abs.(get_contribution.(B_flattened))

get_threshold_by_quantile(B_flattened::Vector{B_flatten_t}, _quantile_=0.95) = 
    StatsBase.quantile(get_abs_contributions(B_flattened), _quantile_)

function make_nz_dict_from_B(B_flattened::Vector{B_flatten_t})
    nz_dict = nz_dict_t()
    for b in B_flattened
        entry = CartesianIndex{2}(b.pos, b.fil)
        haskey(nz_dict, b.seq) ? push!(nz_dict[b.seq], entry) : (nz_dict[b.seq] = [entry])
    end
    return nz_dict
end

#=
code_mag_thresh: the threshold for the code magnitude; 
    so that the pattern is not too weak; i.e. noise
contribution_quantile_thresh: the quantile threshold for the contribution; 
    so that the pattern is actually influencing the outcomes; i.e. not noise
    this is only for finding the coalitions (configurations)


    high-level summary:
        nz_dict is processed to save time for motif discovery;
            - use code-magnitude threshold to get all the conserved patterns
            - use contribution threshold to get the positive and negative contributing patterns

        df_B_seq is used to calculate the power indices
            no thresholding is done here; all the contributing components are considered

    more detailed summary:
        nz_dict_positive: used to obtain the positive contributing motifs
        nz_dict_negative: used to obtain the negative contributing motifs
        nz_dict_all: used to obtain the motif discovery -- 
            the motif discovery needs to consider all the motifs occurrences 
            that exceeds the code-magnitude threshold
        
        df_B_seq: the DataFrame grouped by sequence index
            used to calculate the power indices
            needs to consider all the players --- so needs to consider 
                all the motifs that are below the code threshold as well.
                i.e. all the non-zero contributing elements needs to be considered here
=#

round2(x)=round(x,digits=2)

function extract_nz_dicts_and_df_B_seq(
        B_flattened, 
        last_sequence_index; 
        code_mag_quantile_thresh=0.9,
        contribution_quantile_thresh=0.9)
    # B_flattened, last_sequence_index = get_flattened_B(B)
    # get the code-thresholded version of B_flattened
    code_mag_thresh = 
        StatsBase.quantile(map(x->x.code_mag, B_flattened), code_mag_quantile_thresh)
    @info "code_mag_thresh ($code_mag_quantile_thresh quantile): $(round2(code_mag_thresh))"
    B_flattened_code_threshed = filter(x->x.code_mag > code_mag_thresh, B_flattened)
    q = get_threshold_by_quantile(B_flattened_code_threshed, contribution_quantile_thresh)
    @info "contribution_thresh ($contribution_quantile_thresh quantile): $(round2(q))"

    # make a new B_flattened by filtering B_flattened using the quantile q
    B_flattened_q_positive = filter(x->x.contribution > q, B_flattened)
    B_flattened_q_negative = filter(x->x.contribution < -q, B_flattened)
    # map to nz_dict for CMS operations
    nz_dict_positive = make_nz_dict_from_B(B_flattened_q_positive) # seq -> (pos, fil)
    nz_dict_negative = make_nz_dict_from_B(B_flattened_q_negative) # seq -> (pos, fil)

    nz_dict_all = make_nz_dict_from_B(B_flattened_code_threshed)
    df_B = DataFrame(B_flattened)
    df_B_seq = groupby(df_B, :seq)

    @info "nz_dict_all mean size $(round2(StatsBase.mean(length.(values(nz_dict_all)))))"
    @info "df_B_seq mean rows: $(round2(StatsBase.mean(nrow(group) for group in df_B_seq)))"
    @assert df_B_seq[end].seq[1] == last_sequence_index "sequence index mismatch"
    return nz_dict_positive, nz_dict_negative, nz_dict_all, df_B, df_B_seq
end

function get_configs(nz_dict, hp;
    config_num_fil=2, 
    max_config_size = 250,
    cms_delta = 0.001,
    cms_epsilon = 0.0001
    )
    configs, min_count = nothing, 1
    while true
        configs = ShaneGPUCountMinSketch.obtain_enriched_configurations(
            nz_dict, 
            config_num_fil,
            hp.pfm_len;
            min_count=min_count,
            delta=cms_delta, epsilon=cms_epsilon
        )
        if length(configs) > max_config_size
            min_count += 1
        else 
            break
        end
    end
    return configs
end
