# calculate the power index

const fil_pos_tup_t = @NamedTuple{fil::Int, pos::Int}

#=
Get the list of namedtuples (fil, position) so can query 
the dataframe entries of retrival results (df_B_seq)
=#
function get_list_of_fil_pos(start_pos, config, hp)
    num_fils_config = length(config) ÷ 2 + 1
    list_of_fil_pos = Vector{fil_pos_tup_t}(undef, num_fils_config)
    list_of_fil_pos[1] = (fil = config[1], pos = start_pos)
    for i = 2:num_fils_config
        start_pos += hp.pfm_len + config[2*(i-1)]
        list_of_fil_pos[i] = (fil = config[2*(i-1)+1], pos = start_pos)
    end
    return list_of_fil_pos
end

#=
get the contributions of the target motif and the rest 
=#
function get_target_and_nontarget_contribution(
        list_of_fil_pos, df_B_seq, sequence_index; 
        excluded_fils::Union{Vector{Int},Nothing} = nothing
    )
    this_df = df_B_seq[(sequence_index,)]
    fil_pos_indicator = falses(size(this_df, 1))
    for fil_pos in list_of_fil_pos
        fil, pos = fil_pos.fil, fil_pos.pos
        fil_pos_indicator .|= ((this_df.fil .== fil) .& (this_df.pos .== pos))
    end
    
    # take into account the excluded filters
    excluded_fils_indicator = falses(size(this_df, 1))
    if !isnothing(excluded_fils)
        for fil_pos in list_of_fil_pos
            @assert !(fil_pos.fil ∈ excluded_fils) "excluded fil in list_of_fil_pos"
        end
        for fil in excluded_fils
            excluded_fils_indicator .|= (this_df.fil .== fil)
        end
    end

    target_contributions = this_df[fil_pos_indicator, :contribution]
    non_target_contributions = this_df[(.!fil_pos_indicator) .| excluded_fils_indicator, :contribution]
    # non_target_contributions = this_df[(.!fil_pos_indicator), :contribution]
    @assert length(target_contributions) == length(list_of_fil_pos) "length mismatch"
    return target_contributions, non_target_contributions
end

function obtain_contributions(start_pos, config, sequence_index, df_B_seq, hp; 
    excluded_fils::Union{Vector{Int},Nothing} = nothing
)
    list_of_fil_pos = get_list_of_fil_pos(start_pos, config, hp)
    target_contributions, non_target_contributions = 
        get_target_and_nontarget_contribution(
            list_of_fil_pos, df_B_seq, sequence_index; excluded_fils = excluded_fils)
    return target_contributions, non_target_contributions
end

#=
Calculate the approximate marginal contribution
=#
function calculate_approx_power_index(
        m, # model type
        target_contributions, 
        non_target_contributions;
        sample_size=12500)

    indicator_mat = CUDA.rand(float_type, (length(non_target_contributions), sample_size));
    indicator_mat = indicator_mat .> 0.5
    target_sum = target_contributions |> sum

    nt = indicator_mat .* cu(non_target_contributions);
    coalition_wo = vec(output_function.(m.beta[1] .* sum(nt, dims=1)));
    coalition_w = vec(output_function.(m.beta[1] .* (sum(nt, dims=1) .+ target_sum)));
    power_index = sum(coalition_w .- coalition_wo) / sample_size
    return power_index
end