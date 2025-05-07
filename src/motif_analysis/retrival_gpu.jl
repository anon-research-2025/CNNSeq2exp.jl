

#=
fill the placement matrices 
    i.e. each filter's start and end position in the sequence
and return the start and end position of the motif instance

Input:
    pos: start position of the motif instance
    config: configuration
    ds: number of distances in the configuration
    expansions: expansions of the filters

placements: matrix of size (x, 2)
    x: number of filters
    column1 => start position
    column2 => end position
    # TODO change l ? 
=#
function fil_placements_and_return_start_end!(
    pos,
    config,
    placements;
    l = 7,
)
    start_pos, end_pos = nothing, nothing
    num_fils_config = length(config) ÷ 2 + 1
    for i = 1:num_fils_config
        start_pos_f = pos
        end_pos_f = pos + l - 1 
        i == 1 && (start_pos = copy(start_pos_f))
        i == num_fils_config && (end_pos = copy(end_pos_f))
        placements[i, 1] = start_pos_f
        placements[i, 2] = end_pos_f
        if i < num_fils_config
            # no inc because config[2*i] is the distance between 
            # the filters (not including expansions)
            pos = end_pos_f + config[2*i] + 1
        end
    end
    return start_pos, end_pos
end

function obtain_x_plet_motifs!(
    spf,
    onehotarr,
    xplet_configs,
    xplet_motifs,
    df_B_seq,
    m, 
    hp;
    pos_info = nothing,
    power_info = nothing,
)
    isempty(xplet_configs) && return
    L = size(onehotarr, 2)

    # placeholder variable: to keep track of the placements: (start, end) of each motif
    num_placement_rows = infer_size(xplet_configs)
    # num_placement_rows: number of filters in the configuration
    # this placement matrix is pre-allocted to avoid memory allocation
    placements = zeros(Int, (num_placement_rows, 2))
    ds = num_placement_rows - 1 # number of distances in the configuration
    for config in xplet_configs
        xplet_name = config2namedtup_fil(config)
        xplet_dist = config2namedtup_d(config)
        for seq in keys(spf)
            for pos in keys(spf[seq])
                for fil in spf[seq][pos]
                    if fil == config[1] # config[1] is the first filter's index in the configuration
                        if check_existance(spf, seq, pos, config, ds; i = 1, l = hp.pfm_len)
                            start_pos, end_pos = fil_placements_and_return_start_end!(
                                pos,
                                config,
                                placements;
                                l = hp.pfm_len,
                            )
                            (start_pos < 1 || end_pos > L) && continue
                            onehot_subarray = obtain_onehot_subarray(
                                onehotarr,
                                start_pos,
                                end_pos,
                                seq,
                            )
                            insert_!(
                                xplet_motifs,
                                xplet_name,
                                xplet_dist,
                                onehot_subarray,
                            )
                            add_pos_info!(
                                pos_info,
                                xplet_name,
                                xplet_dist,
                                seq,
                                placements,
                            )

                            target_contributions, non_target_contributions = 
                                obtain_contributions(start_pos, config, seq, df_B_seq, hp)
                            power_index_here = calculate_approx_power_index(
                                m, target_contributions, non_target_contributions)
                            insert_power_index!(
                                power_info,
                                xplet_name,
                                xplet_dist,
                                power_index_here,
                            )
                        end
                    end
                end
            end
        end
    end
    filter_low_count_configurations!(xplet_motifs)
end


# values( (m1 = "2", m2 = "2"))
function update_contrib_w_exclusion!(
    pos_info,
    power_info,
    i, # i is the number of filters in the motifs
    df_B_seq,
    opposite_contributors,
    m,
    hp;
    )
    for c in keys(pos_info) # iter through filter configs
        dict_here = Dict{eltype(pos_info[c].keys), Vector{float_type}}()
        for d in keys(pos_info[c]) # iter through distance
            config = to_config(c, d)
            power_indices = Vector{float_type}()
            for j = 1:i:length(pos_info[c][d]) # iter through sequences
                sequence_index = pos_info[c][d][j][1]+1 # pos_info has (seq, pos_start, pos_end, is_reverse_complemented)
                start_pos = pos_info[c][d][j][2]+1      # both +1 to adjust back from 0-based indexing

                target_contributions, non_target_contributions = 
                    obtain_contributions(start_pos, config, sequence_index, df_B_seq, hp; 
                        excluded_fils = opposite_contributors[c])
                power_index_here = calculate_approx_power_index(
                    m, target_contributions, non_target_contributions)

                push!(power_indices, power_index_here)      
            end           
            dict_here[d] = power_indices
        end
        @info "c: $c"
        power_info[c] = (excluded=opposite_contributors[c], contrib=dict_here)
    end
end

# insert_power_index!(
#     power_info,
#     c,
#     d,
#     power_index_here,
# )



