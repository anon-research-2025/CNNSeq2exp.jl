

function get_start_and_end_pos(i, nz_inds, hp)
    # get the position and filter
    start_pos, fi = nz_inds[i][1], nz_inds[i][2]
    end_pos = start_pos + hp.pfm_len - 1
    return start_pos, end_pos, "$fi"
end

#=
c_mat_num_rows: number of rows in the crosslinking matrix
    - determined by how many datasets are there
=#
function obtain_singleton_motifs!(
    nz_dict_all,
    onehotarr,
    hp,
    singleton_motifs,
    df_B_seq,
    m;
    pos_info = nothing,
    power_info = nothing
)
    L = size(onehotarr, 2)
    for sequence_index in keys(nz_dict_all)
        nz_inds = nz_dict_all[sequence_index] # CartesianIndex{2} (pos, fil)
        len_here = length(nz_inds)
        len_here == 0 && continue

        if len_here ≥ 1
            for i = 1:len_here
                start_pos, end_pos, fi_ind = get_start_and_end_pos(i, nz_inds, hp)
                (start_pos < 1 || end_pos > L) && continue
                onehot_subarray = obtain_onehot_subarray(onehotarr, start_pos, end_pos, sequence_index)

                if haskey(singleton_motifs, fi_ind)
                    singleton_motifs[fi_ind] += onehot_subarray
                else
                    singleton_motifs[fi_ind] = onehot_subarray
                end

                if !isnothing(pos_info)
                    tup = (sequence_index - 1, start_pos - 1, end_pos, 0)
                    haskey(pos_info, fi_ind) ? push!(pos_info[fi_ind], tup) : (pos_info[fi_ind] = [tup])
                end

                if !isnothing(power_info)
                    target_contributions, non_target_contributions = 
                        obtain_contributions(start_pos, (parse(Int, fi_ind),), sequence_index, df_B_seq, hp)
                    power_index_here = calculate_approx_power_index(
                        m, target_contributions, non_target_contributions)

                    power_info_fil = get!(power_info, fi_ind, Vector{float_type}())
                    push!(power_info_fil, power_index_here)
                end
            end
        end
    end
    # filter motifs with low counts
    for k in keys(singleton_motifs)
        if sum((@view singleton_motifs[k][:, 1])) < config_count_threshold
            delete!(singleton_motifs, k)
        end
    end
    return singleton_motifs
end