
less_than(x, y) = isless(x, y)
greater_than(x, y) = isless(y, x)

#=
keep only the top (keep) number of opposite contributors
=#
function trim_opp_contributors!(opp_contributors; keep=100, comparitor=less_than)
    for k in keys(opp_contributors)
        len_here = length(opp_contributors[k])
        opp_contributors[k] = sort(opp_contributors[k], lt=comparitor)[1:min(len_here, keep)]
    end
end

#=
i: number of filters in the config 
df_B_seq: DataFrame of B_flattened; contains all the contributions
=#
function get_opposite_contributors(_m_, i, df_B_seq; 
    num_opposite_contributors = 3, keep=2500, opposite="negative")

    @assert opposite in ["negative", "positive"] "opposite must be either negative or positive"
    comparitor = opposite == "negative" ? less_than : greater_than;
    opposite_contributors = Dict{eltype(keys(_m_.pos_contribs[i])), Vector{Int}}() 

    # c: config with filters, e.g. i=2 and (m1="2", m2="2)
    for c in keys(_m_.pos_contribs[i])
        # find the sequence indices under the configurations (same filters, but relax distance)
        sequence_indices_here = Set{Int}()
        for d in keys(_m_.pos[i][c])
            for tup in _m_.pos[i][c][d]
                push!(sequence_indices_here, tup[1]+1) # adjust back from 0 based indexing
            end
        end
        sequence_indices_here = sequence_indices_here |> collect

        # gather all the opposite contributors' contributions
        opp_contributors = Dict{Int, Vector{float_type}}()
        for sequence_index in sequence_indices_here
            df_here = df_B_seq[sequence_index]
            for (fil, contribution) in zip(df_here.fil, df_here.contribution)
                # if comparitor(contribution, 0)
                if contribution < 0
                    haskey(opp_contributors, fil) ? push!(opp_contributors[fil], contribution) : (opp_contributors[fil] = [contribution])
                end
            end
        end

        trim_opp_contributors!(opp_contributors; keep=keep)
        # take the top most (num_opposite_contributors)-opposite-contributors
        k = min(num_opposite_contributors, opp_contributors.count)
        selected_opposite_contributors_indices = 
            partialsortperm(mean.(values(opp_contributors)), 1:k; lt=comparitor)

        selected_opposite_contributors = collect(keys(opp_contributors))[selected_opposite_contributors_indices]
        # make sure the selected opposite contributors are not in the positive contributors
        fil_inds = get_namedtup_ints(c)
        filter!(x->!(x in fil_inds), selected_opposite_contributors)
        opposite_contributors[c] = selected_opposite_contributors
    end

    return opposite_contributors
end


# using HypothesisTests


# quantile(opp_contributors[14], 0.75)
# quantile(opp_contributors[6], 0.75)

# ChisqTest(opp_contributors[14])
