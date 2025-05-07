"""
pfms: A vector where...
        {1st, 2nd,...} component :
        key is the {1st, 2nd,...} entry in keys_t
        value is the {1st, 2nd,...} entry in values_t
        stores the pfms for the training set
similarly for pfms_train_shuffled, pfms_test, pfms_test_shuffled
"""
mutable struct motifs
    pfms::Vector{freq_matrix_t}
    pos_contribs::Vector{contrib_t}
    neg_contribs::Vector{contrib_t}
    pos_contribs_w_exclusions::Vector{contrib_t_w_exclusions}
    neg_contribs_w_exclusions::Vector{contrib_t_w_exclusions}
    avg_contribs::Vector{avg_contrib_t}
    is_pos_contrib::Vector{is_pos_contrib_t}
    pos::Vector{freq_pos_t}
    upto::Int
    """
    data_seq_range: 
        - Vector of UnitRange{Int}, showing 
          different datasets' sequence indices
        - e.g. [1:100, 101:200, 201:300, 301:400]  
        - length(data_seq_range) = how many datasets
    upto: 
        - how many n-ary motifs to consider
        - e.g. upto=3 means singleton, pair, triplet motifs
    """
    function motifs(; upto = 3)
        keys2use = keys_t[1:upto]
        values2use = values_t[1:upto]
        contrib2use = _contrib_t_[1:upto]
        contrib_w_exclusions2use = _contrib_t_w_exclusions_[1:upto]
        avg_contribs2use = _avg_contrib_t_[1:upto]
        is_pos_contrib_t2use = _is_pos_contrib_t_[1:upto]
        pos2use = values_pos_t[1:upto]
        pfms = [Dict{k,v}() for (k, v) in zip(keys2use, values2use)]
        pos_contribs = [Dict{k,v}() for (k, v) in zip(keys2use, contrib2use)]
        neg_contribs = [Dict{k,v}() for (k, v) in zip(keys2use, contrib2use)]
        pos_contribs_w_exclusions = [Dict{k,v}() for (k, v) in zip(keys2use, contrib_w_exclusions2use)]
        neg_contribs_w_exclusions = [Dict{k,v}() for (k, v) in zip(keys2use, contrib_w_exclusions2use)]
        avg_contribs = [Dict{k,v}() for (k, v) in zip(keys2use, avg_contribs2use)]
        is_pos_contrib = [Dict{k,v}() for (k, v) in zip(keys2use, is_pos_contrib_t2use)]
        pos =  [Dict{k,v}() for (k, v) in zip(keys2use, pos2use)]
        
        return new(
            pfms,
            pos_contribs,
            neg_contribs,
            pos_contribs_w_exclusions,
            neg_contribs_w_exclusions,
            avg_contribs,
            is_pos_contrib,
            pos,
            upto
        )
    end
end


# m = motifs(;unto=3)