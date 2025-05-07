get_pval(pvals, k, d) = replace(@sprintf("%.1e", pvals[k][d]), "e-0" => "e-", "e+0" => "e+")
get_pval(pvals, k) = replace(@sprintf("%.1e", pvals[k]), "e-0" => "e-", "e+0" => "e+")
get_pval(pval) = replace(@sprintf("%.1e", pval), "e-0" => "e-", "e+0" => "e+")

get_descriptive_str(k, d::Integer) = "pattern $(k.m1) and $(k.m2) are $(d) nucleotides apart"

get_counts(x_motifs, k, d) = sum(@view x_motifs[k][d][:, 1]) |> Int

function get_descriptive_str(keys_, ds::NamedTuple)
    ks = [(keys_[i], keys_[i+1]) for i = 1:length(keys_)-1]
    str = ""
    for ((m1, m2), d) in zip(ks, ds)
        str *= "pattern $(m1) and $(m2) are $(d) nucleotides apart<br>"     
    end
    return str
end

function make_save_paths(; saved_results_dir = saved_results_dir)
    saved_pairs_dir = joinpath(saved_results_dir, x_dir[2])
    saved_triplets_dir = joinpath(saved_results_dir, x_dir[3])
    saved_quadruplets_dir = joinpath(saved_results_dir, x_dir[4])
    mkpath(saved_pairs_dir)
    mkpath(saved_triplets_dir)
    mkpath(saved_quadruplets_dir)
end

function init_json_dict()
    json_x_motifs =
        Dict{String,Dict{String,Union{String,Vector{String},Vector{Vector{String}}}}}()
    return json_x_motifs
end

function init_dict_for_html_render()
    html_dict = Dict(
        tag_div_img_id => String[],
        tag_i => String[],
        tag_img_src => String[],
        tag_img_alt => String[],
        tag_div_text_id => String[],
        tag_p_id1_default => String[],
        tag_p_id2_default => String[],
        tag_p_id3_default => String[],
        tag_p_id4_default => String[],
        tag_p_id5_default => String[],
        tag_p_id6_default => String[],
        tag_div_slide_id => String[],
        tag_max_comb => String[],
    )
    return html_dict
end


round2percentage(x) = round(round(x, digits = 2) * 100, digits = 0)

################################### fill the cluster message ###################################
function get_this_mode_pval_dist_percentage(m, x::Int, k)
    @assert sum(m.pfms_mode_weights[x][k]) ≈ 1 "Weights do not sum to 1"
    this_mode_dist = Int.(m.pfms_mode_distance[x][k])
    this_mode_percentage = Int.(round2percentage.(m.pfms_mode_weights[x][k]))
    # round.(round.(m.pfms_mode_weights[x][k], digits=2) .* 100, digits=0))
    this_mode_pval = m.pfms_mode_pvals[x][k]
    @assert size(this_mode_dist, 1) == length(this_mode_percentage) "Lengths of distances and percentages do not match"

    sorted_indices = sortperm(this_mode_percentage)
    len = length(sorted_indices)
    this_mode_percentage = @view this_mode_percentage[sorted_indices]
    this_mode_dist = @view this_mode_dist[sorted_indices, :]
    this_mode_colors = PlotPWM.arrow_color_palette[1:len] # note that the arrow plot sort weights from low to high
    return this_mode_pval, this_mode_dist, this_mode_percentage, this_mode_colors
end

function get_pval_msg(pval)
    return "P-value: $((SparseCodeMotifs.get_pval(pval)))"
end

function get_cluster_msg(mode_percentage, mode_here, distances_here, this_mode_colors, ind)
    _percentage_ = "Cluster $ind: Occuring approximately $mode_percentage% of the time with"
    mode_len = length(mode_here)
    _in_bt_ = Array{String}(undef, mode_len - 1)
    for i = 1:mode_len-1
        ip1 = i + 1
        _in_bt_[i] = "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; $(distances_here[i])-bp in between pattern $(mode_here[i]) and $(mode_here[ip1])"
    end
    color_here = this_mode_colors[ind]
    return "<p style=\\'color:$color_here\\'>" *
           join([_percentage_, _in_bt_...], "<br>") *
           "</p>"
end

function get_mode_msg(m, x, k)
    this_mode_pval, this_mode_dist, this_mode_percentage, this_mode_colors =
        get_this_mode_pval_dist_percentage(m, x, k)
    pval_msg = get_pval_msg(this_mode_pval)
    descriptions = join(
        [
            get_cluster_msg(p, k, view(this_mode_dist, ind, :), this_mode_colors, ind)
            for (ind, p) in enumerate(this_mode_percentage)
        ],
        "",
    )
    return join([pval_msg, descriptions], "<br><br>")
end
################################################################################################

fill_csv_link(csv_link) =
    "<a href=\"#\" onclick=\"openHighlightPage(\'" * csv_link * "\')\">string highlight</a>"

fill_cluster_link(img_link, text) =
    "<a href=\"#2\" onclick=\"openHtmlWindow(\'" *
    "$img_link\', \'$text" *
    "\')\">clustered view</a>"

fill_concensus_link(text) = 
    "<a href=\"#3\" onclick=\"openHtmlWindowText(\'$text\')\">consensus string</a>"

fill_gap_img_ink(img_link) =
    "<a href=\"#4\" onclick=\"openHtmlWindowImg(\'$img_link\')\">gap histogram</a>"

fill_box_plot_kd_link(img_link, avg_contribution) =
    "<a href=\"#5\" onclick=\"openHtmlWindowImg(\'$img_link\')\"><span class=\"putBar\">contribution</span>: $avg_contribution</a>"
   

function trim_dash(arr::Vector{Char})
    # Find the first and last indices that are not '-'
    start_idx = findfirst(x -> x != _placeholder_char_, arr)
    end_idx = findlast(x -> x != _placeholder_char_, arr)
     # Return the trimmed array
     if start_idx === nothing || end_idx === nothing
        return ['-','-','-']  # Return a three dash array if all elements are '-'
    else
        return arr[start_idx:end_idx]
    end
end

function get_relaxed_consensus_str(pfm; dash = _placeholder_char_, prob_thresh = 0.5, rna = false)
    argmax_inds = reshape(argmax(pfm, dims = 1), (size(pfm, 2),))
    char_array =
        rna ? [_ind2dna_str_rna[i[1]] for i in argmax_inds] :
        [_ind2dna_str_[i[1]] for i in argmax_inds]
    char_array[findall((@view pfm[argmax_inds]) .< prob_thresh)] .= dash
    return join(trim_dash(char_array))
end

function countmat2consensus(
    countmat::AbstractMatrix;
    pseudo_count = float_type(0.1),
    rna = false,
)
    _countmat_ = countmat .+ pseudo_count
    pfm = _countmat_ ./ sum(_countmat_, dims = 1)
    consensus_str = get_relaxed_consensus_str(pfm; rna = rna)
    return consensus_str
end

function update_json_and_dict!(
    m,
    json_x_motifs,
    html_dict;
    t = 1,
    sorted_configs = nothing,
    rna = false,
)
    x_motifs = m.pfms[t]
    x_avg_contribs = m.avg_contribs[t]
    keys2iter = isnothing(sorted_configs) ? keys(x_motifs) : sorted_configs

    for (ind, k) in enumerate(keys2iter)
        mode_ind = "mode_$ind"
        json_x_motifs[mode_ind] = Dict{String,Vector{String}}()
        json_x_motifs[mode_ind][pwms_str] = Vector{String}()
        json_x_motifs[mode_ind][labels_str] = Vector{String}()
        json_x_motifs[mode_ind][texts_str] = Vector{Vector{String}}()

        if t == 1 # singleton motifs
            # @info "k: $k"
            push!(json_x_motifs[mode_ind][pwms_str], get_pwm_singleton_save_path(k))
            push!(json_x_motifs[mode_ind][labels_str], "pattern $(k)")
            counts_here = (@view x_motifs[k][:, 1]) |> sum |> Int
            meme_link = "<a href=\"$(x_dir[1])/$k.meme\">.meme file</a>"
            consensus_str = countmat2consensus(x_motifs[k]; rna = rna)
            csv_str = fill_csv_link(get_csv_singleton_save_path(k))
            push!(
                json_x_motifs[mode_ind][texts_str],
                [
                    fill_box_plot_kd_link(get_box_singleton_save_path(k), x_avg_contribs[k]),
                    "# occurrences: $counts_here",
                    meme_link,
                    csv_str,
                    consensus_str,
                    ""
                ],
            )
        else
            keys_here = collect(keys(x_motifs[k]))
            for d in sort(keys_here)
                push!(json_x_motifs[mode_ind][pwms_str], get_pwm_save_path(k, d; t = t))
                push!(json_x_motifs[mode_ind][labels_str], get_descriptive_str(k, d))
                counts_here = get_counts(x_motifs, k, d)
                meme_link = "<a href=\"$(get_pwm_save_folder(k))/$(get_d_str(d)).meme\">.meme file</a>"
                consensus_str = fill_concensus_link(countmat2consensus(x_motifs[k][d]; rna = rna))
                gapped_histogram = length(keys_here) == 1 ? "--" : fill_gap_img_ink(get_mode_gap_hist_save_path(k))
                csv_str = fill_csv_link(get_csv_save_path(k, d; t = t))
                push!(
                    json_x_motifs[mode_ind][texts_str],
                    [
                        fill_box_plot_kd_link(get_mode_boxplot_save_path(k, d), x_avg_contribs[k][d]),
                        "# occurrences: $counts_here",
                        meme_link,
                        csv_str,
                        consensus_str,
                        gapped_histogram
                    ],
                )
            end
        end
        push!(html_dict[tag_i], "$ind")
        push!(html_dict[tag_div_img_id], "imageContainer$ind")
        push!(html_dict[tag_img_src], json_x_motifs[mode_ind][pwms_str][1])
        push!(html_dict[tag_img_alt], json_x_motifs[mode_ind][labels_str][1])
        push!(html_dict[tag_div_text_id], "textContainer$ind")
        push!(html_dict[tag_p_id1_default], json_x_motifs[mode_ind][texts_str][1][1])
        push!(html_dict[tag_p_id2_default], json_x_motifs[mode_ind][texts_str][1][2])
        push!(html_dict[tag_p_id3_default], json_x_motifs[mode_ind][texts_str][1][3])
        push!(html_dict[tag_p_id4_default], json_x_motifs[mode_ind][texts_str][1][4])
        push!(html_dict[tag_p_id5_default], json_x_motifs[mode_ind][texts_str][1][5])
        push!(html_dict[tag_p_id6_default], json_x_motifs[mode_ind][texts_str][1][6])
        push!(html_dict[tag_div_slide_id], "slideContainer$ind")
        push!(html_dict[tag_max_comb], "$(length(json_x_motifs[mode_ind][pwms_str])-1)")
    end
end