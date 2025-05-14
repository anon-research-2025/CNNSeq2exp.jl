function make_cross_color_df(crosslink_signitures::Vector{String})
    list = []
    for (ind, s) in enumerate(crosslink_signitures)
        push!(list, (crosslink_signiture = s, color = EntroPlots.crosslink_palette["$ind"]))
    end
    return DataFrame(list)
end

function render_jth_html(
    json_x_motifs,
    html_dict,
    real_upto;
    saved_results_dir = saved_results_dir,
    sequence_paths::Vector{String},
    singleton = false,
    protein_name = "",
    meta_str = "",
    j = 1
)
    # save the JSON file
    open(joinpath(saved_results_dir, "data$j.json"), "w") do io
        JSON3.pretty(io, json_x_motifs)
    end
    # copy the CSS file
    if !isfile(joinpath(saved_results_dir, "styles.css"))
        out = Mustache.render(template_css)
        io = open(joinpath(saved_results_dir, "styles.css"), "w")
        print(io, out)
        close(io)
    end
    # render the HTML file
    df = html_dict |> DataFrame

    if singleton
        html_template_rendered = Mustache.render(
            html_template_singleton,
            DF = df,
            protein_name = protein_name,
            j = j,
        )
    else
        html_template_rendered = Mustache.render(
            html_template,
            DF = df,
            protein_name = protein_name,
            j = j,
        )
    end

    # hover_on = "class=\"hover-window\" id=\"hoverWindow\""
    # hover = Mustache.render(html_hover_default,
    #         hover_on = hover_on,
    #         meta_str = meta_str
    # )

    # add hover and </body></html>
    html_template_rendered = html_template_rendered * html_end

    script_template_rendered = Mustache.render(
        script_template,
        mode_counts = size(df, 1),
        j = j,
        upto = real_upto,
        sequence_file_paths = sequence_paths,
    )

    io = open(saved_results_dir * "/index$j.html", "w")
    print(io, html_template_rendered)
    close(io)
    io = open(saved_results_dir * "/scripts$j.js", "w")
    print(io, script_template_rendered)
    close(io)
end

#=
Obtain the maximum number of motifs to display
    This is done to save time and space when rendering the motifs

    Args:
        x_motifs_test: Dict{NamedTuple, Dict{String, Array{float_type, 2}}}
            The motifs to display
    returns:
        sorted_configs: Array{NamedTuple}
            The sorted configurations 
            i.e. keys of x_motifs_test sorted by the average number of motifs
            (counts from high to low)

    note: max_display_motifs_x (integer) is a constant set in src/const.jl
=#
function obtain_display_list(m::motifs, x::Int)
    sorted_configs = sort(collect(m.pfms_mode_pvals[x]), by = x -> x[2])
    sorted_configs =
        [i[1] for i in sorted_configs[1:min(length(sorted_configs), max_display_motifs_x)]]
    # filter step: make sure that configs are in the test set
    filter!(y -> (y in keys(m.pfms_test[x])), sorted_configs)
    return sorted_configs
end

# _m_.pos_contribs[2]

function obtain_mode_avg_contributions(_m_::motifs, mode_i::Int)
    k_pos_avg_contribs = Dict{eltype(keys(_m_.pos_contribs[mode_i])), float_type}()
    for k in keys(_m_.pos_contribs[mode_i])
        mean_here = Vector{float_type}(undef, (length(_m_.pos_contribs[mode_i][k])))
        for (index, d) in enumerate(keys(_m_.pos_contribs[mode_i][k]))
            mean_here[index] = StatsBase.mean(_m_.pos_contribs[mode_i][k][d])
        end
        k_pos_avg_contribs[k] = StatsBase.mean(mean_here)
    end
    k_neg_avg_contribs = Dict{eltype(keys(_m_.neg_contribs[mode_i])), float_type}()
    for k in keys(_m_.neg_contribs[mode_i])
        mean_here = Vector{float_type}(undef, (length(_m_.neg_contribs[mode_i][k])))
        for (index, d) in enumerate(keys(_m_.neg_contribs[mode_i][k]))
            mean_here[index] = StatsBase.mean(_m_.neg_contribs[mode_i][k][d])
        end
        k_neg_avg_contribs[k] = StatsBase.mean(mean_here)
    end
    
end


# get the number of x-plet motifs to display
function get_real_upto(m::motifs)
    real_upto = 0
    for i in eachindex(m.pfms)
        if !isempty(m.pfms[i])
            real_upto += 1
        end
    end
    return real_upto
end

function get_sorted_config(_m_; t=2)
    if t == 1
        sorted_indices = sortperm(_m_.avg_contribs[t] |> values |> collect, by=x->abs(x), rev=true)
        sorted_configs = collect(keys(_m_.avg_contribs[t]))[sorted_indices]
        return sorted_configs
    else
        avg_contribs_i = Dict(k=>float_type(0) for k in keys(_m_.avg_contribs[t]))
        for k in keys(_m_.avg_contribs[t])
            place_holder = Vector{float_type}(undef, length(_m_.avg_contribs[t][k]))
            for (index, d) in enumerate(keys(_m_.avg_contribs[t][k]))
                place_holder[index] = _m_.avg_contribs[t][k][d]
            end
            avg_contribs_i[k] = StatsBase.mean(place_holder)
        end
        sorted_indices = sortperm(collect(values(avg_contribs_i)), by=x->abs(x), rev=true)
        sorted_configs = collect(keys(avg_contribs_i))[sorted_indices]
        return sorted_configs
    end
end


#=
save_path_where::String where to save the rendered HTML and logo files
t = 1 for pair motifs, t = 2 for triplet motifs, and so on
j = 1 for the first HTML file, j = 2 for the second HTML file, and so on
=#
function render_sorted_list(
    _m_::motifs, hp,
    save_path_where::String;
    sequence_paths::Vector{String}=["seqs.fa"],
    t = 2,
)
    # sorted_configs = obtain_display_list(m, t + 1)
    # sorted_configs = nothing
    sorted_configs = get_sorted_config(_m_; t=t)

    make_save_paths(saved_results_dir = save_path_where)

    save_logos!(
        _m_, hp;
        sorted_configs = sorted_configs,
        t = t, 
        saved_results_dir = save_path_where
    )

    save_boxplots(_m_, save_path_where; t=t)

    json_x_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()

    update_json_and_dict!(
        _m_,
        json_x_motifs,
        html_dict;
        t = t,
        sorted_configs = sorted_configs,
    )

    render_jth_html(
        json_x_motifs,
        html_dict,
        get_real_upto(_m_);
        j = t,
        saved_results_dir = save_path_where,
        sequence_paths = sequence_paths,
    )
end

function render_singletons(
    _m_::motifs,
    save_path_where::String;
    sequence_paths::Vector{String}=["seqs.fa"],
    rna = false,
)
    # sorted_configs = [i[1] for i in sort(pvals |> collect, by = x -> x[2])]
    # sorted_configs = nothing
    sorted_configs = get_sorted_config(_m_; t=1)
    save_logos_singleton(
        _m_;
        saved_results_dir = save_path_where,
    )

    save_boxplots(_m_, save_path_where; t=1)

    json_x_motifs = init_json_dict()
    html_dict = init_dict_for_html_render()
    update_json_and_dict!(
        _m_,
        json_x_motifs,
        html_dict;
        sorted_configs = sorted_configs,
        t = 1,
        rna = rna,
    )
    # println(html_dict);
    render_jth_html(
        json_x_motifs,
        html_dict,
        get_real_upto(_m_);
        j = 1,
        saved_results_dir = save_path_where,
        singleton = true,
        sequence_paths = sequence_paths,
    )
end

# copy the readme files to save_path_where
function copy_readme(save_path_where)
    source_dir = joinpath(Base.pathof(SparseCodeMotifs) |> dirname, "render/readme")
    @assert isdir(source_dir) "The source path does not exist"
    dest_dir = joinpath(save_path_where, "readme")
    mkpath(dest_dir)
    for file in readdir(source_dir)
        src_file = joinpath(source_dir, file)
        dest_file = joinpath(dest_dir, file)
        if isfile(dest_file)
        else
            cp(src_file, dest_file)
        end
    end
end