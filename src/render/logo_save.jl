get_d_str(d::NamedTuple) = join([i for i in d], "_")
get_d_str(d::Integer) = string(d)
get_d_str(d::String) = d
get_pwm_singleton_save_path(k) = joinpath(x_dir[1], "$(k).png")
get_csv_singleton_save_path(k) = joinpath(x_dir[1], "$(k).csv")
get_box_singleton_save_path(k) = joinpath(x_dir[1], "box_$(k).png")

get_k_mode_str(k) = k isa String ? k : join(["$i" for i in k], "_")
get_pwm_save_folder(k; t = 1) = joinpath(x_dir[t], get_k_mode_str(k) )

get_mode_pwm_save_path(k; t = 1) = joinpath(get_pwm_save_folder(k; t = t), "mode.png")

get_mode_k(k) = k isa String ? 1 : length(k)
get_mode_gap_hist_save_path(k) = joinpath(get_pwm_save_folder(k; t = get_mode_k(k)), "gaps.png")
get_mode_boxplot_save_path(k) = joinpath(get_pwm_save_folder(k; t = get_mode_k(k)), "box_$(get_k_mode_str(k)).png")
get_mode_boxplot_save_path(k, d) = joinpath(get_pwm_save_folder(k; t = get_mode_k(k)), "box_$(get_d_str(d)).png")

get_pwm_save_path(k, d::Integer; t = 1) =
    joinpath(get_pwm_save_folder(k; t = t), "$(d).png")
get_csv_save_path(k, d::Integer; t = 1) =
    joinpath(get_pwm_save_folder(k; t = t), "$(d).csv")
get_pwm_save_path(k, d::NamedTuple; t = 1) =
    joinpath(get_pwm_save_folder(k; t = t), "$(get_d_str(d)).png")
get_csv_save_path(k, d::NamedTuple; t = 1) =
    joinpath(get_pwm_save_folder(k; t = t), "$(get_d_str(d)).csv")

# get_filter_index(k_m::String) = parse(Int, k_m[1]);
is_reversed_complement(k_m::String) = k_m[end] == 'r'

function make_highlighted_regions(k, d::Integer, hp)
    # get the expansion (left decrement, right increment)
    f1_endpos = hp.pfm_len
    f2_startpos = f1_endpos + d + 1
    f2_endpos = f2_startpos + hp.pfm_len - 1
    first_region = 1:f1_endpos
    second_region = f2_startpos:f2_endpos
    return [first_region, second_region]
end

function make_highlighted_regions(k, d::NamedTuple, hp)
    highlighted_regions = Vector{UnitRange{Int}}()
    next_start_pos = 1
    for (ind, f) in enumerate(k)
        f_endpos = next_start_pos + hp.pfm_len - 1
        push!(highlighted_regions, next_start_pos:f_endpos)
        ind != length(k) && (next_start_pos = f_endpos + d[ind] + 1)
    end
    return highlighted_regions
end

normalize_countmat(countmat) = countmat ./ sum(countmat, dims = 1)

function plot_horizontal_bars(ys::AbstractVector{<:Number}, 
                              y_names::Union{AbstractVector{String}, AbstractVector{<:Number}}, 
                              savepath::String;
                              plot_height_min=200,
                              plot_height_max=1200,
                              bar_offset_factor = 0.025
                              )
    # Ensure the lengths of ys and y_names match
    if length(ys) != length(y_names)
        error("Lengths of `ys` and `y_names` must match.")
    end

    # Determine dynamic text size
    n_bars = length(ys)
    plot_height = clamp(n_bars * 30, plot_height_min, plot_height_max)  # Adjust plot height dynamically (min: 200px, max: 1200px)
    text_size = PlotPWM.Plots.clamp(Int(ceil(plot_height / n_bars / 15)), 7, 20)  # Scale text size, limit between 7 and 20
    ylabel_size = clamp(Int(ceil(plot_height / 100)), 8, 20)  # Scale ylabel size dynamically based on plot height

    # Define the plot
    p = PlotPWM.Plots.bar(1:n_bars, ys,
        orientation = :horizontal,
        legend = false, 
        yrange = (0.5, n_bars + 0.5),  # Ensure the range fits the bars
        palette = repeat([:skyblue], n_bars),
        xticks = false,       # Turn off x-axis ticks
        yticks = (1:n_bars, y_names),  # Assign positions and names to y-axis ticks
        ylabel = "Gap size (number of nucleotides)",     # Y-axis title
        yguidefontsize = ylabel_size, # Adjust the font size of the ylabel dynamically
        margin = 10PlotPWM.Plots.mm,   # Add space for annotations
        size = (650, plot_height)  # Dynamically set the height of the plot
    )

    # Annotate the quantities on the bars
    annotation_offset = maximum(ys) * bar_offset_factor
    for (ind, yi) in enumerate(ys)
        text_here = PlotPWM.Plots.text("$yi", :black, text_size, :left)
        # @info "$text_here"
        PlotPWM.Plots.annotate!(p, yi + annotation_offset, ind, text_here)  # Offset for annotations
    end

    # Save the plot to the specified path
    PlotPWM.Plots.savefig(p, savepath)
end

function save_as_meme(this_pfm, save_name::String; name = "", num_sites = 100)
    # write as meme file as well
    io = open(save_name, "w")
    print(io, "MEME version 4\n\n")
    print(io, "ALPHABET= ACGT\n\n")
    print(io, "strands: + -\n\n")
    print(io, "Background letter frequencies\n")
    print(io, "A 0.25 C 0.25 G 0.25 T 0.25\n\n")
    print(io, "MOTIF $name \n")
    print(
        io,
        "letter-probability matrix: alength= 4 w= $(size(this_pfm,2)) nsites= $num_sites E= 0\n",
    )
    for i in axes(this_pfm, 2)
        print(io, " ")
        for j = 1:4
            print(io, "$(this_pfm[j,i]) ")
        end
        print(io, "\n")
    end
    close(io)
end


function save_pos_info_as_csv(pos_info, save_path)
    df = DataFrame(pos_info, csv_header)
    CSV.write(save_path, df)
end

function save_boxplots(_m_, saved_results_dir; t=2)
    x_motifs = _m_.pfms[t]
    for k in keys(x_motifs)
        # global plot
        logo_paths, contribs, parent_save_path  = obtain_logo_paths_and_contribs(_m_, k, saved_results_dir)
        # @info "k: $k"
        fig, ax_lefts, ax_rights, grid = setup_fig_and_axis(length(logo_paths))
        setup_plots!(ax_lefts, ax_rights, grid, logo_paths, contribs);
        name_here = get_k_mode_str(k)
        save(joinpath(parent_save_path, "box_$(name_here).png"), fig)

        if t > 1
            for d in keys(x_motifs[k])
                # local plot
                logo_paths, contribs, parent_save_path  = obtain_logo_paths_and_contribs(_m_, k, d, saved_results_dir)
                fig, ax_lefts, ax_rights, grid = setup_fig_and_axis(length(logo_paths))
                setup_plots!(ax_lefts, ax_rights, grid, logo_paths, contribs);
                save(joinpath(parent_save_path, "box_$(get_d_str(d)).png"), fig)
            end
        end
    end
end

function save_pfm_as_transfac(pfm, fp::String; count_default=1000)
    io = open(fp, "w")
    println(io, "ID\t")
    println(io, "XX\t")
    println(io, "BF\t")
    println(io, "XX\t")
    println(io, "P0\tA\tC\tG\tT")
    q = Int.(floor.(pfm .* count_default)); # make it a count matrix
    for j in axes(q, 2)
        cur_rows = j < 10 ? string(Int(floor(j)))*"$j" : string(j);
        println(io, cur_rows*"\t$(q[1,j])\t$(q[2,j])\t$(q[3,j])\t$(q[4,j])")
    end
    println(io, "XX\t")
    close(io)
end

#=
t = 1: pair_motifs
t = 2: triplet_motifs
... and so on
=#
function save_logos!(
    _m_::motifs, hp; sorted_configs = nothing,
    dpi = pwm_dpi,
    saved_results_dir = saved_results_dir,
    t = 2,
)
    x_motifs = _m_.pfms[t]
    pos_info = _m_.pos[t]
    ispos_contrib = _m_.is_pos_contrib[t]

    keys2iter = isnothing(sorted_configs) ? keys(x_motifs) : sorted_configs
    for k in keys2iter
        save_folder_here = joinpath(saved_results_dir, get_pwm_save_folder(k; t = t))
        mkpath(save_folder_here)

        gap_names = values.(keys(x_motifs[k]))
        gap_counts = map(x->sum((@view x[:,1]), dims=1)[1] |> Int, values(x_motifs[k]))
        gap_names_sorted_ind = sortperm(gap_names)
        gap_names_sorted = string.(@view gap_names[gap_names_sorted_ind])
        gap_counts_sorted = @view gap_counts[gap_names_sorted_ind]

        # @info "k: $k"
        plot_horizontal_bars(gap_counts_sorted, gap_names_sorted, joinpath(save_folder_here, "gaps.png"))
        
        for d in keys(x_motifs[k]) # k is the motifs tuple, d is the distance type
            pfm = normalize_countmat(x_motifs[k][d])
            save_as_meme(pfm, joinpath(save_folder_here, "$(get_d_str(d)).meme"))
            highlighted_regions = make_highlighted_regions(k, d, hp)
            PlotPWM.chk_overlap(highlighted_regions) && (@info "k: $k, d:$d")
            save_path_pfm = joinpath(save_folder_here, "$(get_d_str(d)).png")
            save_logoplot(
                pfm,
                save_path_pfm;
                dpi = dpi,
                highlighted_regions = highlighted_regions, 
                uniform_color = false, pos = ispos_contrib[k][d]
            )

            ######### for publication display only ##########
            # save as transfac file
            fp = joinpath(save_folder_here, "$(get_d_str(d)).transfac")
            save_pfm_as_transfac(pfm, fp)
            # pwm_path = joinpath(save_folder_here, "$(get_d_str(d))_.png")
            # save web_logo
            # Base.run(`weblogo 
            #     -D transfac 
            #     -f $fp 
            #     -n 100 
            #     -F png_print 
            #     -s large 
            #     --title " "
            #     --title-fontsize 12
            #     --fineprint ""
            #     --errorbars NO 
            #     --resolution 800 
            #     --fontsize 21 
            #     --number-fontsize 14
            #     --color-scheme classic 
            #     --aspect-ratio 3.75
            #     --stack-width 24
            #     --number-interval 60 
            #     -o $pwm_path`
            #     );
            ######### end of for publication display only ##########
            # save pos_info
            save_pos_info_as_csv(pos_info[k][d], joinpath(save_folder_here, "$(get_d_str(d)).csv"))
        end
    end
end

function save_logos_singleton(
    m::motifs;
    dpi = pwm_dpi,
    saved_results_dir = saved_results_dir,
    contrib_threshold = 0.01,
    transparent = 0.4,
    alpha = 1.0 
)
    singleton_motifs = m.pfms[1]
    pos_info = m.pos[1]
    _contribtions_ = m.pos_contribs[1]

    save_folder_here = joinpath(saved_results_dir, x_dir[1])
    mkpath(save_folder_here)

    for k in keys(singleton_motifs)
        pfm = normalize_countmat(singleton_motifs[k])
        save_path_pfm = joinpath(save_folder_here, "$k.png")
        mean_contribution = mean(_contribtions_[k])
        # alpha = abs(mean_contribution) > contrib_threshold ? 1.0 : transparent
        pos = mean_contribution > 0 ? true : false
        save_logoplot(pfm, save_path_pfm; dpi = dpi, alpha = alpha, uniform_color=false, pos=pos)

        ######### for publication display only ##########
        # save as transfac file
        fp = joinpath(save_folder_here, "$(k).transfac")
        save_pfm_as_transfac(pfm, fp)
        # pwm_path = joinpath(save_folder_here, "$(k)_.png")
        # save web_logo
    #    Base.run(`weblogo 
    #             -D transfac 
    #             -f $fp 
    #             -n 100 
    #             -F png_print 
    #             -s large 
    #             --title " "
    #             --title-fontsize 12
    #             --fineprint ""
    #             --errorbars NO 
    #             --resolution 800 
    #             --fontsize 21 
    #             --number-fontsize 14
    #             --color-scheme classic 
    #             --aspect-ratio 3.75
    #             --stack-width 24
    #             --number-interval 60 
    #             -o $pwm_path`
    #             );
        ######### end of for publication display only ##########

        save_path_csv = joinpath(saved_results_dir, x_dir[1], "$k.csv")
        save_pos_info_as_csv(pos_info[k], save_path_csv)
        save_as_meme(
            pfm, joinpath(saved_results_dir, x_dir[1], "$(get_d_str(k)).meme"),
        )    
    end
end