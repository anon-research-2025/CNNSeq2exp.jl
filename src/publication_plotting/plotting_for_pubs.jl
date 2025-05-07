
function plot_image!(row_layout, i_logo_paht; occupy_start=1, occupy_till=5, 
        margin=nothing)
    ax_logo = Axis(row_layout[1, occupy_start:occupy_till], title="")
    hidespines!(ax_logo)
    hidedecorations!(ax_logo)
    logo = load(i_logo_paht) |> rotr90
    image!(ax_logo, logo)
    if !isnothing(margin)
        ax_logo.xautolimitmargin = (margin, margin)
    end
end


# Create a function to remove trailing zeros, dots after 0/1, and the leading zero in fractions
function format_ticks(value::Real)
    # Format the number to a maximum of 3 decimal places
    formatted = @sprintf("%.2f", value)
    formatted = chomp(formatted)  # Remove trailing newline

    # Remove the dot after 0 and 1 (if value is whole)
    if formatted == "0.00"
        return "0" 
    elseif formatted == "1.00"
        return "1"
    elseif formatted == "-1.00"
        return "-1"
    else 
        return ""
    end
    # @info "formatted: $formatted"
    # Remove trailing zeros after the decimal point
    formatted = replace(formatted, r"0*$" => "")  # Remove trailing zeros

    # Remove leading zero in fractions
    if startswith(formatted, "0.")
        formatted = replace(formatted, r"^0" => "")
    elseif startswith(formatted, "-0")
        formatted = replace(formatted, r"^-0" => "-")
    end
    return formatted
end

function prep_rain_axis(row_layout, ax_scatter_prev; 
    ax_scatter_x_range=(-0.25, 0.5), 
    occupy_start=6, 
    occupy_end=14,
    )
    rain_grid = row_layout[1, occupy_start:occupy_end] = GridLayout()
    ax_text = Axis(rain_grid[1:5, 1], limits = (nothing, (0, 0.5)))
    ax_box = Axis(rain_grid[5:6, 1])
    ax_scatter = Axis(rain_grid[7:9, 1], xticklabelsize=18,
        limits = (ax_scatter_x_range, nothing))

    ax_scatter.xtickformat = values -> [format_ticks(value) for value in values]
    linkaxes!(ax_text, ax_box, ax_scatter)
    if !isnothing(ax_scatter_prev)
        linkaxes!(ax_scatter_prev, ax_text)
    end
    rowgap!(rain_grid, 4)
    return ax_text, ax_box, ax_scatter
end

function plot_rain_clouds!(ax_text, ax_box, ax_scatter, pts_i; 
    last_index=false, 
    spine_as_well=false,
    jitter_width=0.1,
    markersize = 0.4,
    lw = 1,
    text_size=21,
    scatter_spinewidth=2.5,
    )
    # text
    med = round(StatsBase.median(pts_i), digits=2)
    text!(ax_text, med .+ 0.1, -0.8, text = "$med", 
        align = (:right, :center), fontsize = text_size, color = :black)
    # box plot
    Makie.boxplot!(ax_box, fill(-0.8, length(pts_i)), pts_i; 
        orientation=:horizontal,
        color = (:orange, 0.98),
        whiskercolor = :black, 
        mediancolor = :black, 
        strokecolor = :black, 
        whiskerwidth = 1,  # width of the T (top)
        whiskerlinewidth = lw, # width of the T (bottom)
        strokewidth = lw, 
        medianlinewidth = lw,
        gap = 0.5,
        show_outliers=false,
        label = "horizontal");
    # scatter
    CairoMakie.scatter!(ax_scatter, 
        pts_i, jitter_width .* randn(length(pts_i)), 
            color = :maroon, 
            alpha=0.8,
            # markersize = 2, # for singletons
            markersize = markersize, # for pairs
            # markersize = 4,
            fxaa = true,
            marker=:circle
            );
    # spine
    hide_everything!(ax_text; spine_as_well=spine_as_well) # hide x-axis as well
    hide_everything!(ax_box; spine_as_well=spine_as_well) # hide x-axis as well
    if !last_index
        hide_everything!(ax_scatter; spine_as_well=spine_as_well) # hide x-axis as well
    else
        hide_everything_except_xaxis!(ax_scatter; spine_as_well=spine_as_well)
        ax_scatter.spinewidth = scatter_spinewidth
        ax_scatter.xticksize = 4
        ax_scatter.xtickwidth = 2.5
    end
end

function plot_singleton_logos(logo_paths, contribs_each; 
        fig_size=(1000, 345),
        padding = 1,
        img_start=1,
        img_end=5,
        rain_start=6, 
        rain_end=14,
        jitter_width=0.1,
        markersize=0.4,
        )
    fig = Figure(size=fig_size, figure_padding = padding);
    ga = fig[1, 1] = GridLayout()  # GridLayout for precise control
    
    @assert length(logo_paths) == length(contribs_each) "Length mismatch"

    ax_scatter = nothing
    for i in axes(contribs_each, 1)
        row_layout = ga[i, 1] = GridLayout()
        plot_image!(row_layout, logo_paths[i]; 
            occupy_start=img_start, occupy_till=img_end)    

        ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
            ax_scatter_x_range=(-1, 1),
            occupy_start=rain_start, occupy_end=rain_end)

        plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_each[i]; 
                last_index=(i == length(contribs_each)),
                jitter_width=jitter_width,
                markersize=markersize,
                spine_as_well=true,
                )

        Box(row_layout[1, 1:end], 
                strokewidth = 1, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :gray)

    end
    colgap!(ga, 1)  # Minimize column gap
    rowgap!(ga, 5)  # Increase gap between rows for better spacing
    return fig
end


function plot_singleton_logos_matrix(logo_paths_mat, contribs_mat, solid_indices; 
        fig_size=(1000, 1200),
        img_end=5,
        padding = 1,
        jitter_width=0.1,
        markersize=0.4,
        rain_start=6, 
        rain_end=14,
        c_and_r_gap=10,
        box_stroke_width=1,
        x_margin_img=0.05,
        text_size=16,
        scatter_spinewidth=4.0,
        )
    @assert size(logo_paths_mat) == size(contribs_mat) "Dimension mismatch"
    
    rows, cols = size(logo_paths_mat)
    fig = Figure(size=fig_size, figure_padding = padding);
    grid_layouts = [fig[1, j] = GridLayout() for j in 1:cols]  # One GridLayout per column
    
    for j in 1:cols
        ga = grid_layouts[j]  # Current column layout
        
        ax_scatter = nothing
        for i in 1:rows

            row_layout = ga[i, 1] = GridLayout()
            plot_image!(row_layout, logo_paths_mat[i, j]; 
                occupy_till=img_end, margin=x_margin_img, 
                )
    
            ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
                ax_scatter_x_range = (-1,1),
                occupy_start=rain_start, occupy_end=rain_end)
            ax_scatter.xticklabelsize = 32

            plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_mat[i, j]; 
                jitter_width=jitter_width,
                markersize=markersize,
                last_index=(i == rows),
                spine_as_well=true,
                text_size=text_size,
                scatter_spinewidth=scatter_spinewidth,
            )
            colgap!(row_layout, -5)  # Minimize column gap
                        
            # top_bot_cond = (i == 1 && j == 1) || (i==rows && j == cols)

            solid_cond = (i,j) in solid_indices

            Box(row_layout[1, 1:end], 
                strokewidth = box_stroke_width, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :black, 
                linestyle = solid_cond ? :solid : :dot,  # Solid for top and bottom, dotted for others
                )
        end
        colgap!(ga, 0)  # Minimize column gap
        rowgap!(ga, c_and_r_gap)  # Increase gap between rows for better spacing
    end
    
    colgap!(fig.layout, c_and_r_gap)  # Space between columns
    return fig
end

function plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=(1400, 1200),
        img_start=1,
        img_end=5,
        padding = 17,
        jitter_width=0.1,
        markersize=2.4,
        rain_start=6, 
        rain_end=8,
        x_margin_img=0.05,
        text_size=16,
    )
    @assert size(logo_paths_mat) == size(contribs_mat) "Dimension mismatch"

    rows, cols = size(logo_paths_mat)
    fig = Figure(size=fig_size, figure_padding = padding);
    grid_layouts = [fig[1, j] = GridLayout() for j in 1:cols]  # One GridLayout per column

    for j in 1:cols
        ga = grid_layouts[j]  # Current column layout
        
        ax_scatter = nothing
        
        for i in 1:rows
            row_layout = ga[i, 1] = GridLayout()
            plot_image!(row_layout, logo_paths_mat[i, j]; 
                occupy_start=img_start, occupy_till=img_end, margin=x_margin_img)

            ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
                ax_scatter_x_range= (-1,1),
                occupy_start=rain_start, occupy_end=rain_end)

            plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_mat[i, j]; 
                jitter_width=jitter_width,
                markersize= markersize isa AbstractVector ? markersize[j] : markersize,
                last_index=(i == rows),
                spine_as_well=true,
                text_size=text_size
            )
            colgap!(row_layout, -5)  # Minimize column gap
            
            # top_bot_cond = (i == 1 && j == 1) || (i==rows && j == cols)

            Box(row_layout[1, 1:end], 
                strokewidth = 1, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :black)
        end
        colgap!(ga, 25)  # Minimize column gap
        rowgap!(ga, 5)  # Increase gap between rows for better spacing
    end

    colgap!(fig.layout, 5)  # Space between columns
    return fig
end


function sort_by_median(contribs_here, logo_paths_here)
    contribs_median = median.(contribs_here)
    sortperm_contribs = sortperm(contribs_median, rev=true)
    logo_paths_here_sorted = logo_paths_here[sortperm_contribs]
    contribs_here_sorted = contribs_here[sortperm_contribs]
    contribs_median_sorted = contribs_median[sortperm_contribs]
    return logo_paths_here_sorted, contribs_here_sorted, contribs_median_sorted
end

function get_singleotons_paths_sorted(_m_, folder)
    ### do the plot for all the singletons
    logo_paths_here = String[]
    contribs_here = Vector{float_type}[]
    # for i = 1:15
    for i in keys(_m_.pfms[1])
        logo_paths, contribs, parent_save_path = 
            obtain_logo_paths_and_contribs(_m_, i, folder)
        push!(logo_paths_here, logo_paths[1])
        push!(contribs_here, contribs[1])
    end

    logo_paths_here_sorted, contribs_here_sorted, contribs_median_sorted = 
        sort_by_median(contribs_here, logo_paths_here);
    
    return logo_paths_here_sorted, contribs_here_sorted, contribs_median_sorted
end

function get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=5)
    logo_paths = String[]
    contribs = Vector{float_type}[]

    for k in keys_here
        d_sorted = sort(collect(keys(_m_.pfms[length(k)][k])))
        logo_paths_here = String[]
        contribs_here = Vector{float_type}[]
        for d in d_sorted
            _logo_paths_, _contribs_, _ = 
                obtain_logo_paths_and_contribs(_m_, k, d, folder)
            push!(logo_paths_here, _logo_paths_[1])
            push!(contribs_here, _contribs_[1])
        end
        selected_indices = Int.([floor(quantile(1:length(d_sorted), i)) for i in 0:(1/(take_how_many-1)):1])
        logo_paths_here = logo_paths_here[selected_indices]
        contribs_here = contribs_here[selected_indices]
        append!(logo_paths, logo_paths_here)
        append!(contribs, contribs_here)
    end
    return logo_paths, contribs
end

##################################################################################


function get_logos_left_done!(fig, single_row_range, logo_paths, contribs_each, 
        jitter_width,
        markersize_singletons,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end, 
        box_spinewidth,
    )
    ga = fig[single_row_range, 1:2] = GridLayout()  # GridLayout for precise control    
    @assert length(logo_paths) == length(contribs_each) "Length mismatch"
    ax_scatter = nothing
    for i in axes(contribs_each, 1)
        row_layout_0 = ga[i, 1:6] = GridLayout()
        ax_text_0 = Axis(row_layout_0[1,1]); hide_everything!(ax_text_0)
        # place text
        text!(ax_text_0, 1, 1, text = "filter $i", 
            align = (:center, :center), fontsize = 24, color = :black, 
            rotation = pi/2)

        row_layout = ga[i, 8:99] = GridLayout()
        plot_image!(row_layout, logo_paths[i]; 
            occupy_start=img_start, occupy_till=img_end,  margin=x_margin_img) 

        ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
            ax_scatter_x_range=(-1,1),
            occupy_start=rain_start, occupy_end=rain_end)
        ax_scatter.xticklabelsize = 28

        plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_each[i]; 
                last_index=(i == length(contribs_each)),
                jitter_width=jitter_width,
                markersize=markersize_singletons,
                spine_as_well=true,
                )

        Box(row_layout[1, 1:end], 
                strokewidth = box_spinewidth, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :black)

    end
    colgap!(ga, 1)  # Minimize column gap
    rowgap!(ga, 10)  # Increase gap between rows for better spacing

    # empty axis
    empty_ax = Axis(ga[:, 100], title="")
    empty_ax.spinewidth = 0
    hidedecorations!(empty_ax)
end

function get_logos_right_done!(fig, logo_paths_mat, contribs_mat, 
        jitter_width,
        markersize_plets,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end, 
        box_spinewidth,
        col_names_here;
        ax_scatter_x_range=(-1,1),
    )

    rows_xplets, _ = size(logo_paths_mat)
    # pair, triplets, ...    
    grid_layouts = [fig[1:6, (3*j):(3*j+2)] = GridLayout() 
        for j in axes(contribs_mat,2)]  # One GridLayout per column

    for j in axes(contribs_mat,2)
        gb = grid_layouts[j]  # Current column layout
        Label(gb[:, :, Top()], col_names_here[j], valign = :bottom,
            fontsize= 24, color = :black,
            padding = (0, 0, 5, 0))
        ax_scatter = nothing
        for i in 1:rows_xplets
            row_layout = gb[i, 1] = GridLayout()
            plot_image!(row_layout, logo_paths_mat[i, j]; 
                occupy_start=img_start, occupy_till=img_end, margin=x_margin_img) 

            ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
                ax_scatter_x_range= ax_scatter_x_range,
                occupy_start=rain_start, occupy_end=rain_end)
            ax_scatter.xticklabelsize = 28
            
            plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_mat[i, j]; 
                jitter_width=jitter_width,
                markersize= markersize_plets isa AbstractVector ? markersize_plets[j] : markersize_plets,
                last_index=(i == rows_xplets),
                spine_as_well=true,
            )
            colgap!(row_layout, -5)  # Minimize column gap
            
            Box(row_layout[1, 1:end], 
                strokewidth = box_spinewidth, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :black)
        end
        colgap!(gb, 25)  # Minimize column gap
        rowgap!(gb, 10)  # Increase gap between rows for better spacing
        @info "here $j done !"
    end
    @info "here 2 !"
    colgap!(fig.layout, 10)  # Space between columns
end


function plot_singleton_logos_extras(
        logo_paths, contribs_each,
        logo_paths_mat,
        contribs_mat; 
        fig_size=(1000, 345),
        padding = 1,
        img_start=1,
        img_end=5,
        rain_start=6, 
        rain_end=14,
        jitter_width=0.1,
        markersize_singletons=0.4,
        markersize_plets=1.0,
        x_margin_img=1.0,
        single_row_range = 2:5,
        box_spinewidth=1.5,
        col_names_here=fill("", length(logo_paths_mat)),
        )
    fig = Figure(size=fig_size, figure_padding = padding);

    get_logos_left_done!(fig, single_row_range, logo_paths, contribs_each, 
        jitter_width,
        markersize_singletons,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end, 
        box_spinewidth,
    )

    get_logos_right_done!(fig, logo_paths_mat, contribs_mat, 
        jitter_width,
        markersize_plets,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end, 
        box_spinewidth,
        col_names_here,
    )

    return fig
end

struct mixed_plt
    l_selected::Vector{String}
    c_selected::Vector{Vector{float_type}}
    logo_paths_mat::Matrix{String}
    contribs_mat::Matrix{Vector{float_type}}
    col_names_here::Vector{String}
    function mixed_plt(_m_, folder, l_selected, c_selected, keys_here, col_names_here; 
            take_how_many=5, )
        logo_paths_mat, contribs_mat = 
            get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=take_how_many)
        logo_paths_mat = reshape(logo_paths_mat, (take_how_many,3))[1:3,:]
        contribs_mat = reshape(contribs_mat, (take_how_many,3))[1:3,:]
        return new(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
    end
    function mixed_plt(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
        return new(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
    end
end

function plot_singleton_logos_extras_0(mp::mixed_plt;
        fig_size=(1400, 275) .* 1.5,
        padding = 25,
        single_row_range = 2:5,
        markersize_singletons=0.4,
        markersize_plets=[3,4,3],
        jitter_width=0.1,
        img_start=1,
        img_end=5,
        rain_start=6, 
        rain_end=8,
        rain_end_single=11,
        x_margin_img=0.025,
        box_spinewidth=1.5,
        ax_scatter_x_range=(-1,1),
    )

    fig = Figure(size=fig_size, figure_padding = padding);

    get_logos_left_done!(fig, single_row_range, 
        mp.l_selected, mp.c_selected, 
        jitter_width,
        markersize_singletons,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end_single, 
        box_spinewidth,
    )

    get_logos_right_done!(fig, 
        mp.logo_paths_mat, 
        mp.contribs_mat, 
        jitter_width,
        markersize_plets,
        img_start, img_end, x_margin_img, 
        rain_start, rain_end, 
        box_spinewidth,
        col_names_here;
        ax_scatter_x_range=ax_scatter_x_range
    )
    return fig
end

function plot_singleton_logos_extras_1(mps::Vector{mixed_plt};
        fig_size=(1400, length(mps) * 275 + 10) .* 1.5,
        padding = 25,
        single_row_range = 2:5,
        markersize_singletons=0.4,
        markersize_plets=[3,4,3],
        jitter_width=0.1,
        img_start=1,
        img_end=5,
        rain_start=6, 
        rain_end=8,
        rain_end_single=11,
        x_margin_img=0.025,
        box_spinewidth=1.5,
        ABCs = ["A", "B", "C", "D", "E", "F", "G", "H"],
        ax_scatter_x_range=(-1,1),

    )

    fig = Figure(size=fig_size, figure_padding = padding);
    for (index, mp) in enumerate(mps)
        g = fig[index, 1] = GridLayout() 
        get_logos_left_done!(g[1,1], 
            single_row_range isa Vector ? single_row_range[index] : single_row_range, 
            mp.l_selected, mp.c_selected, 
            jitter_width,
            markersize_singletons,
            img_start, img_end, x_margin_img, 
            rain_start, rain_end_single, 
            box_spinewidth,
        )
        get_logos_right_done!(g[1,1], 
            mp.logo_paths_mat, 
            mp.contribs_mat, 
            jitter_width,
            markersize_plets,
            img_start, img_end, x_margin_img, 
            rain_start, rain_end, 
            box_spinewidth,
            mp.col_names_here;
            ax_scatter_x_range=ax_scatter_x_range
        )

        Label(g[1, 1, Top()], ABCs[index],
            fontsize = 36,
            font = :bold,
            padding = (0, 5, 5, 0),
            halign = :left)
    end
    rowgap!(fig.layout, 60)  # Increase gap between rows for better spacing
    return fig
end

#########################################################




# TODO make constant later
len2rs = Dict(16=>(4,4), 15=>(5,3), 12=>(4,3), 10=>(2,5));
rs2figsize = Dict((4,4)=>(1200, 350), (5,3)=>(1200, 400), (2,5)=>(1200, 200), (4,3)=>(1200, 300));

# setting 2
len2rs = Dict(16=>(4,4), 15=>(3,5), 12=>(4,3), 10=>(2,5), 9 =>(3,3));
rs2figsize = Dict((4,4)=>(1200, 350), (3,5) => (1200,225), 
    (5,3)=>(1200, 400), (2,5)=>(1200, 200), (4,3)=>(1200, 300), (3,3)=>(800, 265));

# Given number of rows m, define f(z)
function linear_to_rc(z, m)
    row = mod1(z, m)
    col = div(z - 1, m) + 1
    return (row, col)
end


function get_singletons_matrix_fig_inner(logo_paths_mat, contribs_mat, solid_indices; 
    fig_size_test=(1200, 300), )
    fig = plot_singleton_logos_matrix(logo_paths_mat, contribs_mat, solid_indices; 
            fig_size=fig_size_test .* 2,
            padding = 21,
            jitter_width=0.15, 
            markersize=0.6,
            rain_start=6, 
            rain_end=12,
            c_and_r_gap=10,
            box_stroke_width=2,
            text_size=29,
        );
    return fig
end

function get_singletons_matrix_fig(_m_, folder, solid_indices)
    # load the singletons
    logo_paths_here_sorted, contribs_here_sorted, _ = 
        get_singleotons_paths_sorted(_m_, folder)
    rs, cs  = len2rs[length(logo_paths_here_sorted)]

    # map solid indices to matrix indices
    solid_indices_use = linear_to_rc.(solid_indices, rs)

    logo_paths_mat = reshape(logo_paths_here_sorted, (rs,cs))
    contribs_mat = reshape(contribs_here_sorted, (rs,cs))
    fig = get_singletons_matrix_fig_inner(logo_paths_mat, contribs_mat, solid_indices_use; 
        fig_size_test=rs2figsize[(rs,cs)])
    return fig, logo_paths_here_sorted, contribs_here_sorted
end

function get_plets_fig(_m_, folder, keys_here; take_how_many=15, fig_size=(1400, 1200))
    logo_paths, contribs = 
        get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=take_how_many)
    # determine the shape
    logo_paths_mat = reshape(logo_paths, (take_how_many,3))
    contribs_mat = reshape(contribs, (take_how_many,3))
    fig = plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=fig_size,
        markersize=[3,4,3])
    return fig
end

