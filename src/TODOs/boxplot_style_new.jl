
using StatsBase

using DataFrames, PlotPWM, CSV
using Printf, JSON3, Mustache
using CairoMakie, ImageShow
using JLD2
using Flux

include("const.jl")

include("motif_analysis/motif.jl")
include("motif_analysis/helper.jl")
include("render/const.jl")
include("render/onehot2fasta.jl")
include("render/logo_save.jl")
include("render/fill_json.jl")
include("render/template_css.jl")
include("render/template_script.jl")
include("render/templates_html.jl")
include("render/boxplots.jl")
include("render/render.jl")

include("render/boxplot_publication.jl")



colors = [Makie.wong_colors(); Makie.wong_colors()];
# Modify the color with transparency (alpha = 0 for full transparency)
# transparent_colors = [RGBA(c.r, c.g, c.b, 0.0) for c in colors]
# Generate example data

function plot_image!(row_layout, i_logo_paht; occupy_start=1, occupy_till=5)
    ax_logo = Axis(row_layout[1, occupy_start:occupy_till], title="")
    hidespines!(ax_logo)
    hidedecorations!(ax_logo)
    logo = load(i_logo_paht) |> rotr90
    image!(ax_logo, logo)

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
    )
    # text
    med = round(StatsBase.median(pts_i), digits=2)
    text!(ax_text, med .+ 0.1, -0.8, text = "$med", 
        align = (:right, :center), fontsize = 21, color = :black)
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
        ax_scatter.spinewidth = 2.5
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
            ax_scatter_x_range=(-0.25, 0.5),
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


#################################################

@load "save_m/cross_m.jld2" _m_

function sort_by_median(contribs_here, logo_paths_here)
    contribs_median = median.(contribs_here)
    sortperm_contribs = sortperm(contribs_median, rev=true)
    logo_paths_here_sorted = logo_paths_here[sortperm_contribs]
    contribs_here_sorted = contribs_here[sortperm_contribs]
    contribs_median_sorted = contribs_median[sortperm_contribs]
    return logo_paths_here_sorted, contribs_here_sorted, contribs_median_sorted
end

display_size = (1050, 2400)
test_size = display_size ./ 2

folder = "../saved_results/cross_gpro_code_thresh_0.8_c_0.75_2__3"

folder = "../saved_results/aviv_yeast_code_thresh_0.8_c_0.75__3"

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



fig1 = plot_singleton_logos(logo_paths_here_sorted, contribs_here_sorted; 
    fig_size=test_size,
    );

fig1



save("try.png", fig1; px_per_unit=2)


contribs_here[10] |> histogram
