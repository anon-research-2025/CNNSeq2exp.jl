
# load the paths and contribs


function plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=(1000, 1200),
        columns_xrange = nothing,
        img_start=1,
        img_end=5,
        padding = 1,
        jitter_width=0.1,
        markersize=0.4,
        rain_start=6, 
        rain_end=14,
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
                occupy_start=img_start, occupy_till=img_end)    

            ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
                ax_scatter_x_range= columns_xrange[j],
                occupy_start=rain_start, occupy_end=rain_end)

            plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_mat[i, j]; 
                jitter_width=jitter_width,
                markersize= markersize isa AbstractVector ? markersize[j] : markersize,
                last_index=(i == rows),
                spine_as_well=true,
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


# _m_.pfms[2][k]


logo_paths = String[]
contribs = Vector{float_type}[]

for k in keys_here
    d_sorted = sort(collect(keys(_m_.pfms[length(k)][k])))
    logo_paths_here = String[]
    contribs_here = Vector{float_type}[]
    for d in d_sorted
        _logo_paths_, _contribs_, _ = obtain_logo_paths_and_contribs(_m_, k, d, folder)
        push!(logo_paths_here, _logo_paths_[1])
        push!(contribs_here, _contribs_[1])
    end
    selected_indices = Int.([floor(quantile(1:length(d_sorted), i)) for i in 0:0.25:1])
    logo_paths_here = logo_paths_here[selected_indices]
    contribs_here = contribs_here[selected_indices]
    append!(logo_paths, logo_paths_here)
    append!(contribs, contribs_here)
end


##################


keys_here = [(m1 = "11", m2 = "11"), 
             (m1 = "11", m2 = "11", m3 = "11"), 
             (m1 = "1", m2 = "1")]

logo_paths_mat = reshape(logo_paths, (5,3))
contribs_mat = reshape(contribs, (5,3))

median.(contribs_mat)
minimum.(contribs_mat)
maximum.(eachcol(maximum.(contribs_mat)))

columns_xrange = [(-0.25, 1), (0.0, 0.75), (-0.25, 0.25)]


fig_size_test = (1100, 275)

fig = plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=fig_size_test .* 2,
        columns_xrange = columns_xrange,
        img_start=2,
        img_end=40,
        padding = 21,
        jitter_width=0.25, 
        markersize=[3, 6, 3],
        rain_start=41, 
        rain_end=90,
    );

fig


save("cross_results/pairs_and_triplets.png", fig; px_per_unit=2.5)



##########


keys_here = [(m1 = "9", m2 = "9"), 
             (m1 = "9", m2 = "9", m3 = "9"), 
             (m1 = "2", m2 = "2")]

logo_paths_mat = reshape(logo_paths, (5,3))
contribs_mat = reshape(contribs, (5,3))

columns_xrange = [(-0.25, 1), (0.0, 1.0), (-0.5, 0.25)]


fig_size_test = (1100, 275)

fig = plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=fig_size_test .* 2,
        columns_xrange = columns_xrange,
        img_start=2,
        img_end=40,
        padding = 21,
        jitter_width=0.25, 
        markersize=[3, 6, 3],
        rain_start=41, 
        rain_end=90,
    );

fig

save("yeast_result/pairs_and_triplets.png", fig; px_per_unit=2.5)
