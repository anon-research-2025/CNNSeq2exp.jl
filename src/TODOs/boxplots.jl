
function plot_singleton_logos_matrix(logo_paths_mat, contribs_mat; 
        fig_size=(1000, 1200),
        columns_xrange = nothing,
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
            plot_image!(row_layout, logo_paths_mat[i, j]; occupy_till=img_end)    
    
            ax_text, ax_box, ax_scatter = prep_rain_axis(row_layout, ax_scatter; 
                ax_scatter_x_range = columns_xrange[j],
                occupy_start=rain_start, occupy_end=rain_end)
    
            plot_rain_clouds!(ax_text, ax_box, ax_scatter, contribs_mat[i, j]; 
                jitter_width=jitter_width,
                markersize=markersize,
                last_index=(i == rows),
                spine_as_well=true,
            )
            colgap!(row_layout, -5)  # Minimize column gap
            
            # top_bot_cond = (i == 1 && j == 1) || (i==rows && j == cols)

            Box(row_layout[1, 1:end], 
                strokewidth = 1, cornerradius = 2, color = (:tomato, 0.0), 
                strokecolor = :black)
        end
        colgap!(ga, 0)  # Minimize column gap
        rowgap!(ga, 5)  # Increase gap between rows for better spacing
    end
    
    colgap!(fig.layout, 5)  # Space between columns
    return fig
end

# rs, cs = 5, 3
# logo_paths_mat = reshape(logo_paths_here_sorted, (rs,cs))
# contribs_mat = reshape(contribs_here_sorted, (rs,cs))


# fig_size_test = (1100, 275)

# fig_size_test = (900, 275)

# columns_xrange = [(-0.25, 0.75), (-0.25, 0.25), (-0.5, 0.25)]

# fig = plot_singleton_logos_matrix(logo_paths_mat, contribs_mat; 
#         fig_size=fig_size_test .* 2,
#         columns_xrange = columns_xrange,
#         padding = 21,
#         jitter_width=0.25, 
#         markersize=0.6,
#         rain_start=6, 
#         rain_end=10,
#     );
# fig

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