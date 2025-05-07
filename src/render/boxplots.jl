function hide_x_ticks!(ax::Axis)
    ax.xticksvisible = false
    ax.xticklabelsvisible = false
end

function hide_y_ticks!(ax::Axis)
    ax.yticksvisible = false
    ax.yticklabelsvisible = false
end

function hide_ticks!(ax::Axis)
    hide_x_ticks!(ax)
    hide_y_ticks!(ax)
end

function hide_xlabels!(ax::Axis)
    ax.xlabelvisible = false
end

function hide_ylabels!(ax::Axis)
    ax.ylabelvisible = false
end

function hide_xylabels!(ax::Axis)
    hide_xlabels!(ax)
    hide_ylabels!(ax)
end

function hide_spine!(ax::Axis)
    ax.spinewidth = 0
end

function change_x_tick_size!(ax::Axis, howbig)
    ax.xticklabelsize = howbig
end

function change_x_tick_len!(ax::Axis, howbig)
    ax.xticksize = howbig
end

function hide_ticks_and_spine!(ax::Axis)
    hide_ticks!(ax)
    hide_spine!(ax)
end

function hide_grid!(ax::Axis)
    ax.xgridwidth=0
    ax.ygridwidth=0
    ax.xgridstyle=:none 
    ax.ygridstyle=:none 
    ax.xgridvisible=0
    ax.ygridvisible=0
end

hide_everything_except_title!(ax::Axis) = 
    (hide_ticks_and_spine!(ax); hide_xylabels!(ax))

function get_quantiles(x::Int)
    q = range(0, 1, length=x)[2:end-1] |> collect
    return q
end

# get the quantile 0, .25, .5, .75, 1 of the index set
function get_quantile_indices(list; q=get_quantiles(5)) 
    @assert length(list) ≥ 5 "List is too short"
    return vcat([1], Int.([ceil(length(list)*q_here) for q_here in q]), [length(list)])
end

#=
return also the parent save path to save the individual box plot
=#
function obtain_logo_paths_and_contribs(_m_, k, save_path_where; max_display_num=5, min_count=0)
    t = typeof(k) == String ? 1 : length(k)
    if t == 1
        parent_save_path = joinpath(save_path_where, x_dir[1])
        logo_paths = [joinpath(parent_save_path, k*".png")]
        contribs = [_m_.pos_contribs[1][k]]
        return logo_paths, contribs, parent_save_path
    else
        parent_save_path = joinpath(save_path_where, get_pwm_save_folder(k, t=t))
        distances = sort!(collect(keys(_m_.pfms[t][k]))) # k is the motifs tuple, d is the distance type

        logo_paths = Vector{String}(undef, min(length(distances), max_display_num))
        contribs = Vector{Vector{float_type}}(undef, length(logo_paths))

        if length(distances) < max_display_num
            for (index, d) in enumerate(distances)
                logo_paths[index] = joinpath(parent_save_path, "$(get_d_str(d)).png")
                contribs[index] =  _m_.is_pos_contrib[t][k][d] ? _m_.pos_contribs[t][k][d] : _m_.neg_contribs[t][k][d]
            end
        else
            indices2use = get_quantile_indices(collect(1:length(distances)), q=get_quantiles(max_display_num))
            for (index, q_index) in enumerate(indices2use)
                d = distances[q_index]
                logo_paths[index] = joinpath(parent_save_path, "$(get_d_str(d)).png")
                contribs[index] = _m_.is_pos_contrib[t][k][d] ? _m_.pos_contribs[t][k][d] : _m_.neg_contribs[t][k][d]
            end
        end

        selected_indices = findall(length.(contribs) .> min_count)
        return (@view logo_paths[selected_indices]), (@view contribs[selected_indices]), parent_save_path
    end
end

function obtain_logo_paths_and_contribs(_m_, k, d, save_path_where)
    t = length(k)
    parent_save_path = joinpath(save_path_where, get_pwm_save_folder(k, t=t))
    logo_path = joinpath(parent_save_path, "$(get_d_str(d)).png")
    contrib = _m_.is_pos_contrib[t][k][d] ? _m_.pos_contribs[t][k][d] : _m_.neg_contribs[t][k][d]
    return [logo_path], [contrib], parent_save_path
end


function setup_fig_and_axis(how_many_rows; 
    row_width=175, 
    col_width=600)
    # Create figure
    # height_factor = how_many_rows < 4 ? 0.35 : 0.6

    H = row_width*how_many_rows

    # W = max(col_width, row_width*how_many_rows); H = height_factor * W
    title_gap_here = 0.05*row_width
    fig = Figure(size = (col_width, H));

    ax_lefts = Vector{Axis}(undef, how_many_rows)
    ax_rights = Vector{Axis}(undef, how_many_rows)
    grid = fig[1, 1] = GridLayout()
    grid_left = grid[1:how_many_rows, 1]
    grid_right = grid[1:how_many_rows, 2]

    # @info "how_many_rows: $how_many_rows"

    ax_left = Axis(grid_left, title = "Sequence logo", titlesize = 12)
    ax_right = Axis(grid_right, title = "Banzhaf index (expected change)", titlesize = 12)
    # ax_right = Axis(grid_right, title = "Bojo index (X-Y)", titlesize = 16)
    hide_grid!(ax_left); hide_grid!(ax_right)
    hidedecorations!(ax_left); hidedecorations!(ax_right);
    hide_everything_except_title!(ax_left); hide_everything_except_title!(ax_right);
    ax_left.titlegap = title_gap_here; 
    ax_right.titlegap = title_gap_here;
    return fig, ax_lefts, ax_rights, grid
end

function setup_plots!(ax_lefts, ax_rights, grid, logo_paths, contribs; 
    text_size=0.025*row_width)
    how_many_rows = length(logo_paths)
    # Define grid layout
    for r = 1:how_many_rows
        box_data = contribs[r]  # Box plot data
        img = ImageShow.load(logo_paths[r]);
        img = rotr90(img);  # Rotate 90 degrees counterclockwise

        median_box_data = StatsBase.median(box_data)
        x_limit_right = max(1.0, median_box_data+0.25) # so that the word median does not get cut off

        ax_lefts[r] = Axis(grid[r,1], aspect=DataAspect(), spinewidth=0)
        ax_rights[r] = Axis(grid[r,2], limits = ((-1.0, x_limit_right), (0, 1)), 
            spinewidth=0.75, xticksize=1.0)
        image!(ax_lefts[r], img); 
        hide_grid!(ax_lefts[r]); hide_grid!(ax_rights[r])

        CairoMakie.boxplot!(ax_rights[r], fill(0.65, length(box_data)), box_data, 
            width = 0.275, color = (:navy, 0.45),
            whiskerwidth = 0.5,  # Narrower whiskers
            whiskercolor = (:navy, 1), 
            range = 1.5, 
            mediancolor = :navy, markersize = 3, strokecolor = :black,            
            strokewidth = 1, label = "Horizontal", orientation = :horizontal)
        # Scatter beneath the box plot (making it slightly below the plot)
        med = round2(median_box_data)
        text!(ax_rights[r], med .+ 0.35, 0.88, text = "median $med", 
            align = (:right, :center), fontsize = 12, color = :black)
        # text!(ax_rights[r], med .+ 0.225, 0.85, text = "median $med", 
        #     align = (:right, :center), fontsize = 12, color = :black)

        CairoMakie.scatter!(ax_rights[r], 
            box_data, 0.035 .* randn(length(box_data)) .+ 0.35 , 
                color = :maroon, markersize = boxplot_marker_size)
            # box_data, 0.035 .* randn(length(box_data)) .+ 0.35 , color = :maroon, markersize = 5)
        # linkyaxes!([ax_lefts[r], ax_rights[r]])

        hide_y_ticks!(ax_rights[r]);
        hide_ticks_and_spine!(ax_lefts[r])
        if r != how_many_rows
            hide_x_ticks!(ax_rights[r])
            hide_xlabels!(ax_lefts[r])
            hide_xlabels!(ax_rights[r])
        end
        change_x_tick_size!(ax_rights[r], text_size)
        change_x_tick_len!(ax_rights[r],  text_size)
    end
    hide_spine!.(ax_lefts); 
    linkxaxes!(ax_rights)
    colgap!(grid, 0);
end
