
function get_avg_ic(f; alpha = 1)
    mat_squared = exp.(f)
    mat_squared_sum = sum(mat_squared; dims = 1)
    normalize_mats = mat_squared ./ mat_squared_sum
    ic_cols = alpha .+ sum(normalize_mats .* log2.(normalize_mats); dims = 1)
    avg_ic = sum(ic_cols) / (size(mat_squared, 2) * size(mat_squared, 4))
    # return avg_ic
    return 0.5f0 * avg_ic
end

function construct_multihead_model(my_model::model, hp, L, num_heads)
    clen_last = get_last_lvl_pcode_len(hp, L) # (debug) for checking 
    @info "last embedding length $clen_last"
    computational_graph = Chain(
            Parallel(hcat, 
            [Chain(
                x->forward_pass_regression(hp, my_model, x),
                Dense(clen_last*hp.num_img_filters[end]=>num_heads)
                ) for _ in 1:num_heads]...),
            x->sum(x, dims=2),
            x->sigmoid.(x),
            x->reshape(x, (num_heads, size(x,3)))
        ) |> gpu
    return computational_graph
end

function huber_loss(ŷ, y; delta=float_type(0.85), p5=float_type(0.5), p5δ=p5 * delta)
    isnan_mask = Flux.Zygote.@ignore .!isnan.(y)
    num_valid_pts = Flux.Zygote.@ignore isnan_mask |> sum |> float_type
    abs_error = (abs.(ŷ - y))[isnan_mask]
    abs_error_mask = Flux.Zygote.@ignore abs_error .< delta
    n_abs_error_mask = Flux.Zygote.@ignore  .!abs_error_mask
    return sum(p5 .* (abs_error .* abs_error) .* abs_error_mask 
        .+ delta .* (abs_error .- p5δ) .* n_abs_error_mask ) / num_valid_pts
end