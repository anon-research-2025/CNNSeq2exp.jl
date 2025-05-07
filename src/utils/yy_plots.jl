
# using StatsBase
# using Plots

function rsq(ŷ, y)
    ȳ = mean(y)
    SS_tot = sum((y .- ȳ) .^ 2)
    SS_res = sum((y .- ŷ) .^ 2)
    return 1 - SS_res / SS_tot
end

function pcc(ŷ, y)
    y_mean = mean(y)
    yhat_mean = mean(ŷ)
    num = sum((y .- y_mean) .* (ŷ .- yhat_mean))
    den = sqrt(sum((y .- y_mean) .^ 2) * sum((ŷ .- yhat_mean) .^ 2))
    return num / den
end

function get_predicted_and_y_true(data_load, m, hp; make_sparse=false)
    predicted_y = float_type[]
    y_true = float_type[]
    for S in data_load
        seq, y = S[1:end-1][1], S[end]
        seq, y = seq |> gpu, y |> gpu
        ŷ = model_forward_pass(hp, m, seq; make_sparse=make_sparse)
        append!(predicted_y, cpu(ŷ))
        append!(y_true, cpu(y))
    end
    mask = .!isnan.(y_true)
    predicted_y = predicted_y[mask]
    y_true = y_true[mask]
    return predicted_y, y_true
end

function plot_yy(data_load, m, hp; plotyy = false, make_sparse=false)
    predicted_y, y_true = get_predicted_and_y_true(data_load, m, hp, make_sparse=make_sparse)
    r2 = rsq(predicted_y, y_true)
    r = pcc(predicted_y, y_true)
    @info "r2: $r2, r $r"
    if !plotyy
        return r2, r
    else
        Plots.scatter(predicted_y, y_true)
        Plots.annotate!(maximum(y_true)-0.2, 0.1, 
            Plots.text("r = $(round(r, digits=3))\nR² = $(round(r2, digits=3))", 10))
    # savefig("re_now/yy.png")
    end
end
