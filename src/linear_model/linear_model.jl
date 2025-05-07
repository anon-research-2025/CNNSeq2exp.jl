
#=
given, e.g. 
    feature_list = [6, 3, 5, (2,6), (4,6), (6,9), (2,2,6), (2,6,6)]
    where each element is either an int or a tuple
    indicating filter or product of filters
=#
const f_type = Vector{Union{Int, NTuple{2, Int}, NTuple{3, Int}}};


# get the largest sequence index
function get_largest_seq_index(B_flattened)
    largest_seq_index = 0
    for b in B_flattened
        (b.seq > largest_seq_index) && (largest_seq_index = b.seq)
    end
    return largest_seq_index
end

#=
    Partition B_flattened in to a list 
    where the list[seq_index] correspond to 
    all the activations in that sequence
=#
function partition_B_flattened(B_flattened)
    largest_seq_index = get_largest_seq_index(B_flattened)
    partitioned = [Vector{eltype(B_flattened)}() for _ = 1:largest_seq_index]
    for b in B_flattened
        push!(partitioned[b.seq], b)
    end
    return partitioned
end

function get_mags(partition; num_filters=9)
    mags = zeros(float_type, (num_filters,))
    for b in partition
        mags[b.fil] += b.code_mag
    end
    return mags
end

#=
    map the magnitudes of filter activations
=#
map_magnitudes(mags, x) = x isa Int ? mags[x] : map(y->mags[y], x) |> prod

function prep_design_matrix(feature_list::f_type, 
    partitioned;
    num_filters=9, normalized=true, linear_model=nothing)
    
    num_features = length(feature_list) # number of features
    X = zeros(float_type, (length(partitioned), num_features)) # design matrix
    
    for (seq_index, partition) in enumerate(partitioned)
        mags = get_mags(partition; num_filters=num_filters)
        for (feature_index, feature) in enumerate(feature_list)
            X[seq_index, feature_index] = map_magnitudes(mags, feature)
        end    
    end
    # normalize the design matrix
    X_mean, X_std = nothing, nothing
    if normalized 
        if isnothing(linear_model)
            X_mean = StatsBase.mean(X, dims=1)
            X_std = StatsBase.std(X, dims=1)
            X = (X .- X_mean) ./ X_std
        else
            X = (X .- linear_model.X_mean) ./ linear_model.X_std
        end
    end

    # prepend column of ones for the intercept
    X = hcat(ones(float_type, size(X, 1)), X)  # prepend column of ones
    return X, X_mean, X_std
end

struct my_lmodel
    feature_list::f_type
    X::Matrix{float_type} # design matrix
    X_mean::Union{Matrix{float_type}, Nothing} # a 1 x num_filters matrix
    X_std::Union{Matrix{float_type}, Nothing} # a 1 x num_filters matrix
    sse::float_type
    df::Int
    beta::Vector{float_type}
    function my_lmodel(feature_list, B_flattened, data_load, hp; normalized=true)
        # feature_list e.g. [6, 3, 5, (2,6), (4,6), (6,9), (2,2,6), (2,6,6)]
        feature_list = f_type(feature_list)
        partitioned = partition_B_flattened(B_flattened);
        X, X_mean, X_std = prep_design_matrix(feature_list, partitioned; 
            num_filters = hp.num_pfms, normalized = normalized, linear_model=nothing)
        @info "X_mean: $(round.(X_mean, digits=2))"
        @info "X_std: $(round.(X_std, digits=2))"
        y = (data_load.data[2] |> vec)[1:size(X,1)]
        beta = X \ y

        sse = sum((y .- X * beta) .^ 2)
        df = size(X, 1) - size(X, 2)

        new(feature_list, X, X_mean, X_std, sse, df, beta)
    end
end


function t_test_beta(lm::my_lmodel)
    sigma_squared = lm.sse / lm.df
    standard_errors = sqrt.(sigma_squared * (1 ./ diag(lm.X' * lm.X)))
    t_stat = lm.beta ./ standard_errors
    tdist = TDist(lm.df)
    p_values = 2 .* (1 .- cdf(tdist, abs.(t_stat)))
    coef_names = ["Intercept"; ["x$j" for j in 1:size(lm.X, 2)-1]...]  # adjust if needed
    summary_table = DataFrame(
        Term = coef_names,
        Estimate = lm.beta,
        StdError = standard_errors,
        tStat = t_stat,
        pValue = p_values
    )
    return summary_table    
end


struct my_full_and_reduced_lm
    full_lm::my_lmodel
    reduced_lm::my_lmodel
    function my_full_and_reduced_lm(feature_list, B_flattened, data_load, hp; normalized=true)
        reduced_model_feature_indicator = length.(feature_list) .== 1
        full_lm = my_lmodel(feature_list, B_flattened, data_load, hp; normalized=normalized)
        reduced_lm = my_lmodel(feature_list[reduced_model_feature_indicator], 
            B_flattened, data_load, hp; normalized=normalized)
        new(full_lm, reduced_lm)
    end
end

function return_yy(lm, B_flattened_t, dataload_t, hp; normalized=true)
    partitioned = partition_B_flattened(B_flattened_t);
    X,_,_ = prep_design_matrix(lm.feature_list, partitioned; 
        num_filters = hp.num_pfms, normalized=normalized, linear_model=lm)  
    y = (dataload_t.data[2] |> vec)[1:size(X,1)]
    predicted = X * lm.beta
    return predicted, y
end


function pretty_print_tuples(lm)
    @info "b0: $(round(lm.beta[1], digits=2))"
    for (index, (a, b)) in enumerate(zip(lm.feature_list, lm.beta[2:end]))
        a_str = isa(a, Tuple) ? string(a) : string(a)
        b_str = @sprintf("%.2f", b)
        @info "$a_str, b$index: $b_str"
    end
end

#=
 test whether the full model is significantly different from the reduced model
    if p small then the full model works better (check r^2)
    else the reduced model works better (check r^2)
 =#
 function F_test(lms)
    v1 =  (lms.reduced_lm.df - lms.full_lm.df)
    v2 = lms.full_lm.df
    t1 = (lms.reduced_lm.sse - lms.full_lm.sse) / v1
    t2 = lms.full_lm.sse / v2
    f_star = t1 / t2 
    p_value = 1 - cdf(FDist(v1, v2), f_star)
    @info "F-test p_value: $p_value"
end

function plot_scatter_yy(predicted, y, r, r2)
    Plots.scatter(predicted, y; 
    title="Predicted vs Actual",
    xlabel="Predicted",
    ylabel="Actual",)
    Plots.annotate!(maximum(y)-0.2, 0.1, 
    Plots.text("r = $(round(r, digits=3))\nR² = $(round(r2, digits=3))", 10))
end

function print_performance(lms, B_flattened_t, data_load_t, hp; plotyy=false)
    predicted, y = return_yy(lms.full_lm, B_flattened_t, data_load_t, hp;
        normalized=!isnothing(lms.full_lm.X_std))

    r2 = round(rsq(predicted, y), digits=2)
    r =  round(pcc(predicted, y), digits=2)
    @info "full model: r2: $r2, r $r"
    plotyy && plot_scatter_yy(predicted, y, r, r2)

    predicted, y = return_yy(lms.reduced_lm, B_flattened_t, data_load_t, hp;
        normalized=!isnothing(lms.reduced_lm.X_std))
    r2 = round(rsq(predicted, y), digits=2)
    r =  round(pcc(predicted, y), digits=2)
    @info "reduced model: r2: $r2, r $r"
    plotyy && plot_scatter_yy(predicted, y, r, r2)
end


function get_full_model_predicted_and_label(lms, B_flattened_t, data_load_t, hp)
    predicted, y = return_yy(lms.full_lm, B_flattened_t, data_load_t, hp;
        normalized=!isnothing(lms.full_lm.X_std))
    r2 = round(rsq(predicted, y), digits=3)
    r =  round(pcc(predicted, y), digits=3)
    return predicted, y, r, r2
end