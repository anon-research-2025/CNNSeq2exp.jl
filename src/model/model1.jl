struct learned_pwms
    f::AbstractArray{float_type,4} # pfms
    b::AbstractArray{float_type,4} # background
    w::AbstractArray{float_type,4} # pos-neg weights for regression
    etas::AbstractArray{float_type, 1}
    function learned_pwms(;
        q = 6,
        M = 48,
        ϵ = 1e-1,
        cuda = false,
    )
        # pool_len = 0 means no pooling
        f = ϵ .* randn(float_type, (4, q, 1, M))
        b = fill(float_type(0.25), (4, q, 1, 1))
        w = randn(float_type, (1, M, 1, 1))
        etas = [ϵ * rand(float_type)]
        if cuda
            f = cu(f)
            b = cu(b)
            w = cu(w)
        end
        return new(f, b, w, etas)
    end
    function learned_pwms(f, b, w, etas)
        return new(f, b, w, etas)
    end
end


Flux.@layer learned_pwms

function first_lvl_prep_param(lpwms::learned_pwms; reverse_comp = false, use_log = true)
    return prep_pwm(lpwms.f, lpwms.b; reverse_comp = reverse_comp, use_log = use_log),
            prep_pos_hyperparam(lpwms.etas)
end

#=
Init the sparse code according to the learned PWM objective
=#
init_code_learned_pwm(grad, etas) = Flux.NNlib.relu(etas[1] .* grad)

function first_lvl_forward_pass(
    lpwms,
    S;
    return_pwms = false,
    reverse_comp = false,
    use_log = true,
)
    pwms, etas =
        first_lvl_prep_param(lpwms; reverse_comp = reverse_comp, use_log = use_log)
    grad = convolve(S, pwms; pad = 0, flipped = true)
    code = init_code_learned_pwm(grad, etas)
    return_pwms && (return code, pwms)
    return code
end

function (lpwms::learned_pwms)(S; return_pwms = false, use_log = true)
    return first_lvl_forward_pass(lpwms, S; return_pwms = return_pwms, use_log = use_log)
end

####### IMAGE LAYER #######

struct learned_codeimg_filters
    f::AbstractArray{float_type,4}
    eta::float_type
    function learned_codeimg_filters(
        code_img_1_width;
        h = 6,
        num_filters = 24,
        ϵ = 1e-3,
        cuda = false,
    )
        filters = ϵ .* randn(float_type, (h, code_img_1_width, 1, num_filters))
        eta = ϵ .* rand(float_type)
        cuda && (filters = cu(filters))
        return new(filters, eta)
    end
    function learned_codeimg_filters(f, eta)
        return new(f, eta)
    end
end

Flux.@layer learned_codeimg_filters

function i_lvl_prep_param(img_fils_obj; alpha=float_type(500), make_sparse=false)
    img_filters = normalize_img_filters(img_fils_obj.f; alpha = alpha, make_sparse = make_sparse)
    eta = prep_pos_hyperparam(img_fils_obj.eta)
    return img_filters, eta
end

init_code_img_fils(grad) = Flux.NNlib.relu(grad)

function i_lvl_forward_pass(
    img_fils_obj,
    prev_code_img, hp;
    return_filters = false,
    make_sparse=false,
)
    img_filters, _ = i_lvl_prep_param(img_fils_obj; 
        alpha = hp.softmax_strength_img_fil, make_sparse=make_sparse)
    grad = convolve(prev_code_img, img_filters; pad = 0, flipped = true)
    code = init_code_img_fils(grad)    
    return_filters && (return code, img_filters)
    return code
end

function (img_fils::learned_codeimg_filters)(
    prev_code_img, hp;
    return_filters = false,
    make_sparse=false,
)
    return i_lvl_forward_pass(
        img_fils,
        prev_code_img, hp;
        return_filters = return_filters,
        make_sparse=make_sparse,
    )
end

####### MODEL #######

function forward_pass_regression(hp, my_model, S::CuArray; make_sparse=false)
    code = my_model.pwms(S)
    code_img = code_pool(
        code;
        poolsize = get_pool(hp, 0),
        stride = get_stride(hp, 0),
        base_lvl = true,
    )

    for lvl = 1:get_num_lvl_abv_base(hp)
        code = my_model.img_filters[lvl](code_img, hp; make_sparse=make_sparse)
        code_img = code_pool(
            code;
            poolsize = get_pool(hp, lvl),
            stride = get_stride(hp, lvl),
            identity = lvl > hp.pool_lvl_top,
        )
    end
    code_img_reshape = reshape(code_img, (size(code_img, 1) * size(code_img, 2), 1, :))

    return code_img_reshape
end



struct model
    pwms::learned_pwms
    img_filters::Vector{learned_codeimg_filters}
    w::AbstractArray{float_type, 3}
    beta::AbstractArray{float_type, 1}
    function model(hp, seq_len::Integer; ϵ=float_type(5e-1), cuda=true)
        pwms = learned_pwms(;
            q = hp.pfm_len,
            M = hp.num_pfms,
            ϵ = ϵ,
            cuda = cuda,
        )

        img_filters = [
            learned_codeimg_filters(
                w;
                h = h,
                num_filters = nf,
                ϵ = ϵ,
                cuda = cuda,
            ) for (w, h, nf) in zip(
                hp.img_fil_widths,
                hp.img_fil_heights,
                hp.num_img_filters,
            )
        ]

        beta = [float_type(1)]

        clen_last = get_last_lvl_pcode_len(hp, seq_len) # (debug) for checking 
        w = randn(float_type, (1, clen_last*hp.num_img_filters[end], 1)) |> cu

        return new(pwms, img_filters, w, beta)
    end
    function model(pwms, img_filters, w, beta)
        return new(pwms, img_filters, w, beta)
    end
end

Flux.@layer model

model_weights(m::model) = m.w

function model_forward_pass(hp, m::model, seq; make_sparse=false)
    # q is the vectorized last conv-layer's embedding, i.e. size is (l x m, 1, n)
    q = forward_pass_regression(hp, m, seq; make_sparse=make_sparse) # code_img_reshape 
    netowrk_output = output_function.(m.beta[1] .* sum(batched_mul(model_weights(m), q), dims=2))
    output = reshape(netowrk_output, (1, size(netowrk_output, 3)))
    return output
end

(m::model)(hp, seq; make_sparse=false) = 
    model_forward_pass(hp, m, seq; make_sparse=make_sparse)



