function prep_normalized_matrix(
    mat;
    sum_along_dim = 1,
    ϵ = float_type(1e-5),
    reverse_comp = false,
)
    # could be f or b
    mat_squared = mat .^ 2 .+ ϵ
    # mat_squared = exp.(mat) 
    mat_squared_sum = sum(mat_squared; dims = sum_along_dim)
    normalize_mats = mat_squared ./ mat_squared_sum
    if reverse_comp
        normalize_mats_rc = reverse(normalize_mats; dims = (1, 2))
        return cat(normalize_mats, normalize_mats_rc; dims = 4)
    end
    return normalize_mats
end

min_zero(x) = max.(0, x)
min_max_zero_bdd(x; bdd = 25) = min.(bdd, max.(0, x))
prep_pos_hyperparam(x) = min_max_zero_bdd(x .^ 2; bdd = float_type(0.5))

function prep_pwm(f, b; bg = float_type(0.25), reverse_comp = false, use_log = true)
    # return f
    if use_log
        return log2.(prep_normalized_matrix(f; reverse_comp = reverse_comp) ./ bg)
    else
        return prep_normalized_matrix(f; reverse_comp = reverse_comp) ./ bg
    end
end

function prep_pwm(f; bg = float_type(0.25), reverse_comp = false)
    return log.(prep_normalized_matrix(f; reverse_comp = reverse_comp) ./ bg)
end

function normalize_img_filters(f; alpha = 500, make_sparse = false)
    # abs_f = @ignore abs.(f)
    # f = f .* softmax(alpha .* abs_f; dims = 2)
    # f_norm = @ignore sqrt.(sum(f .^ 2; dims = (1, 2)))
    # return f ./ f_norm
    # if regression
    if make_sparse
        # println("normalizing img filters")
        abs_f = @ignore abs.(f)
        f = f .* softmax(alpha .* abs_f; dims = 2)
        f_norm = @ignore sqrt.(sum(f .^ 2; dims = (1, 2)))
        return f ./ f_norm
    else
        f_norm = @ignore sqrt.(sum(f .^ 2; dims = (1, 2)))
        return f ./ f_norm
    end
    
    # else
    #     expf = exp.(f)
    #     expf = expf .* softmax(alpha .* expf; dims = 2)
    #     expf_norm = sqrt.(sum(expf .^ 2; dims = (1, 2)))
    #     return expf ./ expf_norm
    # end
end

# return the stride tuple 
get_pool_i(hp, lvl) = lvl == 0 ? hp.pool_base : hp.poolsize[lvl]
get_pool(hp, lvl) = lvl == 0 ? (hp.pool_base, 1) : (hp.poolsize[lvl], 1)
get_stride_i(hp, lvl) = lvl == 0 ? hp.stride_base : hp.stride[lvl]
get_stride(hp, lvl) = lvl == 0 ? (hp.stride_base, 1) : (hp.stride[lvl], 1)
# ith level's pooled code length
# prime means the code is pooled
# im1: i minus 1
lvl_i_code_len(codelen_im1_prime, fi_len) = codelen_im1_prime - fi_len + 1
function lvl_i_pcode_len(code_len_i, hp, lvl)
    return (code_len_i - get_pool_i(hp, lvl)) ÷ get_stride_i(hp, lvl) + 1
end
base_code_len(hp, seq_len) = seq_len - hp.pfm_len + 1

function get_last_lvl_pcode_len(hp, seq_len)
    ci_code_len = base_code_len(hp, seq_len)
    ci_pcode_len = lvl_i_pcode_len(ci_code_len, hp, 0)
    for lvl = 1:get_num_lvl_abv_base(hp)
        if lvl ≤ hp.pool_lvl_top
            ci_code_len = lvl_i_code_len(ci_pcode_len, hp.img_fil_heights[lvl])
            ci_pcode_len = lvl_i_pcode_len(ci_code_len, hp, lvl)
        else
            ci_pcode_len = lvl_i_code_len(ci_pcode_len, hp.img_fil_heights[lvl])
        end
    end
    return ci_pcode_len
end

# function soft_max_sparse_code(
#     sparse_code;
#     alpha_base = 125.0f0,
#     alpha = 25.0f0,
#     between_layer = true,
#     last_layer = false,
#     base_layer = false,
# )
#     # TODO more principled way for alpha
#     if between_layer
#         return softmax(alpha .* sparse_code; dims = 3)
#     elseif base_layer
#         return softmax(alpha_base .* sparse_code; dims = (1, 3))
#     elseif last_layer
#         return softmax(alpha .* sparse_code; dims = (1, 3))
#     end
#     return error("no layer selected")
# end

pool_dim(input_dim, poolsize, pad, stride) = (input_dim + 2 * pad - poolsize) ÷ stride + 1

function code_img_pool(code_img; poolsize = (2, 1), stride = (1, 1))
    c, M, _, N = size(code_img)
    c_prime = @ignore pool_dim(c, poolsize[1], 0, stride[1]) # dimension after pooling
    pooled_code_img = reshape(
        Flux.NNlib.maxpool(code_img, poolsize; pad = 0, stride = stride),
        (c_prime, M, 1, N),
    )
    return pooled_code_img
end

function get_lvl_size(code; base_lvl = false)
    if base_lvl
        return size(code, 2), size(code, 3), size(code, 4)
    end
    return size(code, 1), size(code, 3), size(code, 4)
end

function code_pool(
    code;
    poolsize = (2, 1),
    stride = (1, 1),
    base_lvl = false,
    identity = false,
    w = nothing,
)
    c, M, N = get_lvl_size(code; base_lvl = base_lvl)
    identity && (return reshape(code, c, M, 1, N))
    code_img = reshape(code, c, M, 1, N)
    !isnothing(w) && (code_img = w .* code_img)
    pooled_code_img = code_img_pool(code_img; poolsize = poolsize, stride = stride)
    return pooled_code_img
end
