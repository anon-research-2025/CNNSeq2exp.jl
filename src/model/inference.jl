
function get_first_code_img(my_model::model, S::CuArray)
    code = my_model.pwms(S)
    return code
end

function img_pass(hp, my_model, code; make_sparse=false)
    _code_img_ = code_pool(
        code;
        poolsize = get_pool(hp, 0),
        stride = get_stride(hp, 0),
        base_lvl = true,
    )

    for lvl = 1:get_num_lvl_abv_base(hp)
        _code_ = my_model.img_filters[lvl](_code_img_, hp; make_sparse=make_sparse)
        _code_img_ = code_pool(
            _code_;
            poolsize = get_pool(hp, lvl),
            stride = get_stride(hp, lvl),
            identity = lvl > hp.pool_lvl_top,
        )
    end
    code_img_reshape = 
        reshape(_code_img_, (size(_code_img_, 1) * size(_code_img_, 2), 1, :))
    return code_img_reshape
end

function obtain_output_sum(hp, m::model, code_img; make_sparse=false)
    q = img_pass(hp, m, code_img; make_sparse=make_sparse) # code_img_reshape 
    output_sum =  m.beta[1] * sum(batched_mul(model_weights(m), q))
    return output_sum
end

function obtain_gradCodeProduct_and_code(hp, m, seq; make_sparse=false)
    code = get_first_code_img(m, seq)
    grad = gradient(x->obtain_output_sum(hp, m, x; make_sparse=make_sparse), code)[1]
    code = reshape(code, (size(code, 2), size(code, 3), size(code, 4)))
    grad = reshape(grad, (size(grad, 2), size(grad, 3), size(grad, 4)))
    gradCodeProduct = grad .* code;
    code, gradCodeProduct = code |> cpu, gradCodeProduct |> cpu
    non_zero_coordinates = findall(code .> 0) # pos, fil, seq
    return code, gradCodeProduct, non_zero_coordinates
end

function push_to_B_flattened!(B_flattened, code, gradCodeProduct, non_zero_coordinates, seq_offset)
    @inbounds for c in non_zero_coordinates
        pos, fil, seq = c[1], c[2], c[3] # pos, fil, seq
        contrib = gradCodeProduct[pos, fil, seq]
        code_mag = code[pos, fil, seq]
        tuple_here = (seq=seq+seq_offset, fil=fil, pos=pos, contribution=contrib, code_mag=code_mag)
        push!(B_flattened, tuple_here)
    end
end

function obtain_B_flattened(hp, m, data_load_noshuffle; make_sparse=false)
    # l = (size(dataload_w_v.data[2],2) ÷ hp.batch_size) * hp.batch_size
    # ys = Vector{float_type}(undef, (l,))

    B_flattened = Vector{B_flatten_t}()
    last_sequence_index = 0
    for (index, S) in enumerate(data_load_noshuffle)
        # ys[last_sequence_index+1:last_sequence_index+hp.batch_size] = S[2]'
        seq = S[1] |> gpu;
        code, gradCodeProduct, non_zero_coordinates = 
            obtain_gradCodeProduct_and_code(hp, m, seq; make_sparse=make_sparse)
        push_to_B_flattened!(B_flattened, code, gradCodeProduct, non_zero_coordinates, last_sequence_index)
        last_sequence_index += hp.batch_size
        @info "seq $(hp.batch_size*index) done"
    end
    return B_flattened, last_sequence_index
end


# B_flattened, last_sequence_index = obtain_B_flattened(hp, m, dataload_all)


function sanity_chk(B_flattened, last_sequence_index, hp, m, data_load_noshuffle;
        make_sparse=false)
    # store the retrieved predictions
    contribs_2_chk = zeros(float_type, (last_sequence_index,))
    for i in axes(B_flattened, 1)
        seq = B_flattened[i].seq;
        contribs_2_chk[seq] += B_flattened[i].contribution;
    end
    contribs_2_chk = tanh.(contribs_2_chk)

    # store the predictions from the model
    seq_offset = 0
    model_outputs = zeros(float_type, (last_sequence_index,))
    for S in data_load_noshuffle
        predicted_outputs = 
            model_forward_pass(hp, m::model, S[1] |> gpu; 
                make_sparse=make_sparse) |> cpu
        model_outputs[seq_offset+1:seq_offset+hp.batch_size] = predicted_outputs
        seq_offset += hp.batch_size
    end
    differences = contribs_2_chk .- model_outputs
    maximum_difference = maximum(abs.(differences))
    approximately_correct_percentage = 
        sum(contribs_2_chk .≈ model_outputs) / length(model_outputs) * 100
    @info "approximately correct percentage: $(round2(approximately_correct_percentage))"
    @info "max difference: $maximum_difference"
end

# sanity_chk(B_flattened, last_sequence_index, hp, m, dataload_all)


# [tanh.(contribs_2_chk)[1:128] cpu(o')]



# o = model_forward_pass(hp, m::model, S[1] |> gpu)
# S[2]

# [o' S[2]']



# get_all_coordinates(code)


# code_cpu = code |> cpu
# grad_code_product_cpu = grad_code_product |> cpu


