
data_load_noshuffle = dataload_all;


for (index, S) in enumerate(data_load_noshuffle)
    seq = S[1] |> gpu;
    code, gradCodeProduct, non_zero_coordinates = 
        obtain_gradCodeProduct_and_code(hp, m, seq)
    push_to_B_flattened!(B_flattened, code, gradCodeProduct, non_zero_coordinates, last_sequence_index)
    last_sequence_index += hp.batch_size
    @info "seq $(hp.batch_size*index) done"
end


S = first(data_load_noshuffle)
seq = S[1] |> gpu;
code, gradCodeProduct, non_zero_coordinates = 
    obtain_gradCodeProduct_and_code(hp, m, seq)

sum(code .== 0) / sum(prod(size(code)))
sum(code .!= 0) / sum(prod(size(code)))



for (index, S) in enumerate(data_load_noshuffle)
    seq = S[1] |> gpu;
    code, gradCodeProduct, non_zero_coordinates = 
        obtain_gradCodeProduct_and_code(hp, m, seq)

    zc_percent = sum(code .== 0) / sum(prod(size(code)))
    nz_percent = sum(code .!= 0) / sum(prod(size(code)))

    @info "non-zero percent: $(nz_percent) zero percent: $(round2(zc_percent))"
end




# check the sparsity in img layers as well
for (index, S) in enumerate(data_load_noshuffle)
    seq = S[1] |> gpu;
    code_imgs =
        obtain_code_imgs_this_batch(hp, m, seq; reverse_comp = false)
    str_here = sparsity_check(code_imgs)
    @info "$(str_here)"
end




S = first(data_load_noshuffle)
seq = S[1] |> gpu;
code_imgs =
    obtain_code_imgs_this_batch(hp, m, seq; reverse_comp = false)

function sparsity_check(code_imgs)
    zc = Vector{float_type}()
    nz = Vector{float_type}()
    for i in 1:length(code_imgs)
        zc_percent = sum(code_imgs[i] .== 0) / sum(prod(size(code_imgs[i])))
        nz_percent = sum(code_imgs[i] .!= 0) / sum(prod(size(code_imgs[i])))
        push!(zc, zc_percent)
        push!(nz, nz_percent)
        # @info "non-zero percent: $(nz_percent) zero percent: $(round2(zc_percent))"
    end
    str_here = ["$(round2(zc_i)) $(round2(n_iz)) -- " for (zc_i, n_iz) in zip(zc, nz)] |> prod
    return str_here
end

    code_imgs[1]
    code_imgs[2]
    code_imgs[3]
    code_imgs[4]
    code_imgs[5]