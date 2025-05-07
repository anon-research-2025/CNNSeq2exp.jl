Base.@kwdef struct HyperParameters
    pfm_len::Int = 10
    num_pfms::Int = 24
    num_img_filters::Vector{Int} = [65, 98, 128, 128, 76, 5]
    img_fil_widths::Vector{Int} = vcat([num_pfms], num_img_filters[1:(end-1)])
    img_fil_heights::Vector{Int} = [6, 6, 6, 6, 6, 5] # yeast
    pool_base::Int = 2
    stride_base::Int = 1
    poolsize::Vector{Int} = [2, 2, 2, 2, 2, 1] # yeast
    stride::Vector{Int} =   [1, 1, 2, 2, 2, 1] # yeast
    pool_lvl_top::Int = 5
    softmax_strength_img_fil::float_type = 500.0
    batch_size::Int = 256
end

# get the number of layers above the base layer
get_num_lvl_abv_base(hp) = length(hp.num_img_filters)

function get_input_sequence_length(data)
    input_seq_len = size(data.data_matrices_full[1], 2)
    return input_seq_len
end

function get_last_embedding_len(data_load, hp)
    embedding_len = 
    len_fil_conv_then_pool_stride(
        size(data_load.data[1],2), hp.pfm_len, 
        hp.pool_base, 
        hp.stride_base)

    for i = 1:get_num_lvl_abv_base(hp)
        embedding_len =
            len_fil_conv_then_pool_stride(
                embedding_len, hp.img_fil_heights[i], 
                hp.poolsize[i], hp.stride[i])
    end
    return embedding_len
end

# tuning 
function random_hyperparameters()
    pfm_len = rand(5:10)  # Random length between 5 and 10
    num_pfms = rand(8:16)  # Random number of PFMs between 8 and 16
    num_layers = 6
    num_img_filters = [rand(64:256) for _ in 1:(num_layers-1)]  # Random filter sizes
    push!(num_img_filters, 5)  # Keep the last filter size fixed

    img_fil_widths = vcat([num_pfms], num_img_filters[1:(end-1)])
    img_fil_heights = [rand(5:10) for _ in 1:num_layers]  # Heights in range 4 to 8
    # img_fil_heights[end] = 5  # Keep the last height fixed

    pool_base = 2
    stride_base = 1
    poolsize = [rand(1:3) for _ in 1:(num_layers-1)]  # Random pool sizes
    push!(poolsize, 1)  # Keep last pool size fixed

    stride = [rand(1:3) for _ in 1:(num_layers-1)]  # Random strides
    push!(stride, 1)  # Keep last stride fixed

    pool_lvl_top = 5 
    batch_size = 128  # Choose batch size from 64, 128, 256

    return HyperParameters(
        pfm_len=pfm_len,
        num_pfms=num_pfms,
        num_img_filters=num_img_filters,
        img_fil_widths=img_fil_widths,
        img_fil_heights=img_fil_heights,
        pool_base=pool_base,
        stride_base=stride_base,
        poolsize=poolsize,
        stride=stride,
        pool_lvl_top=pool_lvl_top,
        softmax_strength_img_fil=500.0,  # Fixed
        batch_size=batch_size
    )
end


# pooling and stride strategy for CNN hyperparameters
struct pooling_strategy
    pool_base::Int        # default pool for base layer
    stride_base::Int      # default stride for base layer
    pool::Int             # default pool for img layers 
    stride::Int           # default stride for img layers  
    num_img_stride_1::Int #= number of stride 1 layers right above the 
                                   base layer =#
    function pooling_strategy()
        new(2, 2, 2, 2, 0)
    end
    function pooling_strategy(pool_base, stride_base, pool, stride, num_img_stride_1)
        new(pool_base, stride_base, pool, stride, num_img_stride_1)
    end
end

# pooling strategies
TwoByOneFirst_Pooling() = pooling_strategy(2, 1, 2, 2, 0)
TwoByOneX_Pooling(x::Int) = pooling_strategy(2, 1, 2, 2, x)
TwoByTwoAllTheWay_Pooling() = pooling_strategy(2, 2, 2, 2, 0)
TwoByOneFirst_KByKAllTheWay_Pooling(k::Int) = pooling_strategy(2, 1, k, k, 0)
KByKAllTheWay_Pooling(k::Int) = pooling_strategy(k, k, k, k, 0)

# get the code length after filter convolution
len_fil_conv(input_len, fil_len) = input_len - fil_len + 1
# get the code length after pooling and stride
len_pool_stride(input_len, pool, stride) = ((input_len - pool) ÷ stride) + 1 # assume no zero padding
# get the code length after filter convolution, pooling, and stride
len_fil_conv_then_pool_stride(input_len, fil_len, pool, stride) =
    len_pool_stride(len_fil_conv(input_len, fil_len), pool, stride)

#= 
    calculate the last embedding length 
    and the number of image layers for constructing 
    the hyperparameters
=#
function get_last_embedding_len(
    input_seq_len::Int, 
    _pooling_strategy_::pooling_strategy,;
    pfm_len = 7,
    img_fil_height = 5,
)
    embedding_len = 
        len_fil_conv_then_pool_stride(
            input_seq_len, pfm_len, 
            _pooling_strategy_.pool_base, 
            _pooling_strategy_.stride_base)
    num_img_layers = 0
    @info "embedding_len: $embedding_len"

    while true
        embedding_len_next =
            len_fil_conv_then_pool_stride(
                embedding_len, img_fil_height, 
                _pooling_strategy_.pool, _pooling_strategy_.stride)
        if embedding_len_next > img_fil_height
            embedding_len = embedding_len_next
            num_img_layers += 1
            @info "embedding_len: $embedding_len"
        else
            break
        end
    end

    # while last_embedding_len > img_fil_height # TODO take pool & stride into account
    #     # stride_here = num_img_layers < _pooling_strategy_.num_img_stride_1 ? 
    #     #     1 : _pooling_strategy_.pool

    #     last_embedding_len =
    #         len_fil_conv_then_pool_stride(
    #             last_embedding_len, img_fil_height, 
    #             _pooling_strategy_.pool, _pooling_strategy_.stride)
    #     @info "last_embedding_len: $last_embedding_len"
    #     num_img_layers += 1
    # end

    return embedding_len, num_img_layers
end

#=
Construct the pool sizes and stride sizes 
for the image layers for hyperparameters construction
=#
function construct_img_layers_pool_stride_hyperparam(
    input_seq_len::Int,
    _pooling_strategy_::pooling_strategy;
    pfm_len = 7,
    img_fil_height = 5
)
    @info "pfm_len: $pfm_len, img_fil_height: $img_fil_height"
    last_embedding_len, num_img_layers = 
        get_last_embedding_len(
            input_seq_len, _pooling_strategy_;
            pfm_len = pfm_len, img_fil_height = img_fil_height)
    @info "embedding before the last has length: $last_embedding_len"

    pool_sizes = fill(_pooling_strategy_.pool, num_img_layers)
    stride_sizes = fill(_pooling_strategy_.stride, num_img_layers)
        # vcat(
        #     ones(Int, _pooling_strategy_.num_img_stride_1),
        #     fill(_pooling_strategy_.stride, 
        #         num_img_layers - _pooling_strategy_.num_img_stride_1)
        # )
    @info "making the last embedding to have 1 row"
    pool_sizes = vcat(pool_sizes, 1)
    stride_sizes = vcat(stride_sizes, 1)
    num_img_layers += 1
    return pool_sizes, stride_sizes, num_img_layers, last_embedding_len
end

function get_pool_lvl_top(pool_sizes, stride_sizes)
    len_poolsize = pool_sizes |> length
    pool_lvl_top = pool_sizes |> length
    for i in len_poolsize:-1:1
        if pool_sizes[i] == 1 && stride_sizes[i] == 1
            pool_lvl_top -=1
        else
            break
        end
    end
    return pool_lvl_top
end

function construct_hyperparam(
    input_seq_len::Int, 
    _pooling_strategy_::pooling_strategy;
    pfm_len = 7,
    img_fil_height = 5,
    num_pfms = 48,
    num_img_fils_hidden = 72,
    num_img_fils_last = 5,
    softmax_strength_img_fil = float_type(500.0)
)
    # input_seq_len = get_input_sequence_length(data)
    @info "input sequence length: $input_seq_len"
    # construct the pool sizes and stride sizes for the image layers
    pool_sizes, stride_sizes, num_img_layers, last_embedding_len = 
        construct_img_layers_pool_stride_hyperparam(
            input_seq_len, _pooling_strategy_;
            pfm_len = pfm_len, img_fil_height = img_fil_height)
    # the last layer uses a different number of filters (num_img_fils_last)
    num_img_fils = vcat(fill(num_img_fils_hidden, num_img_layers - 1), num_img_fils_last)
    # the width of the image filters is the same as the number of filters in the previous layer
    img_fil_widths  = vcat([num_pfms], num_img_fils[1:(end-1)])
    img_fil_heights = vcat(fill(img_fil_height, num_img_layers-1), last_embedding_len)
    HyperParameters(pfm_len, 
                    num_pfms, 
                    num_img_fils, 
                    img_fil_widths, 
                    img_fil_heights,
                    _pooling_strategy_.pool_base,
                    _pooling_strategy_.stride_base,
                    pool_sizes, stride_sizes,
                    get_pool_lvl_top(pool_sizes, stride_sizes),
                    softmax_strength_img_fil
                    )
end
