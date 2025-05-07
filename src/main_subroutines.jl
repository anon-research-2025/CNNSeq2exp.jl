
function objective(hp, m, seq, y; make_sparse=false)
    loss_here = huber_loss(
        model_forward_pass(hp, m, seq; make_sparse=make_sparse), 
        y; delta = 0.5)
    # println("loss: $loss_here")
    return loss_here
end

#=
hp: HyperParameters object
fp: file path to the hdf5 file
=#
function load_data(hp, fp)
    sd = sequence_data(fp);
    data_load = obtain_dataload_ith_mega_batch(sd, 1, "train"; batch_size=hp.batch_size);
    data_load_v = obtain_dataload_ith_mega_batch(sd, 1, "valid";  batch_size=hp.batch_size);
    data_load_t = obtain_dataload_ith_mega_batch(sd, 1, "test";  batch_size=hp.batch_size);
    return sd, data_load, data_load_v, data_load_t
end

function make_data_load(one_hot_arr, labels, batch_size; shuffle=true, partial=false)
    return Flux.DataLoader((one_hot_arr, labels), 
        batchsize=batch_size, shuffle=shuffle, partial=partial)
end

function make_data_load_w_v(data_load, data_load_v, hp; shuffle=true, partial=false)
    one_hot_arr=cat(data_load.data[1],
                data_load_v.data[1], dims=4)
    labels = cat(data_load.data[2],
                data_load_v.data[2],dims=2)
    return make_data_load(one_hot_arr, labels, hp.batch_size; shuffle=shuffle, partial=partial)
end

function make_data_load_all(data_load, data_load_v, data_load_t, hp; shuffle=false, partial=false)
    one_hot_arr=cat(data_load.data[1],
                data_load_v.data[1],
                data_load_t.data[1], dims=4)
    labels = cat(data_load.data[2],
                data_load_v.data[2],
                data_load_t.data[2], dims=2)
    return make_data_load(one_hot_arr, labels, hp.batch_size; shuffle=shuffle, partial=partial)
end

function train_model(hp, sd, data_load; n_epochs=48, make_sparse=false)
    m = model(hp, sd.seq_len);
    opt_state = Flux.setup(Flux.AdaBelief(), m);
    for e = 1:n_epochs
        @info "epoch: $e"
        for S in data_load
            seq, y = S[1:end-1][1], S[end]
            seq, y = seq |> gpu, y |> gpu
            gs = Flux.gradient(
                x -> objective(hp, x, seq, y; make_sparse=make_sparse), m)
            Flux.update!(opt_state, m, gs[1])
        end
    end
    return m
end

function train_model_further!(m, hp, data_load; n_epochs=5, epoch_offset=0, make_sparse=false)
    opt_state = Flux.setup(Flux.AdaBelief(), m);
    for e = 1:n_epochs
        @info "epoch: $(e+epoch_offset)"
        for S in data_load
            seq, y = S[1:end-1][1], S[end]
            seq, y = seq |> gpu, y |> gpu
            gs = Flux.gradient(
                x -> objective(hp, x, seq, y; make_sparse=make_sparse), m)
            Flux.update!(opt_state, m, gs[1])
        end
    end
end


function obtain_motifs(nz_dict_all, 
    nz_dict_positive, nz_dict_negative, df_B_seq,
    onehotarr, hp, m::model; 
    UPTO=UPTO,
    MAX_CONFIG_SIZE=MAX_CONFIG_SIZE)

    # declare motif instance
    _m_ = motifs(;upto=UPTO);
    # create nested dict data structure
    spf = create_spf(nz_dict_all)

    #=
    For singletons, just store the results in contributions in 
        positive contribution arrays
    =#
    obtain_singleton_motifs!(
        nz_dict_all, onehotarr, hp, _m_.pfms[1],
        df_B_seq, m::model; pos_info = _m_.pos[1], 
        power_info = _m_.pos_contribs[1]) 
    #= singleton motifs's contribution are all stored in positive contribution arrays 
    regardless of whether they are positive or negative (for now) =#

    for i = 2:UPTO
        pos_configs = get_configs(nz_dict_positive, hp; config_num_fil=i, max_config_size=MAX_CONFIG_SIZE);
        neg_configs = get_configs(nz_dict_negative, hp; config_num_fil=i, max_config_size=MAX_CONFIG_SIZE);
        # remove configs that are in both pos and neg
        rid_of_ambiguous_configs!(pos_configs, neg_configs)

        obtain_x_plet_motifs!(
            spf,
            onehotarr,
            pos_configs,
            _m_.pfms[i],
            df_B_seq,
            m,
            hp;
            pos_info = _m_.pos[i],
            power_info = _m_.pos_contribs[i]
        )
        
        obtain_x_plet_motifs!(
            spf,
            onehotarr,
            neg_configs,
            _m_.pfms[i],
            df_B_seq,
            m,
            hp;
            pos_info = _m_.pos[i],
            power_info = _m_.neg_contribs[i]
        )
    end
    return _m_
end

function results_render(_m_, onehotarr, labels, hp; save_path_where="yo", seq_save_name="seqs.fa", UPTO=UPTO)
    make_fasta_at_save_path(onehotarr, labels; save_path=save_path_where, save_name=seq_save_name)
    render_singletons(_m_, save_path_where)
    for t = 2:UPTO
        render_sorted_list(_m_, hp, save_path_where; t=t)
    end
end

function save_model(m, hp; save_where = "saved_models", save_name="cnn_yeast_gpro", trial_num=1)
    m_cpu = model2cpu(m)
    @save joinpath(save_where, "$(save_name)_$trial_num.jld2") m_cpu
    @save joinpath(save_where, "$(save_name)_hp_$trial_num.jld2") hp
end

function load_model(load_name; load_where="saved_models", trial_num=1)
    @load joinpath(load_where, "$(load_name)_$trial_num.jld2") m_cpu
    @load joinpath(load_where, "$(load_name)__hp_$trial_num.jld2") hp
    m = model2gpu(m_cpu)
    return m, hp
end

# function get_data_loads(fp, hp)
#     sd = sequence_data(fp)
#     data_load = obtain_dataload_ith_mega_batch(sd, 1, "train"; batch_size=hp.batch_size)
#     data_load_v = obtain_dataload_ith_mega_batch(sd, 1, "valid";  batch_size=hp.batch_size)
#     data_load_t = obtain_dataload_ith_mega_batch(sd, 1, "test";  batch_size=hp.batch_size)
#     dataload_w_v = make_data_load_w_v(data_load, data_load_v, hp)
#     dataload_all = make_data_load_all(data_load, data_load_v, data_load_t, hp)
#     return sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all
# end


function get_data_loads(fp, hp; w_v_shuffle=false)
    sd = sequence_data(fp)
    data_load = obtain_dataload_ith_mega_batch(sd, 1, "train"; batch_size=hp.batch_size)
    data_load_v = obtain_dataload_ith_mega_batch(sd, 1, "valid";  batch_size=hp.batch_size)
    data_load_t = obtain_dataload_ith_mega_batch(sd, 1, "test";  batch_size=hp.batch_size, 
        shuffle=false)
    dataload_w_v = make_data_load_w_v(data_load, data_load_v, hp; shuffle=w_v_shuffle)
    dataload_all = make_data_load_all(data_load, data_load_v, data_load_t, hp; shuffle=w_v_shuffle)
    return sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all
end


function obtain_B_flattened_both( hp, m, dataload_w_v, data_load_t; 
    do_sanity_chk = true, make_sparse=true)

    # make sure the data coming in is the shuffled off version
    # e.g. do 
    # _, _, _, _, _dataload_w_v_, _ = get_data_loads(fp, hp; w_v_shuffle=false)


    B_flattened, last_sequence_index = 
        obtain_B_flattened(hp, m, dataload_w_v; make_sparse=make_sparse)
    
    if do_sanity_chk
        sanity_chk(B_flattened, last_sequence_index, hp, m, dataload_w_v; 
            make_sparse=make_sparse)
    end

    B_flattened_t, _ = 
        obtain_B_flattened(hp, m, data_load_t; make_sparse=make_sparse);
        
    return B_flattened, last_sequence_index, B_flattened_t
end

"""
Inputs: 
    fp: file path
    hp: hyperparameters
        if nothing, then random hyperparameters are generated
    last_embedding_len_max: max length of the last embedding
    n_epochs: number of epochs to train the model
Outputs:
    hp: hyperparameters
    m: trained model
"""
function train_model(fp; 
        hp=nothing, 
        last_embedding_len_max=5, 
        make_sparse=true,
        n_epochs=48
        )
    if isnothing(hp)
        hp = random_hyperparameters()
        sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = 
            get_data_loads(fp, hp; w_v_shuffle=false)
        while get_last_embedding_len(data_load, hp) < 1 || get_last_embedding_len(data_load, hp) > last_embedding_len_max
            hp = random_hyperparameters()
        end
    else 
        sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = 
            get_data_loads(fp, hp; w_v_shuffle=false)
    end
    m = train_model(hp, sd, dataload_w_v; n_epochs=n_epochs, make_sparse=make_sparse);
    r2, r = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=make_sparse)
    @info "r2: $(r2), r: $(r)"
    return hp, m
end

"""
Render the motifs 
"""
function render_motifs(fp, hp, m, save_path; 
        code_mag_thresh=0.8, 
        contribution_thresh=0.75, 
        make_sparse=true,
        UPTO=4,
        seq_save_name="seqs.fa"
        )

    _, _, _, data_load_t, _dataload_w_v_, _dataload_all_ = 
        get_data_loads(fp, hp; w_v_shuffle=false)

    B_flattened, last_sequence_index, _ = 
        obtain_B_flattened_both(hp, m, _dataload_w_v_, data_load_t; 
            make_sparse=make_sparse);
    
    nz_dict_positive, nz_dict_negative, nz_dict_all, df_B, df_B_seq = 
        extract_nz_dicts_and_df_B_seq(B_flattened, last_sequence_index; 
            code_mag_quantile_thresh=code_mag_thresh, 
                contribution_quantile_thresh=contribution_thresh);

    @info "mean nz vals len: $(StatsBase.mean(length.(values(nz_dict_positive))))"

    onehotarr, labels = _dataload_all_.data

    _m_ = obtain_motifs(nz_dict_all, nz_dict_positive, 
        nz_dict_negative, df_B_seq, onehotarr, hp, m; UPTO=UPTO);
    
    update_avg_contrib_and_is_pos_contrib_all!(_m_)

    print_mean_contribution_singletons(_m_)

    results_render(_m_, onehotarr, labels, hp; 
        save_path_where=save_path, 
        seq_save_name=seq_save_name,
        UPTO=UPTO)
end


function get_linear_model(fp, hp, m; 
        w_v_shuffle=false, 
        normalized=true, 
        make_sparse=true
        )
    _, _, _, data_load_t, _dataload_w_v_, _ = 
        get_data_loads(fp, hp; w_v_shuffle=w_v_shuffle)

    B_flattened, _, B_flattened_t = 
        obtain_B_flattened_both(hp, m, _dataload_w_v_, data_load_t; make_sparse=make_sparse);
    
    flist = collect(1:hp.num_pfms);

    lms = my_full_and_reduced_lm(flist, B_flattened, _dataload_w_v_, hp; normalized=normalized)
    
    print_performance(lms, B_flattened_t, data_load_t, hp; plotyy=false);

    return lms
end


"""
fp: file path for the formatted hdf5 file;
    that contains the sequence and label information,
save_where: a folder path to save the results;
    note that the folder will be created if it does not exist.
UPTO: maximum number of filters in the coalition to be considered.
"""
function train_and_render(fp::String, 
    save_where::String;
    UPTO::Int=4 
    )
    hp, m = train_model(fp);
    render_motifs(fp, hp , m, save_where; UPTO=UPTO)
end

"""
fp: file path for the formatted hdf5 file;
    that contains the sequence and label information,
hp_load_path: path to the hyperparameters file;
    this file should be in the jld2 format
model_load_path: path to the model file;
    this file should be in the jld2 format
save_where: a folder path to save the results;
    note that the folder will be created if it does not exist.
UPTO: maximum number of filters in the coalition to be considered.
"""
function take_pretrained_and_render(
    fp::String,
    hp_load_path::String,
    model_load_path::String,
    save_where::String;
    UPTO::Int=4
)
    @load hp_load_path hp
    @load model_load_path m_cpu
    m = model2gpu(m_cpu)
    render_motifs(fp, hp , m, save_where; UPTO=UPTO)
end
