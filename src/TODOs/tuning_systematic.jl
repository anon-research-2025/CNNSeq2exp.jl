
using HDF5
using Random
using Flux
using Zygote: @ignore
using CUDA
using StatsBase
using ShaneGPUCountMinSketch

using DataFrames, PlotPWM, CSV
using Printf, JSON3, Mustache
using CairoMakie, ImageShow
using JLD2

const load_how_many = 50000
const convolve = Flux.NNlib.conv

include("const.jl")
include("model/hp1.jl")
include("model/helper1.jl")
include("model/model1.jl")
include("model/forward.jl")
include("model/inference.jl")
include("load.jl")
include("yy_plots.jl")

include("contribution.jl")
include("traverse.jl")

include("motif_analysis/motif.jl")
include("motif_analysis/helper.jl")
include("render/const.jl")
include("render/onehot2fasta.jl")
include("render/logo_save.jl")
include("render/fill_json.jl")
include("render/template_css.jl")
include("render/template_script.jl")
include("render/templates_html.jl")
include("render/boxplots.jl")
include("render/render.jl")

include("save_model.jl")

include("main.jl")

regression_data_folder = "/home/shane/Desktop/seq2exp/seq2exp_datasets/processed_data_tanh"

function get_data_loads(fp, hp)
    sd = sequence_data(fp)
    data_load = obtain_dataload_ith_mega_batch(sd, 1, "train"; batch_size=hp.batch_size)
    data_load_v = obtain_dataload_ith_mega_batch(sd, 1, "valid";  batch_size=hp.batch_size)
    data_load_t = obtain_dataload_ith_mega_batch(sd, 1, "test";  batch_size=hp.batch_size)
    dataload_w_v = make_data_load_w_v(data_load, data_load_v, hp)
    dataload_all = make_data_load_all(data_load, data_load_v, data_load_t, hp)
    return sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all
end


fp_paths = [joinpath(regression_data_folder, "Aviv/yeast/yeast.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cgan_wanglab/ecoli_cgan.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cross_species_wanglab/ecoli_cross.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/ELF1/ELF1.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF4A/HNF4A.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF1A/HNF1A.hdf5")]


fp_names = ["cnn_yeast", 
            "ecoli_cgan", 
            "ecoli_cross", 
            "ELF1", 
            "HNF4A", 
            "HNF1A"]

r2_bests = float_type.([0, 0, 0, 0, 0, 0])
r2_bests_sparse = float_type.([0, 0, 0, 0, 0, 0])
r_bests =  float_type.([0, 0, 0, 0, 0, 0])
r_bests_sparse =  float_type.([0, 0, 0, 0, 0, 0])
result_tuple = []

for _ = 1:10000

    # randomly choose a file path
    rand_fp_index = rand(1:length(fp_paths))
    f_name = fp_names[rand_fp_index]
    @info "f_name: $(f_name)"
    fp = joinpath(regression_data_folder, fp_paths[rand_fp_index]);

    hp = random_hyperparameters();
    sd, data_load, data_load_v, data_load_t = load_data(hp, fp);
    while get_last_embedding_len(data_load, hp) < 1 ||  get_last_embedding_len(data_load, hp) > 5
        hp = random_hyperparameters();
        # sd, data_load, data_load_v, data_load_t = load_data(hp, fp) # batch size is fixed; no concerns
    end
    @info "last embedding len $(get_last_embedding_len(data_load, hp))"

    m = train_model(hp, sd, data_load; n_epochs=48, make_sparse=false);
    r2_ns, r_ns = plot_yy(data_load_v, m, hp; plotyy=false, make_sparse=false)

    # for i = 1:16
    #     train_model_further!(m, hp, data_load; n_epochs=1, epoch_offset=i)
    #     r2_ns, r_ns = plot_yy(data_load_v, m, hp; plotyy=false, make_sparse=false)
    # end

    if r2_ns > r2_bests[rand_fp_index] && r_ns > r_bests[rand_fp_index]
        m_cpu = model2cpu(m)
        @save "saved_models_2/$(f_name)_best.jld2" m_cpu
        @save "saved_models_2/$(f_name)_best_hp.jld2" hp
        r2_bests[rand_fp_index] = r2_ns
        r_bests[rand_fp_index] = r_ns 
    end

    m = train_model(hp, sd, data_load; n_epochs=48, make_sparse=true);
    r2_s, r_s = plot_yy(data_load_v, m, hp; plotyy=false, make_sparse=true)

    # for i = 1:16
    #     train_model_further!(m, hp, data_load; n_epochs=1, epoch_offset=i, make_sparse=true)
    #     r2_s, r_s = plot_yy(data_load_v, m, hp; plotyy=false, make_sparse=true)
    # end

    if r2_s > r2_bests_sparse[rand_fp_index] && r_s > r_bests_sparse[rand_fp_index]
        m_cpu = model2cpu(m)
        @save "saved_models_2/$(f_name)_best_sparse.jld2" m_cpu
        @save "saved_models_2/$(f_name)_best_hp_sparse.jld2" hp
        r2_bests_sparse[rand_fp_index] = r2_s
        r_bests_sparse[rand_fp_index] = r_s 
    end

    result_here = (hp=hp, f_name=f_name, r2_ns=r2_ns, r_ns=r_ns, r2_s=r2_s, r_s=r_s)
    push!(result_tuple, result_here)

end

@save "saved_models_2/result_tuples.jld2" result_tuple
@save "saved_models_2/names.jld2" fp_names
@save "saved_models_2/r2_bests.jld2" r2_bests
@save "saved_models_2/r2_bests_sparse.jld2" r2_bests_sparse
@save "saved_models_2/r_bests.jld2" r_bests
@save "saved_models_2/r_bests_sparse.jld2" r_bests_sparse


@load "saved_models_2/result_tuples.jld2" result_tuple

# b1, b2 = r2_bests[1], r_bests[1]

# df[df.f_name .== "cnn_yeast",:].r2_s |> maximum
# df[df.f_name .== "cnn_yeast",:].r2_ns |> maximum
# df_here = df[df.f_name .== "cnn_yeast",:]
# max_ind = df_here.r2_s |>  argmax
# hp = df_here[max_ind,:].hp
# @save "saved_models_2/cnn_yeast_best_hp_sparse.jld2" hp


# @load "saved_models_2/cnn_yeast_best_hp.jld2" hp

# @load "saved_models_2/cnn_yeast_best_hp_sparse.jld2" hp


# fp = fp_paths[1];

# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=true);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=true)
# m_cpu = model2cpu(m);

# @save "saved_models_2/further_trained/yeast_best_sparse.jld2" m_cpu



# fp = joinpath(regression_data_folder, fp_paths[2]);

# @load "saved_models_2/ecoli_cgan_best_hp.jld2" hp

# @load "saved_models_2/ecoli_cgan_best_hp_sparse.jld2" hp

# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=false);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=false)
# m_cpu = model2cpu(m);

# @save "saved_models_2/further_trained/ecoli_cgan_best.jld2" m_cpu




# fp = joinpath(regression_data_folder, fp_paths[3]);

# @load "saved_models_2/ecoli_cross_best_hp.jld2" hp

# @load "saved_models_2/ecoli_cross_best_hp_sparse.jld2" hp

# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=true);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=true);
# m_cpu = model2cpu(m);

# @save "saved_models_2/further_trained/ecoli_cross_best_sparse.jld2" m_cpu



# fp = joinpath(regression_data_folder, fp_paths[4]);

# @load "saved_models_2/ELF1_best_hp.jld2" hp

# @load "saved_models_2/ELF1_best_hp_sparse.jld2" hp
# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=true);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=true)
# m_cpu = model2cpu(m);
# @save "saved_models_2/further_trained/ELF1_best_sparse.jld2" m_cpu





# fp = joinpath(regression_data_folder, fp_paths[5]);
# @load "saved_models_2/HNF4A_best_hp.jld2" hp
# @load "saved_models_2/HNF4A_best_hp_sparse.jld2" hp
# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=false);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=false)
# m_cpu = model2cpu(m);
# @save "saved_models_2/further_trained/HNF4A_best.jld2" m_cpu


# fp = joinpath(regression_data_folder, fp_paths[6]);
# @load "saved_models_2/HNF1A_best_hp.jld2" hp
# @load "saved_models_2/HNF1A_best_hp_sparse.jld2" hp
# sd, data_load, data_load_v, data_load_t, dataload_w_v, dataload_all = get_data_loads(fp, hp)
# m = train_model(hp, sd, dataload_w_v; n_epochs=20, make_sparse=false);
# r2_ns, r_ns = plot_yy(data_load_t, m, hp; plotyy=false, make_sparse=false)
# m_cpu = model2cpu(m);
# @save "saved_models_2/further_trained/HNF1A_best.jld2" m_cpu





# df[df.f_name .== "ecoli_cgan",:].r2_s |> maximum
# df[df.f_name .== "ecoli_cgan",:].r2_ns |> maximum
# df_here = df[df.f_name .== "ecoli_cgan",:]
# max_ind = df_here.r2_s |>  argmax
# hp = df_here[max_ind,:].hp
# @save "saved_models_2/ecoli_cgan_best_hp_sparse.jld2" hp

# df[df.f_name .== "ecoli_cross",:].r2_s |> maximum
# df[df.f_name .== "ecoli_cross",:].r2_ns |> maximum
# df_here=df[df.f_name .== "ecoli_cross",:]
# max_ind = df_here.r2_s |>  argmax
# hp = df[max_ind,:].hp
# @save "saved_models_2/ecoli_cross_best_hp_sparse.jld2" hp

# df[df.f_name .== "ELF1",:].r2_s |> maximum
# df_here = df[df.f_name .== "ELF1",:]
# max_ind = df_here.r2_s |>  argmax
# hp = df_here[max_ind,:].hp
# @save "saved_models_2/ELF1_best_hp_sparse.jld2" hp

# df[df.f_name .== "HNF4A",:].r2_s |> maximum
# df_here = df[df.f_name .== "HNF4A",:]
# max_ind = df_here.r2_s |>  argmax
# hp = df_here[max_ind,:].hp
# @save "saved_models_2/HNF4A_best_hp_sparse.jld2" hp

# df_here = df[df.f_name .== "HNF4A",:]
# df[df.f_name .== "HNF1A",:].r2_s |> maximum
# max_ind = df[df.f_name .== "HNF1A",:].r2_s |>  argmax
# hp = df_here[max_ind,:].hp
# @save "saved_models_2/HNF1A_best_hp_sparse.jld2" hp


# inds = findall(df.r2_ns .== b1)
# df[inds, :]


result_tuple |> length

df = DataFrame(result_tuple)
findall(df.r2_s .≈ b1)
inds = findall(df.r_s .≈ b2)
df[inds, :]
df[findall(df.r2_s .≈ b1), :]
df[596,:]
df[596,:]


using DataFrames
df = DataFrame(result_tuple)



describe(df)

gdf = groupby(df, :f_name) 

for g in gdf
   describe(g)
end

gdf[1] |> describe
gdf[2] |> describe
gdf[3] |> describe
gdf[4] |> describe
gdf[5] |> describe
gdf[6] |> describe



i = 3
d1 = gdf[i].r_s .- gdf[i].r_ns
fd1 = filter(x->isfinite(x) && !isnan(x), d1)

using CairoMakie



name_map = Dict(
    "cnn_yeast" => 1, 
    "ecoli_cgan" => 2,
    "ecoli_cross" => 3,
    "ELF1" => 4,
    "HNF1A" => 5,
    "HNF4A" => 6
)

names_there = [
    "Yeast", 
    "E.coli cgan", 
    "E.coli cross", 
    "ELF1", 
    "HNF1A", 
    "HNF4A"
]

categories_1 = map(x-> name_map[x], df.f_name) 
categories_2 = map(x-> name_map[x], df.f_name) 

categories = vcat(categories_1, categories_2)

values_1 = df.r_ns
values_2 = df.r_s

values__1 = df.r2_ns    
values__2 = df.r2_s

vals = vcat(values_1, values_2)
vals_= vcat(values__1, values__2)
dodge = vcat(fill(1, length(values_1)), fill(2, length(values_2)))

values_keep_indices = map(x->isfinite(x) && !isnan(x) ? true : false, vals);
# values_keep_indices_ = map(x->isfinite(x) && !isnan(x) ? true : false, vals_);
categories = categories[values_keep_indices]
vals = vals[values_keep_indices]
dodge = dodge[values_keep_indices]
vals_ = vals_[values_keep_indices]

# Makie.scatter(0.001 .* randn(length(fd1)), fd1, markersize=5, color=:red)
# Makie.scatter(fill(1, length(fd1)), fd1, markersize=5, color=:red)


elems = [[MarkerElement(color = col, marker=:circle, markersize = 15,
          strokecolor = :black)] for col in ccolors]

# categories = rand(1:3, 1000)
# values = randn(1000)
# dodge = rand(1:2, 1000)

fig = Figure(size = (1000, 400), 
    );

ax = Axis(fig[1,1:2], 
    title = "Performance on the validation set",
    spinewidth = 0,
    titlegap = 10
    )
hidedecorations!(ax)

ax_r = Axis(fig[1, 1], 
    xticks = (1:6, names_there),
    ylabel = "Pearson correlation coefficient",
    xticklabelrotation=pi/6,
    )

ax_r2 = Axis(fig[1, 2], 
    xticks = (1:6, names_there),
    ylabel = "R squared",
    xticklabelrotation=pi/6,
    )

Makie.boxplot!(ax_r, categories, vals, dodge = dodge, 
    color = map(d->d==1 ? :orange : :lightblue, dodge),
    show_outliers = false,
    )

Makie.boxplot!(ax_r2, categories, vals_, dodge = dodge, 
    color = map(d->d==1 ? :orange : :lightblue, dodge),
    show_outliers = false,
    )

linkyaxes!(ax_r, ax_r2)

# Create marker elements for the legend
elems = [
    MarkerElement(color = :orange, marker = :circle, markersize = 15, strokecolor = :black),
    MarkerElement(color = :lightblue, marker = :circle, markersize = 15, strokecolor = :black)
];


# Add legend to the figure
axislegend(ax_r2, elems, ["normal filter weights", "sparse filter weights"]; position = :rb);

fig

save("hyperparam_performance.png", fig)




###### 

r2_bests = float_type.([0, 0, 0, 0, 0, 0])
r2_bests_sparse = float_type.([0, 0, 0, 0, 0, 0])
r_bests =  float_type.([0, 0, 0, 0, 0, 0])
r_bests_sparse =  float_type.([0, 0, 0, 0, 0, 0])

for r in result_tuple
    hp = r.hp
    f_name = r.f_name
    r2_ns = r.r2_ns
    r_ns = r.r_ns
    r2_s = r.r2_s
    r_s = r.r_s

    rand_fp_index = name_map[f_name]

    if r2_ns > r2_bests[rand_fp_index] && r_ns > r_bests[rand_fp_index]
        # m_cpu = model2cpu(m)
        @save "saved_models_2/$(f_name)_best.jld2" m_cpu
        @save "saved_models_2/$(f_name)_best_hp.jld2" hp
        r2_bests[rand_fp_index] = r2_ns
        r_bests[rand_fp_index] = r_ns 
    end

    if r2_s > r2_bests_sparse[rand_fp_index] && r_s > r_bests_sparse[rand_fp_index]
        # m_cpu = model2cpu(m)
        @save "saved_models_2/$(f_name)_best_sparse.jld2" m_cpu
        @save "saved_models_2/$(f_name)_best_hp_sparse.jld2" hp
        r2_bests_sparse[rand_fp_index] = r2_s
        r_bests_sparse[rand_fp_index] = r_s 
    end
end

[r_bests r_bests_sparse]

[r2_bests r2_bests_sparse]