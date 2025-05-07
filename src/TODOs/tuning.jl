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

const float_type = Float32
const load_how_many = 30000
const convolve = Flux.NNlib.conv

include("model/hp1.jl")
include("model/helper1.jl")
include("model/model1.jl")
include("model/forward.jl")
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
# fp = joinpath(regression_data_folder, "gpro/yeast_110_evolution_aviv/yeast.hdf5")
fp = joinpath(regression_data_folder, "gpro/ecoli_165_cross_species_wanglab/ecoli_cross.hdf5")
fp = joinpath(regression_data_folder, "gpro/ecoli_165_cgan_wanglab/ecoli_cgan.hdf5")

fp = joinpath(regression_data_folder, "chaolin_psi_1/chaolin_psi_1.hdf5")
fp = joinpath(regression_data_folder, "gpro/ecoli_50_wgan_diffusion_wanglab/ecoli_wgan.hdf5")

fp = joinpath(regression_data_folder, "gpro/ecoli_165_cgan_wanglab/ecoli_cgan.hdf5")


fp = joinpath(regression_data_folder, "Aviv/yeast/yeast.hdf5")

fp = joinpath(regression_data_folder, "deepTFBU/HNF1A/HNF1A.hdf5")

fp = joinpath(regression_data_folder, "deepTFBU/HNF4A/HNF4A.hdf5")

fp = joinpath(regression_data_folder, "deepTFBU/ELF1/ELF1.hdf5")


save_name_here = "cnn_cross"

result_tuples = []

r2_best, r_best = 0,0 

# @load "saved_models/cnn_cross_gpro_best_hp.jld2" hp

# hp.batch_size = 32

# hp = HyperParameters()

for _ = 1:15
    hp = random_hyperparameters()
    sd, data_load, data_load_v, data_load_t = load_data(hp, fp)
    while get_last_embedding_len(data_load, hp) < 1 ||  get_last_embedding_len(data_load, hp) > 5
        hp = random_hyperparameters()
    end
    @info "last embedding len $(get_last_embedding_len(data_load, hp))"
    m = train_model(hp, sd, data_load; n_epochs=48, make_sparse=true);
    r2, r = plot_yy(data_load_v, m, hp; plotyy=false, make_sparse=true)

    # m_ = train_model(hp, sd, data_load; make_sparse=true, n_epochs=16);
    # r2_, r_ = plot_yy(data_load_v, m_, hp; plotyy=false)
    # @info "r2: $(r2), r: $(r), r2_s: $(r2_), r_s: $(r_)"

    if r2 > r2_best && r > r_best
        m_cpu = model2cpu(m)
        @save "saved_models/$(save_name_here)_best.jld2" m_cpu
        @save "saved_models/$(save_name_here)_best_hp.jld2" hp
        r2_best, r_best = r2, r
    end
end




@save "result_tuples_yeast.jld2" result_tuples



# m_cpu = model2cpu(m)
# @save "saved_models/cnn_yeast_gpro.jld2" m_cpu
# @save "saved_models/cnn_yeast_gpro_hp.jld2" hp
# m_gpu = model2gpu(m_cpu)

result_tuples |> length

r2_means = map(x->mean(y.r2 for y in x), result_tuples)
argmax_ind = argmax(r2_means)
r2_means[argmax(r2_means)]

r_means = map(x->mean(y.r for y in x), result_tuples)
r_means[isnan.(r_means)] .= 0
argmax_ind = argmax(r_means)
r_means[argmax_ind]

