
using Pkg

# Pkg.develop(PackageSpec(path = "/home/shane/dev/ShaneGPUCountMinSketch"))
# Pkg.add("ShaneGPUCountMinSketch")

using ShaneGPUCountMinSketch
using HDF5
using Random
using Flux
using Zygote: @ignore
using CUDA
using StatsBase
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
include("utils/load.jl")
include("utils/yy_plots.jl")

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

include("utils/save_model.jl")

include("main_subroutines.jl")


using Plots, LinearAlgebra, Distributions
include("linear_model/linear_model.jl")

regression_data_folder = "/home/shane/Desktop/seq2exp/seq2exp_datasets/processed_data_tanh"

SAVE_NAMEs = ["cnn_yeast", "ecoli_cgan", "ecoli_cross", "ELF1", "HNF4A", "HNF1A"]
fp_paths = [joinpath(regression_data_folder, "Aviv/yeast/yeast.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cgan_wanglab/ecoli_cgan.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cross_species_wanglab/ecoli_cross.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/ELF1/ELF1.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF4A/HNF4A.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF1A/HNF1A.hdf5")]

load_paths = ["saved_models/cnn_yeast_best_hp_sparse.jld2", 
              "saved_models/ecoli_cgan_best_hp_sparse.jld2", 
              "saved_models/ecoli_cross_best_hp_sparse.jld2",
              "saved_models/ELF1_best_hp_sparse.jld2",
              "saved_models/HNF4A_best_hp_sparse.jld2",
              "saved_models/HNF1A_best_hp_sparse.jld2"]

hp, m = train_model(fp);

# lms = get_linear_model(fp, hp, m); # test linear model

render_motifs(fp, hp , m, "save_there"; UPTO=4)