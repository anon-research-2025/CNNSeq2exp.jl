
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

# TODO select which dataset to check
i = 1
# TODO setup the name of the folder to save the result
save_name = "yeast_result_test"

regression_data_folder = "datasets/processed"

fp_paths = [joinpath(regression_data_folder, "Aviv/yeast/yeast.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cgan_wanglab/ecoli_cgan.hdf5"), 
    joinpath(regression_data_folder, "gpro/ecoli_165_cross_species_wanglab/ecoli_cross.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/ELF1/ELF1.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF4A/HNF4A.hdf5"),
    joinpath(regression_data_folder, "deepTFBU/HNF1A/HNF1A.hdf5")]

hp_load_paths = ["saved_models/aviv_yeast_hp.jld2", 
                 "saved_models/ecoli_cgan_hp.jld2", 
                 "saved_models/ecoli_cross_hp.jld2",
                 "saved_models/ELF1_hp.jld2",
                 "saved_models/HNF4A_hp.jld2",
                 "saved_models/HNF1A_hp.jld2"]

model_load_paths = ["saved_models/aviv_yeast_model.jld2", 
                    "saved_models/ecoli_cgan_model.jld2", 
                    "saved_models/ecoli_cross_model.jld2",
                    "saved_models/ELF1_model.jld2",
                    "saved_models/HNF4A_model.jld2",
                    "saved_models/HNF1A_model.jld2"]

fp, hp_load_path, model_load_path = 
    fp_paths[i], hp_load_paths[i], model_load_paths[i]


@load hp_load_path hp
@load model_load_path m_cpu

m = model2gpu(m_cpu);

lms = get_linear_model(fp, hp, m); # test linear model

render_motifs(fp, hp , m, save_name; UPTO=4)