module CNNSeq2exp


using HDF5
using Random
using Flux
using Zygote: @ignore
using CUDA
using StatsBase
using DataFrames, CSV
using Printf, JSON3, Mustache
using CairoMakie, ImageShow
using JLD2, ShaneGPUCountMinSketch
using Plots, LinearAlgebra, Distributions
using EntroPlots


include("const.jl")
include("model/hp1.jl")
include("model/helper1.jl")
include("model/model1.jl")
include("model/forward.jl")
include("model/inference.jl")
include("utils/load.jl")
include("utils/yy_plots.jl")

include("motif_analysis/motif.jl")
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


export train_and_render, 
       take_pretrained_and_render,
       model, 
       HyperParameters,
       learned_pwms,
       learned_codeimg_filters


end
