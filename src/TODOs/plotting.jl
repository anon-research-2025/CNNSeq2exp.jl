
# using HDF5
# using Random
# using Flux
# using Zygote: @ignore
# using CUDA
# using StatsBase
# using ShaneGPUCountMinSketch

# using DataFrames, PlotPWM, CSV
# using Printf, JSON3, Mustache
# using CairoMakie, ImageShow
# using JLD2

# include("const.jl")
# include("model/hp1.jl")
# include("model/helper1.jl")
# include("model/model1.jl")
# include("model/forward.jl")
# include("load.jl")
# include("yy_plots.jl")

# include("contribution.jl")
# include("traverse.jl")

# include("motif_analysis/motif.jl")
# include("motif_analysis/helper.jl")
# include("render/const.jl")
# include("render/onehot2fasta.jl")
# include("render/logo_save.jl")
# include("render/fill_json.jl")
# include("render/template_css.jl")
# include("render/template_script.jl")
# include("render/templates_html.jl")
# include("render/boxplots.jl")
# include("render/render.jl")

# include("save_model.jl")

# include("main.jl")
# include("render/boxplot_publication.jl")

# include("plotting_for_pubs.jl")

# SAVE_WHERE = "pub_figs"


# TODO take these to be a small program



struct mixed_plt
    l_selected::Vector{String}
    c_selected::Vector{Vector{float_type}}
    logo_paths_mat::Matrix{String}
    contribs_mat::Matrix{Vector{float_type}}
    col_names_here::Vector{String}
    function mixed_plt(_m_, folder, l_selected, c_selected, keys_here, col_names_here; 
            take_how_many=5, )
        logo_paths_mat, contribs_mat = 
            get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=take_how_many)
        logo_paths_mat = reshape(logo_paths_mat, (take_how_many,3))[1:3,:]
        contribs_mat = reshape(contribs_mat, (take_how_many,3))[1:3,:]
        return new(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
    end
    function mixed_plt(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
        return new(l_selected, c_selected, logo_paths_mat, contribs_mat, col_names_here)
    end
end

len2rs = Dict(16=>(4,4), 15=>(5,3), 12=>(4,3), 10=>(2,5));
rs2figsize = Dict((4,4)=>(1200, 350), (5,3)=>(1200, 400), (2,5)=>(1200, 200), (4,3)=>(1200, 300));

# setting 2
len2rs = Dict(16=>(4,4), 15=>(3,5), 12=>(4,3), 10=>(2,5), 9 =>(3,3));
rs2figsize = Dict((4,4)=>(1200, 350), (3,5) => (1200,225), 
    (5,3)=>(1200, 400), (2,5)=>(1200, 200), (4,3)=>(1200, 300), (3,3)=>(800, 265));


# Given number of rows m, define f(z)
function linear_to_rc(z, m)
    row = mod1(z, m)
    col = div(z - 1, m) + 1
    return (row, col)
end

include("plotting_for_pubs.jl")
include("plotting_for_pubs_2.jl")

function get_singletons_matrix_fig_inner(logo_paths_mat, contribs_mat, solid_indices; 
    fig_size_test=(1200, 300), )
    fig = plot_singleton_logos_matrix(logo_paths_mat, contribs_mat, solid_indices; 
            fig_size=fig_size_test .* 2,
            padding = 21,
            jitter_width=0.15, 
            markersize=0.6,
            rain_start=6, 
            rain_end=12,
            c_and_r_gap=10,
            box_stroke_width=2,
            text_size=29,
        );
    return fig
end

function get_singletons_matrix_fig(_m_, folder, solid_indices)
    # load the singletons
    logo_paths_here_sorted, contribs_here_sorted, _ = 
        get_singleotons_paths_sorted(_m_, folder)
    rs, cs  = len2rs[length(logo_paths_here_sorted)]

    # map solid indices to matrix indices
    solid_indices_use = linear_to_rc.(solid_indices, rs)

    logo_paths_mat = reshape(logo_paths_here_sorted, (rs,cs))
    contribs_mat = reshape(contribs_here_sorted, (rs,cs))
    fig = get_singletons_matrix_fig_inner(logo_paths_mat, contribs_mat, solid_indices_use; 
        fig_size_test=rs2figsize[(rs,cs)])
    return fig, logo_paths_here_sorted, contribs_here_sorted
end

function get_plets_fig(_m_, folder, keys_here; take_how_many=15, fig_size=(1400, 1200))
    logo_paths, contribs = 
        get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=take_how_many)
    # determine the shape
    logo_paths_mat = reshape(logo_paths, (take_how_many,3))
    contribs_mat = reshape(contribs, (take_how_many,3))
    fig = plot_singleton_logos_matrix_plets(logo_paths_mat, contribs_mat; 
        fig_size=fig_size,
        markersize=[3,4,3])
    return fig
end


########################## yeast separate ##############################


@load "save_m/cnn_yeast_3.jld2" _m_
folder = "../saved_results/cnn_yeast_code_thresh_0.8_c_0.75__3"

solid_indices = collect(1:15)

fig, l_yeast, c_yeast = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "yeast_filters.png"), fig, px_per_unit=2.5)


####################################### yeast ################################################

@load "../reproduce_1/save_m/cnn_yeast_3.jld2" _m_
folder = "../saved_results_p/cnn_yeast_code_thresh_0.8_c_0.75__3"

solid_indices = [1, 15]

fig, l_yeast, c_yeast = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "yeast", "singletons_matrix.png"), fig, px_per_unit=1.25)

# get the keys of the plets
keys_here = [(m1 = "12", m2 = "12"), 
             (m1 = "12", m2 = "12", m3 = "12"), 
             (m1 = "4", m2 = "4")]

fig = get_plets_fig(_m_, folder, keys_here; take_how_many=15)

save(joinpath(SAVE_WHERE, "yeast", "xplets_matrix_full.png"), fig, px_per_unit=2.5)


# determine what are the extreme singletons
l_yeast_selected = [l_yeast[1], l_yeast[end]]
c_yeast_selected = [c_yeast[1], c_yeast[end]]

col_names_here = ["2 x filter 1 ", "3 x filter 1 ", "2 x filter 2"]

mp_yeast = mixed_plt(_m_, folder, 
    l_yeast_selected, c_yeast_selected, 
    keys_here, col_names_here; 
    take_how_many=5)

fig = plot_singleton_logos_extras_0(mp_yeast)

save(joinpath(SAVE_WHERE, "yeast", "brief_summary.png"), fig, px_per_unit=2)


################################ cross ################################################


@load "../reproduce_1/save_m/ecoli_cross_6.jld2" _m_
folder = "../saved_results_p/ecoli_cross_code_thresh_0.8_c_0.75__6"

solid_indices = [1, 16]

fig, l_cross, c_corss = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "cross", "singletons_matrix.png"), fig, px_per_unit=1.25)

# get the keys of the plets
keys_here = [(m1 = "11", m2 = "11"), 
             (m1 = "11", m2 = "11", m3 = "11"), 
             (m1 = "1", m2 = "1")]

fig = get_plets_fig(_m_, folder, keys_here; take_how_many=15)

save(joinpath(SAVE_WHERE, "cross", "xplets_matrix_full.png"), fig, px_per_unit=2.5)

l_cross_selected = [l_cross[1], l_cross[end]]
c_cross_selected = [c_corss[1], c_corss[end]]

col_names_here = ["2 x filter 1 ", "3 x filter 1 ", "2 x filter 2"]

mp_cross = mixed_plt(_m_, folder, 
    l_cross_selected, c_cross_selected, 
    keys_here, col_names_here; 
    take_how_many=5)

fig = plot_singleton_logos_extras_0(mp_cross)
save(joinpath(SAVE_WHERE, "cross", "brief_summary.png"), fig, px_per_unit=2)


############################ yeast and cross ############################
include("plotting_for_pubs_2.jl")

fig = plot_singleton_logos_extras_1([mp_yeast, mp_cross]; ABCs=["Yeast", "E.coli"])
save(joinpath(SAVE_WHERE, "additive_effect.png"), fig, px_per_unit=2)


####################################### ELF1 ################################################


@load "../reproduce_1/save_m/ELF1_deepTFBU_1.jld2" _m_
folder = "../saved_results_p/ELF1_deepTFBU_code_thresh_0.8_c_0.75__1"

# solid_indices = [1, 10, 14, 16]
solid_indices = [1, 10, 9, 8]

fig, l_elf, c_elf = 
    get_singletons_matrix_fig(_m_, folder, solid_indices);

save(joinpath(SAVE_WHERE, "elf1", "singletons_matrix.png"), fig, px_per_unit=1.25)

# get the keys of the plets
# keys_here = [(m1 = "12", m2 = "9"), 
#              (m1 = "6", m2 = "9"),
#              (m1 = "16", m2 = "9")]
keys_here = [(m1 = "12", m2 = "9"), 
             (m1 = "5", m2 = "9"),
             (m1 = "11", m2 = "9")]

fig = get_plets_fig(_m_, folder, keys_here; take_how_many=4, fig_size=(1300, 350))

save(joinpath(SAVE_WHERE, "elf1", "xplets_matrix_full.png"), fig, px_per_unit=2.5)


l_elf_selected = [l_elf[i] for i in solid_indices]
c_elf_selected = [c_elf[i] for i in solid_indices]

col_names_here = ["filter 1 and 2", "filter 3 and 2 ", "filter 4 and 2"]

mp_elf = mixed_plt(_m_, folder, 
    l_elf_selected, c_elf_selected, 
    keys_here, col_names_here; 
    take_how_many=4)

fig = plot_singleton_logos_extras_0(mp_elf; single_row_range=1:6, fig_size=(1400, 325) .* 1.5,
    markersize_plets=[3,3,3])

save(joinpath(SAVE_WHERE, "elf1", "brief_summary.png"), fig, px_per_unit=2)


####################################### HNF1A ################################################

@load "save_m/HNF1A_deepTFBU_1.jld2" _m_
folder = "../saved_results/HNF1A_deepTFBU_code_thresh_0.8_c_0.75__1"

solid_indices = [1, 2, 11, 15]

fig, l_hnf1a, c_hnf1a = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "hnf1a", "singletons_matrix.png"), fig, px_per_unit=1.25)


# get the keys of the plets
keys_here = [(m1 = "4", m2 = "10"),
             (m1 = "12", m2 = "11"),
             (m1 = "12", m2 = "12")]


fig = get_plets_fig(_m_, folder, keys_here; take_how_many=5, fig_size=(1300, 425))

save(joinpath(SAVE_WHERE, "hnf1a", "xplets_matrix_full.png"), fig, px_per_unit=2.5)
             
col_names_here = ["filter 2 and 1", "filter 4 and 3 ", "filter 4 and 4"]

l_hnf1a_selected = [l_hnf1a[i] for i in solid_indices]
c_hnf1a_selected = [c_hnf1a[i] for i in solid_indices]

mp_hnf1a = mixed_plt(_m_, folder, 
    l_hnf1a_selected, c_hnf1a_selected, 
    keys_here, col_names_here; 
    take_how_many=5)

fig = plot_singleton_logos_extras_0(mp_hnf1a; single_row_range=1:6, fig_size=(1400, 325) .* 1.5,
markersize_plets=[3,3,3])

save(joinpath(SAVE_WHERE, "hnf1a", "brief_summary.png"), fig, px_per_unit=2)


####################################### HNF4A ################################################

@load "../reproduce_1/save_m/HNF4A_deepTFBU_1.jld2" _m_
folder = "../save_results_selected/HNF4A_deepTFBU_code_thresh_0.8_c_0.75__1"


solid_indices = [1, 9, 5, 10]

fig, l_hnf4a, c_hnf4a = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "hnf4a", "singletons_matrix.png"), fig, px_per_unit=1.25)

# get the keys of the plets
keys_here = [(m1 = "2", m2 = "10"), 
             (m1 = "6", m2 = "10"), 
             (m1 = "6", m2 = "7")]

fig = get_plets_fig(_m_, folder, keys_here; take_how_many=4, fig_size=(1300, 400))

save(joinpath(SAVE_WHERE, "hnf4a", "xplets_matrix_full.png"), fig, px_per_unit=2.5)


l_hnf4a_selected = [l_hnf4a[i] for i in solid_indices]
c_hnf4a_selected = [c_hnf4a[i] for i in solid_indices]

col_names_here = ["filter 2 and 1", "filter 3 and 1 ", "filter 3 and 4"]

mp_hnf4a = mixed_plt(_m_, folder, 
    l_hnf4a_selected, c_hnf4a_selected, 
    keys_here, col_names_here; 
    take_how_many=4)

fig = plot_singleton_logos_extras_0(mp_hnf4a; single_row_range=1:5, 
    fig_size=(1400, 355) .* 1.5, 
    markersize_plets=[3,3,3])

save(joinpath(SAVE_WHERE, "hnf4a", "brief_summary.png"), fig, px_per_unit=2)


#### ############################ elf1 and hnf4a ############################

fig = plot_singleton_logos_extras_1([mp_elf, mp_hnf4a]; single_row_range=[1:6, 1:6], 
    fig_size=(1400, 2 * 315 + 10) .* 1.5,
    markersize_plets=[3,3,3],
    ax_scatter_x_range=(-1, 1), ABCs=["ELF1", "HNF4A"])

save(joinpath(SAVE_WHERE, "conditional_effect.png"), fig, px_per_unit=2)


####################################### HNF4A redo ################################################



@load "save_m_new/HNF4A_sparse.jld2" _m_
folder = "../saved_results_new_sparse/HNF4A_code_thresh_0.8_c_0.75"


solid_indices = [1, 2, 4]

fig, l_hnf4a, c_hnf4a = 
    get_singletons_matrix_fig(_m_, folder, solid_indices)

save(joinpath(SAVE_WHERE, "hnf4a", "singletons_matrix.png"), fig, px_per_unit=1.25)

# get the keys of the plets
keys_here = [(m1 = "4", m2 = "6"), 
             (m1 = "2", m2 = "6"), 
             (m1 = "6", m2 = "9")]

fig = get_plets_fig(_m_, folder, keys_here; take_how_many=10, fig_size=(1300, 800))

save(joinpath(SAVE_WHERE, "hnf4a", "xplets_matrix_full.png"), fig, px_per_unit=2.5)


# get the keys of the plets
# keys_here = [(m1 = "2", m2 = "6", m3 = "6"), 
#              (m1 = "4", m2 = "6", m3 = "6"), 
#              (m1 = "3", m2 = "7", m3 = "6")]



l_hnf4a_selected = [l_hnf4a[i] for i in solid_indices]
c_hnf4a_selected = [c_hnf4a[i] for i in solid_indices]


# get the keys of the plets
keys_here = [(m1 = "4", m2 = "6"), 
             (m1 = "2", m2 = "6"), 
             (m1 = "2", m2 = "6", m3 = "6")]

             take_how_many = 8

mp_hnf4a = mixed_plt(_m_, folder, 
    l_hnf4a_selected, c_hnf4a_selected, 
    keys_here, col_names_here; 
    take_how_many=take_how_many)

# adjust

logo_paths_mat, contribs_mat = 
    get_multi_paths_sorted(_m_, folder, keys_here; take_how_many=take_how_many)
logo_paths_mat = reshape(logo_paths_mat, (take_how_many,3))
contribs_mat = reshape(contribs_mat, (take_how_many,3))

i=2
logo_paths_mat[:,i]
median.(contribs_mat[:,i])

first_col_inds = [1,2,5]
second_col_inds = [2,5,8]
third_col_inds = [1,5,8]

logo_paths_mat = hcat(logo_paths_mat[first_col_inds,1], 
                      logo_paths_mat[second_col_inds,2],
                      logo_paths_mat[third_col_inds,3])

contribs_mat = hcat(contribs_mat[first_col_inds,1], 
                    contribs_mat[second_col_inds,2],
                    contribs_mat[third_col_inds,3])

col_names_here = ["filter 2 and 1", "filter 3 and 1 ", "filter 2, 1, and 1"]

mp_hnf4a = mixed_plt(l_hnf4a_selected, c_hnf4a_selected, 
    logo_paths_mat, contribs_mat, col_names_here)

fig = plot_singleton_logos_extras_0(mp_hnf4a; 
    single_row_range=1:5, fig_size=(1400, 315) .* 1.5, 
    markersize_plets=[3,3,4],
    ax_scatter_x_range=(-1, 1.2))

save(joinpath(SAVE_WHERE, "hnf4a", "brief_summary.png"), fig, px_per_unit=2)

