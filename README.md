# CNNSeq2exp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anon-research-2025.github.io/CNNSeq2exp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anon-research-2025.github.io/CNNSeq2exp.jl/dev/)
[![Build Status](https://github.com/anon-research-2025/CNNSeq2exp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anon-research-2025/CNNSeq2exp.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anon-research-2025/CNNSeq2exp.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anon-research-2025/CNNSeq2exp.jl)




The available datasets that have been processed, the pretrained models, and
hyperparameters, e.g. for the yeast dataset, is
```
https://github.com/anon-research-2025/datasets/tree/main/processed/Aviv/yeast/yeast.hdf5 
https://github.com/anon-research-2025/saved_models/blob/main/aviv_yeast_model.jld2
https://github.com/anon-research-2025/saved_models/blob/main/aviv_yeast_hp.jld2
```

Once you have the files in place, and e.g. set 
```
fp = <yeast hdf5 file path>
hp_load_path = <hyperparameter jld2 file path>
model_load_path = <model jld2 file path>
```
You can then run the following command to use the pretrained model to render the results:
```
save_where = "save_here"
take_pretrained_and_render(fp, hp_load_path, model_load_path, save_where)
```
where `save_where` is the folder path that stores the rendered results.


