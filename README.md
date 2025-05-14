# CNNSeq2exp

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://anon-research-2025.github.io/CNNSeq2exp.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://anon-research-2025.github.io/CNNSeq2exp.jl/dev/)
[![Build Status](https://github.com/anon-research-2025/CNNSeq2exp.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/anon-research-2025/CNNSeq2exp.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/anon-research-2025/CNNSeq2exp.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/anon-research-2025/CNNSeq2exp.jl)




The available processed datasets, along with the corresponding pretrained
models and hyperparameters (for example, those for the yeast dataset), are:
```
https://github.com/anon-research-2025/datasets/tree/main/processed/Aviv/yeast/yeast.hdf5 
https://github.com/anon-research-2025/saved_models/blob/main/aviv_yeast_model.jld2
https://github.com/anon-research-2025/saved_models/blob/main/aviv_yeast_hp.jld2
```

Once the necessary files are in place, set the following paths accordingly:
```
fp = <yeast hdf5 file path>
hp_load_path = <hyperparameter jld2 file path>
model_load_path = <model jld2 file path>
```
To use the pre-trained model and render the results, run:
```
save_where = "save_here"
take_pretrained_and_render(fp, hp_load_path, model_load_path, save_where)
```
Here, `save_where` specifies the directory where the rendered results will be saved.

If you prefer to train a new model instead of using a pre-trained one, simply run:
```
train_and_render(fp, save_where)
```

