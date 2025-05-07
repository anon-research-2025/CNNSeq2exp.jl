
const float_type = Float32
const output_function = Flux.NNlib.tanh

# motif analysis
const config_count_threshold = 4 # a singleton motif must have at least this many counts
const UPTO = 4; # used in main.jl 
const MAX_CONFIG_SIZE = 200; # used in main.jl

const load_how_many = 50000
const convolve = Flux.NNlib.conv