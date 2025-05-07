
"""
Get the number of sequence segments and heads from the hdf5 file.
    fid: hdf5 file, obtain via, e.g., 
        h5open(<filepath>, "r")
"""
function get_num_seq_segments_and_heads(fid::HDF5.File)
    num_seq_segments = 
        [i[1] .== "seq" for i in split.(collect(keys(fid)), "_")] |> sum
    @assert num_seq_segments > 0 "There must be sequences in the hdf5 file"
    num_heads = fid["expr/train"] |> size |> first
    @assert num_heads > 0 "There must be at least one y value"
    return num_seq_segments, num_heads
end

"""
fid: hdf5 file
seq_len: length of the sequence (assume for now that all sequences have the same length)
num_seq_segments: number of sequence segments (number of input sequences)
num_heads: number of heads (number of differny y values)
num_seqs: number of sequences in the dataset
permute_indices_train: indices for permuting the sequence (for training)
"""
mutable struct sequence_data
    fid::HDF5.File
    seq_len::Int
    num_seq_segments::Int
    num_heads::Int
    num_seqs_train::Int
    num_seqs_valid::Int
    num_seqs_test::Int
    permuted_indices_train::Vector{Vector{Int}}
    function sequence_data(hdf5_filepath::String; load_how_many=load_how_many)
        fid = h5open(hdf5_filepath, "r")
        num_seq_segments, num_heads = get_num_seq_segments_and_heads(fid)
        seq_len = size(fid["seq_1/train"],2)
        num_seqs_train = fid["seq_1/train"] |> size |> last
        num_seqs_valid = fid["seq_1/valid"] |> size |> last
        num_seqs_test  = fid["seq_1/test"] |> size |> last
        permuted_indices = shuffle(1:num_seqs_train)
        permuted_indices_train = 
            [permuted_indices[i:min(i+load_how_many-1, num_seqs_train)] 
                for i in 1:load_how_many:num_seqs_train]
        new(fid, 
            seq_len,
            num_seq_segments, 
            num_heads, 
            num_seqs_train, 
            num_seqs_valid, 
            num_seqs_test, 
            permuted_indices_train)
    end
end

close_hdf5!(sd::sequence_data) = close(sd.fid)


"""
Make new permuted indices for the training data.
sd: sequence_data
load_how_many: how many data entries to load each time
"""
function make_new_permuted_indices!(sd::sequence_data; 
    load_how_many=load_how_many)
    permuted_indices = shuffle(1:sd.num_seqs_train)
    sd.permuted_indices_train = 
        [permuted_indices[i:min(i+load_how_many-1, sd.num_seqs_train)] 
            for i in 1:load_how_many:sd.num_seqs_train]
end

"""
Get the number of mega batches for training.
    num_mega_batches: number of mega batches
    sd: sequence_data
"""
num_mega_batches_train(sd::sequence_data) = length(sd.permuted_indices_train)

function obtain_dataload_ith_mega_batch(
    sd::sequence_data, ith_mega_batch::Int, set_name="train";
    batch_size=128, shuffle=false
)
    if set_name == "train"
        permuted_indices = sd.permuted_indices_train[ith_mega_batch]
    elseif set_name == "valid"        
        permuted_indices = collect(1:sd.num_seqs_valid) 
    elseif set_name == "test"
        permuted_indices = collect(1:sd.num_seqs_test) 
    else
        error("set_name must be either 'train', 'valid', or 'test'")
    end

    sequences = [Array{float_type,4}(undef, 
        (4, sd.seq_len, 1, length(permuted_indices))) for _ = 1:sd.num_seq_segments]

    @inbounds for (j,p) in enumerate(permuted_indices)
        for i in 1:sd.num_seq_segments
            sequences[i][:,:,1,j] = 
                reshape(sd.fid["seq_$(i)/$(set_name)"][:,:,p], (4, sd.seq_len, 1))
        end
    end
    
    labels = Array{float_type, 2}(undef,
        (sd.num_heads, length(permuted_indices)))
    @inbounds for (j,p) in enumerate(permuted_indices)
        labels[:,j] = sd.fid["expr/$(set_name)"][:,p]
    end

    data_load = Flux.DataLoader(
        (sequences..., labels);
        batchsize = batch_size,
        shuffle = shuffle,
        partial = false,
    )
    return data_load
end
