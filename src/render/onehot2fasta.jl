
"""
Take an one-hot encoded sequence and convert it to ACGT string
"""
function one_hot_to_seq(one_hot::AbstractArray{T, 2}) where T <: Real
    nucleotides = "ACGT"
    seq = [nucleotides[argmax(@view one_hot[:, j])] for j in axes(one_hot, 2)]
    return join(seq)
end
"""
Take in a vector of (vector of two strings)
    1st string: header
    2nd string: sequence
and then make a fasta file
"""
function make_fasta_file(reads::Vector{T}, labels; 
    save_path="ok", save_name="seqs.fa"
    ) where T <: AbstractString
    mkpath(save_path)
    open(joinpath(save_path, save_name), "w") do fasta_file
        for (index, read) in enumerate(reads)
            println(fasta_file, "> $(round2(labels[index]))")
            println(fasta_file, read)
        end
    end
end

function make_fasta_at_save_path(onehotarr, labels; save_path="yo", save_name="seqs.fa")
    reads = Vector{String}(undef, (size(onehotarr, 4),))
    for n in axes(onehotarr,4)
        reads[n] = one_hot_to_seq(@view onehotarr[:,:,1,n])
    end    
    make_fasta_file(reads, labels; save_path=save_path, save_name=save_name)
end