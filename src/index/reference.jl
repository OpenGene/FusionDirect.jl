
# load chromosome reference sequence as an OpenGene::Sequence
function load_chr(ref_file::AbstractString, chr::AbstractString)
    io = fasta_open(ref_file)
    while (fa = fasta_read(io))!=false
        if fa.name == chr
            return fa.sequence
        end
    end
    return false
end