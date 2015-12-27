export load
function load(ref_file::ASCIIString, chr::ASCIIString, from::Int64, to::Int64)
    io = fasta_open(ref_file)
    while (fa = fasta_read(io))!=false
        if fa.name == chr
            return fa.sequence[from:to]
        end
    end
    return false
end