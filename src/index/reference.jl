
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

function load_single_file_reference(ref_file::AbstractString, ref::Dict{AbstractString, Sequence})
    io = fasta_open(ref_file)
    while (fa = fasta_read(io))!=false
        ref[fa.name] = fa.sequence
    end
end

function load_folder_reference(ref_folder::AbstractString, ref::Dict{AbstractString, Sequence})
    files = readdir(ref_folder)
    for file in files
        if contains(file, ".fa") || contains(file, ".fa.gz")
            load_single_file_reference(joinpath(ref_folder, file), ref)
        end
    end
end

# load the whold reference to a Dict{AbstractString, Sequence}
# the ref_path can be either a single reference fa file, or a folder contains reference files
# all reference files can be gzipped
function load_reference(ref_path::AbstractString)
    ref = Dict{AbstractString, Sequence}()
    if isdir(ref_path)
        load_folder_reference(ref_path, ref)
    else
        load_single_file_reference(ref_path, ref)
    end
    return ref
end