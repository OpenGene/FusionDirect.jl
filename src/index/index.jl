const DNA2BIT = Dict('A'=>0, 'T'=>1, 'C'=>2, 'G'=>3)

export Coord
type Coord
    contig::Int16
    pos::Int32

    Coord(contig, pos) = new(Int16(contig), Int32(pos))
end

# generate a coord represents a kmer key collision
function dupmark()
    return Coord(-1, 0)
end

export Index
type Index
    data::Dict{Int64, Coord}
    Index() = new(Dict{Int64, Coord}())
end

export kmer2key
function kmer2key(seq::Sequence)
    kmer2key(seq.seq)
end

# convert a kmer seq to a int64 key
function kmer2key(str::ASCIIString)
    str = uppercase(str)
    key = 0
    for c in str
        key = key*4 + DNA2BIT[c]
    end
    return key
end

export add
function add(index::Index, seq::Sequence, coord::Coord)
    key = kmer2key(seq)
    if key in keys(index.data)
        index.data[key] = dupmark()
        return false
    end
    index.data[key] = coord
    return true
end

export index_contig
function index_contig(index::Index, contig_seq::Sequence, contig_number::Int; kmer = 16)
    if kmer > 30
        error("kmer should be <= 30")
    end
    for i in 1:length(contig_seq)-kmer+1
        seq = contig_seq[i:i+kmer-1]
        add(index, seq, Coord(contig_number, i))
    end
end

export index_bed
# ref_folder is a folder contains fasta files by chr
# like chr1.fa, chr2.fa ...
function index_bed(ref_folder::AbstractString, bed_file::AbstractString)
    index = Index()
    io = open(bed_file)
    bed = readall(io)
    lines = split(bed, '\n')
    contig_number = 0
    for line in lines
        contig_number += 1
        line = rstrip(line, '\n')
        cols = split(line)
        println(cols)
        chr = ASCIIString(cols[1])
        from = parse(Int64, ASCIIString(cols[2]),)
        to = parse(Int64, ASCIIString(cols[3]))
        contig = ASCIIString(cols[3])
        chr_file = ref_folder * "/" * chr * ".fa"
        contig_seq = load(chr_file, chr, from, to)
        index_contig(index, contig_seq, contig_number)
    end
    return index
end