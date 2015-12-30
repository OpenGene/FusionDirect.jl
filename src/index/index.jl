const DNA2BIT = Dict('A'=>0, 'T'=>1, 'C'=>2, 'G'=>3)
const KMER = 16

type Coord
    # of which contig, usually a chr, a gene, an exon or any sequence
    contig::Int16
    # offset in this contig
    pos::Int32
    # in which strand, 1 means +, -1 means -
    strand::Int16

    Coord(contig, pos, strand = 1) = new(Int16(contig), Int32(pos), Int16(strand))
end

# generate a coord represents a kmer key collision
function dup_coord()
    return Coord(-1, 0)
end

# generate a coord represents a kmer key collision
function unknown_coord()
    return Coord(-2, 0)
end

# distance of two coords
function distance(c1::Coord, c2::Coord)
    if c1.contig != c2.contig || c1.strand != c2.strand
        return Inf
    else
        return c1.pos - c2.pos
    end
end

-(c1::Coord, c2::Coord) = distance(c1::Coord, c2::Coord)

is_dup(coord::Coord) = (coord.contig == -1)
is_unknown(coord::Coord) = (coord.contig == -2)
valid(coord::Coord) = (coord.contig >= 0)

type Index
    data::Dict{Int64, Coord}
    Index() = new(Dict{Int64, Coord}())
end

function kmer2key(seq::Sequence)
    kmer2key(seq.seq)
end

# convert a kmer seq to a int64 key
function kmer2key(str::ASCIIString)
    str = uppercase(str)
    if contains(str, "N")
        return -1
    end
    key = 0
    for c in str
        key = key*4 + DNA2BIT[c]
    end
    return key
end

# add and index a kmer
function add(index::Index, seq::Sequence, coord::Coord)
    key = kmer2key(seq)
    if key in keys(index.data)
        index.data[key] = dup_coord()
        return false
    end
    index.data[key] = coord
    return true
end

# add and index all kmers of a contig
function index_contig(index::Index, contig_seq::Sequence, contig_number::Int)
    len = length(contig_seq)
    for i in 1:len-KMER+1
        seq = contig_seq[i:i+KMER-1]
        add(index, seq, Coord(contig_number, i, 1))
    end
    reverse = ~contig_seq
    for i in 1:len-KMER+1
        seq = reverse[i:i+KMER-1]
        # here we align the position to the + strand for processing pair end sequence easier
        add(index, seq, Coord(contig_number, len-KMER+1, -1))
    end
end

# ref_folder is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function index_bed(ref_folder::AbstractString, bed_file::AbstractString)
    index = Index()
    io = open(bed_file)
    bed_file = readall(io)
    lines = split(bed_file, '\n')
    contig_number = 0
    bed = Dict{Int16, Dict{}}()
    ref = Dict{Int16, Sequence}()
    chr_bed = Dict{ASCIIString, Array{Int}}()
    for line in lines
        line = rstrip(line, '\n')
        cols = split(line)
        if length(cols)<4
            continue
        end
        contig_number += 1
        chr = ASCIIString(cols[1])
        if (chr in keys(chr_bed))==false
            chr_bed[chr]=Array{Int,1}()
        end
        push!(chr_bed[chr], contig_number)
        from = parse(Int64, ASCIIString(cols[2]))
        to = parse(Int64, ASCIIString(cols[3]))
        contig_name = ASCIIString(cols[4])
        bed[contig_number] = Dict("chr"=>chr, "name"=>contig_name, "from"=>from, "to"=>to)
    end
    for (chr,contig_numbers) in chr_bed
        chr_file = ref_folder * "/" * chr * ".fa"
        chr_seq = load_chr(chr_file, chr)
        if chr_seq==false
            error("cannot load data of chromosome $chr")
        end
        for contig_number in contig_numbers
            from = bed[contig_number]["from"]
            to = bed[contig_number]["to"]
            contig_seq = chr_seq[from:to]
            index_contig(index, contig_seq, contig_number)
            ref[contig_number] = contig_seq
        end
    end
    return index, bed, ref
end