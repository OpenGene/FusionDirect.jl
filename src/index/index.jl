# if index format is changed, than increase the format version number
const INDEX_FORMAT_VER = "v2"
const DNA2BIT = Dict('A'=>0, 'T'=>1, 'C'=>2, 'G'=>3)
const KMER = 16
const REF_SAMPLE_STEP = 4

"""
the position of a k-mer sequence is the position of its first base, normalized to the direction of strand +
this is why a kmer sequence and its reverse compelement has a position offset of KMER-1
for example
contig: AAAATTTTCCCCGGG[G]TCATGATTACCAACCAATACCGTGGGATGG
position of the first 16-mer AAAATTTTCCCCGGGG is: 1
position of reverse compelement of the first 16-mer, CCCCGGGGAAAATTTT is: 1+16-1 = 16, the position bracketed
"""
immutable Coord
    # of which contig, usually a chr, a gene, an exon or any sequence
    contig::Int16
    # offset in this contig, the sign of pos is the strand
    pos::Int32

    Coord(contig, pos) = new(Int16(contig), Int32(pos))
    Coord(contig, pos, strand) = new(Int16(contig), Int32(strand * pos))
end

# for debugging
function display_coords(coords::Array{Coord, 1})
    i = 0
    for coord in coords
        if valid(coord)
            print("(", coord.pos>0?"+":"", coord.contig, ":", coord.pos, ")\t")
        elseif coord.contig == -1
            print("(duplicate)\t")
        elseif coord.contig == -2
            print("(unknown)\t")
        end
        i+=1
        if i%10 == 0
            print("\n")
        end
    end
    print("\n")
end

# generate a coord represents a kmer key collision
function dup_coord()
    return Coord(-1, 0)
end

# generate a coord represents a unknown key
function unknown_coord()
    return Coord(-2, 0)
end

# distance of two coords
function distance(c1::Coord, c2::Coord)
    if c1.contig != c2.contig || sign(c1.pos) != sign(c2.pos)
        return Inf
    else
        return abs(c1.pos) - abs(c2.pos)
    end
end

# abs distance of two coords
function abs_distance(c1::Coord, c2::Coord)
    return abs(distance(abs(c1), abs(c2)))
end

-(c1::Coord, c2::Coord) = distance(c1::Coord, c2::Coord)
abs(c::Coord) = Coord(c.contig, abs(c.pos))

is_dup(coord::Coord) = (coord.contig == -1)
is_unknown(coord::Coord) = (coord.contig == -2)
valid(coord::Coord) = (coord.contig >= 0)

typealias KmerCoord Dict{Int64, Coord}
typealias KmerCoordList Dict{Int64, Array{Coord, 1}}

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
    try
        for c in str
            key = key*4 + DNA2BIT[c]
        end
    catch e
        info(e)
        error("Invalid k-mer $str")
    end
    return key
end

# add and index a kmer of a panel
function add_to_panel_index(kmer_coord::KmerCoord, seq::Sequence, coord::Coord)
    key = kmer2key(seq)
    if haskey(kmer_coord, key)
        kmer_coord[key] = dup_coord()
        return false
    end
    kmer_coord[key] = coord
    return true
end


# add and index all kmers of a contig
function index_contig(kmer_coord::KmerCoord, contig_seq::Sequence, contig_number::Int)
    len = length(contig_seq)
    for i in 1:len-KMER+1
        seq = contig_seq[i:i+KMER-1]
        add_to_panel_index(kmer_coord, seq, Coord(contig_number, i))
        add_to_panel_index(kmer_coord, ~seq, Coord(contig_number, - (i+KMER-1)))
    end
end

function load_bed(bed_file::AbstractString)
    io = open(bed_file)
    bed_file = readall(io)
    lines = split(bed_file, '\n')
    contig_number = 0
    panel = Dict{Int16, Dict{}}()
    chr_contigs = Dict{ASCIIString, Array{Int}}()
    for line in lines
        line = rstrip(line, '\n')
        cols = split(line)
        if length(cols)<4
            continue
        end
        contig_number += 1
        chr = ASCIIString(cols[1])
        if !haskey(chr_contigs, chr)
            chr_contigs[chr]=Array{Int,1}()
        end
        push!(chr_contigs[chr], contig_number)
        from = parse(Int64, ASCIIString(cols[2]))
        to = parse(Int64, ASCIIString(cols[3]))
        contig_name = ASCIIString(cols[4])
        panel[contig_number] = Dict("chr"=>chr, "name"=>contig_name, "from"=>from, "to"=>to)
    end
    return panel, chr_contigs
end

function get_ref_files(ref_path::AbstractString)
    ref_files = []
    # if ref_path is a folder, then find all fa files in it
    if isdir(ref_path)
        files = readdir(ref_path)
        for file in files
            if endswith(file, ".fa")
                push!(ref_files, joinpath(ref_path, file))
            end
        end
    else
        push!(ref_files, ref_path)
    end
    return ref_files
end

function index_chr_bed(chrfa::FastaRead, panel_kmer_coord::KmerCoord, chr_contigs, panel, panel_seq)
    contig_numbers = chr_contigs[chrfa.name]
    chr_seq = chrfa.sequence
    for contig_number in contig_numbers
        from = panel[contig_number]["from"]
        to = panel[contig_number]["to"]
        contig_seq = chr_seq[from:to]
        index_contig(panel_kmer_coord, contig_seq, contig_number)
        panel_seq[contig_number] = contig_seq
    end
end

# ref_path is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function make_panel_index(ref_path::AbstractString, bed_file::AbstractString)
    panel_kmer_coord = KmerCoord()
    panel_seq = Dict{Int16, Sequence}()
    panel, chr_contigs = load_bed(bed_file)
    ref_files = get_ref_files(ref_path)

    ref = []
    # for each reference file, get its chromosome and index them
    for ref_file in ref_files
        io = fasta_open(ref_file)
        pos = position(io)
        while (chrfa = fasta_read(io))!=false
            if haskey(chr_contigs, chrfa.name)
                index_chr_bed(chrfa, panel_kmer_coord, chr_contigs, panel, panel_seq)
                info("index_chr_bed " * chrfa.name)
            end
            push!(ref, Dict("file"=>ref_file, "position"=>pos, "length"=>length(chrfa)))
            pos = position(io)
        end
    end

    ref_kmer_coords = make_kmer_coord_list(ref, panel_kmer_coord)
    return Dict("panel"=>panel, "seq"=>panel_seq, "kmer_coord"=>panel_kmer_coord, "ref_kmer_coords"=>ref_kmer_coords)
end

function make_kmer_shared_array(panel_kmer_coord::KmerCoord)
    kmer_buffer_size = div(4^KMER, 8)
    shared_kmer = SharedArray(UInt8, kmer_buffer_size, init=true)
    for k in keys(panel_kmer_coord)
        offset = div(k, 8) + 1
        shift = k % 8
        shared_kmer[offset] |= (0b00000001 << shift)
    end
    return shared_kmer
end

# make an index with the reference data and a panel
# for each kmer of the panel, create an array, and store the coordinations of same kmer
function make_kmer_coord_list(ref, panel_kmer_coord::KmerCoord)
    ref_index = KmerCoordList()
    total = 0
    for k in keys(panel_kmer_coord)
        ref_index[k]=Array{Coord, 1}()
    end

    shared_kmer = make_kmer_shared_array(panel_kmer_coord)

    # create parallel tasks
    tasks = []
    for chrid in 1:length(ref)
        chrinfo = ref[chrid]
        task = Dict("chrid"=>chrid, "chrinfo"=>chrinfo, "shared_kmer"=>shared_kmer)
        push!(tasks, task)
    end
    if length(workers()) <= 1
        #addprocs()
    end

    # run parallel for indexing
    info("tasks:" * string(length(tasks)))
    info("start worker processes")
    results = pmap(make_kmer_coord_list_chr, tasks)
    #results = []
    #for task in tasks
    #    push!(results, make_kmer_coord_list_chr(task))
    #end

    # merge the result index
    i = 0
    for result in results
        i += 1
        info("merging " * string(i))
        if isa(result, RemoteException)
            info(result)
            continue
        end
        for (k, v) in result
            append!(ref_index[k], v)
        end
        gc()
    end

    info("merge done")

    # destroy worker processes
    worker_procs = workers()
    if length(worker_procs)>1 || (length(worker_procs) == 1 && worker_procs[1]!=1)
        rmprocs(worker_procs)
    end

    return ref_index
end

function inkmer(shared_kmer::SharedArray{UInt8, 1}, k::Int64)
    offset = div(k, 8) + 1
    shift = k % 8 
    return (shared_kmer[offset] & (0b00000001 << shift))>0
end

# run in parallel
# make an index with a chromosome of the reference data and a panel
function make_kmer_coord_list_chr(task)
    chrid=task["chrid"]
    chrinfo=task["chrinfo"]
    io = fasta_open(chrinfo["file"])
    chrname = basename(chrinfo["file"])
    seek(io, chrinfo["position"])
    len = chrinfo["length"]
    info("indexing $chrid:" * chrinfo["file"])
    fa = fasta_read(io)
    chrseq = fa.sequence
    shared_kmer=task["shared_kmer"]
    ref_index = KmerCoordList()
    total = 0
    for i in 1:REF_SAMPLE_STEP:len-KMER+1
        if i%10000000 < REF_SAMPLE_STEP
            info("$chrid-$chrname:$i/$len:$total")
        end
        seq = chrseq[i:i+KMER-1]
        key = kmer2key(seq)
        if inkmer(shared_kmer, key)
            if !haskey(ref_index, key)
                ref_index[key]=Array{Coord, 1}()
            end
            push!(ref_index[key], Coord(Int16(chrid), i))
            total+=1
        end
        key = kmer2key(~seq)
        if inkmer(shared_kmer, key)
            if !haskey(ref_index, key)
                ref_index[key]=Array{Coord, 1}()
            end
            push!(ref_index[key], Coord(Int16(chrid), -(i+KMER-1)))
            total+=1
        end
    end
    gc()
    info("finished $chrid:" * chrinfo["file"])
    return ref_index
end

function get_cache_path(ref_path::AbstractString, bed_file::AbstractString)
    ref_name = basename(ref_path) * "." * string(filesize(ref_path))
    bed_name = basename(bed_file) * "." * string(filesize(bed_file))
    cache_name = bed_name * "_" * ref_name * "_" *  INDEX_FORMAT_VER * ".idx"
    cache_path = joinpath(dirname(bed_file), cache_name)
    return cache_path
end

function index_bed(ref_path::AbstractString, bed_file::AbstractString)
    cache_path = get_cache_path(ref_path, bed_file)
    loaded = false
    # load the index from a cache file
    if isfile(cache_path)
        try
            io = open(cache_path)
            index = deserialize(io)
            if !isa(index, Dict)
                error("invalid index file, index should be a serialized Dict")
            elseif !haskey(index, "panel") || !haskey(index, "seq") || !haskey(index, "kmer_coord") || !haskey(index, "ref_kmer_coords")
                error("invalid index file, wrong Dict keys")
            end

            loaded = true
            return index
        catch(e)
            warn(e)
            warn("Failed to load the pre-built index ($cache_path), maybe it is not completed when saving it. Attemping to delete it now!")
            try
                rm(cache_path)
            catch(e)
                warn(e)
                warn("Failed to delete $cache_path, please do it manually")
                warn("FusionDirect exited")
                exit(-1)
            end
        end
    end

    if !loaded
        info("Index doesn't exist, indexing now, it may take several minutes to a few hours")
        info("After the index is created, loading it will be very fast")
        index = make_panel_index(ref_path, bed_file)
        # save the index to a cache file
        try
            io = open(cache_path, "w")
            serialize(io, index)
        catch(e)
            warn("Unable to save index to $cache_path. Please check if the folder is writable for current user, if the index is not cached, it will be built every time you run FusionDirect, which is slow.")
        end
        return index
    end
end