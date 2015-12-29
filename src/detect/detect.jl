const THRESHOLD = 30

# ref_folder is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function detect(ref_folder::AbstractString, bed_file::AbstractString, r1fq::AbstractString, r2fq::AbstractString="")
    index, bed, ref = index_bed(ref_folder, bed_file)
    if r2fq != ""
        io = fastq_open_pair(r1fq, r2fq)
        i = 0
        while (pair = fastq_read_pair(io)) != false
            i += 1
            matches = detect_pair(index, pair)
            if length(matches)>1
                verified = verify(index, bed, ref, pair)
                if verified
                    name = ">$i-"
                    for m in matches
                        name = name * bed[m] * "-"
                    end
                    print(name, "R1\n",pair.read1.sequence.seq,"\n")
                    print(name, "R2\n",pair.read2.sequence.seq,"\n")
                end
            end
        end
    end
end

function verify(index, bed, ref, pair)
    offset, len, distance = overlap(pair)
    # this pair is overlapped, so merged it and segment it
    if len>0 && distance<3
        seq = simple_merge(pair.read1.sequence, pair.read2.sequence, offset, len)
        seg = segment(index, seq, ref)
        if length(seg) > 1
            return true
        end
    end
    return false
end

# detect fusion and find the segmentations
function segment(index::Index, seq::Sequence, ref)
    counts, coords = stat(index, seq)
    clusters = cluster(coords)
    regions = Dict()
    for cluster in clusters
        # skip small clusters
        if sizeof(cluster) < THRESHOLD
            continue
        end
        left, right = span(cluster)
        regions[left]=right
    end
    seg = clean_inside(regions)
end

# map to gene

# remove some inside regions
function clean_inside(regions)
    #find the biggest region
    left = 0
    right = 0
    for (l, r) in regions
        if r-l>right-left
            left = l
            right = r
        end
    end
    clean_region = Dict(left=>right)
    for (l, r) in regions
        if r-right>10 || left-l>10
            clean_region[l]=r
        end
    end
    sorted = sort(collect(clean_region))
    return sorted
end


# span the cluster on the ref to find the region
function span(cluster)
    left = Inf
    right = -Inf
    for point in cluster
        left = min(left, point)
        right = max(right, point)
    end
    return left, right
end

# clustering by distance
function cluster(coords::Array{Coord, 1})
    points = Set(1:length(coords))
    clusters = []
    while !isempty(points)
        # create a new cluster
        seed = pop!(points)
        cluster = Set(seed)
        # create a new working set
        working = Set(seed)
        while !isempty(working)
            current = pop!(working)
            for p in points
                if consistent(coords, current, p)
                    pop!(points, p)
                    push!(working, p)
                end
            end
            # current is done, push it to current cluster
            push!(cluster, current)
        end
        push!(clusters, cluster)
    end
    return clusters
end

function consistent(coords::Array{Coord, 1}, p1, p2)
    dis = coords[p2] - coords[p1]
    if dis > 1000
        return false
    end

    # the line of these two points be near 45 degree
    if abs(p2-p1)<40 && abs((p2-p1) - dis*coords[p1].strand) < 4
        return true
    end

    return false
end

# detect fusion in a pair of reads
function detect_pair(index::Index, pair::FastqPair)
    match1 = detect_read(index, pair.read1.sequence)
    match2 = detect_read(index, pair.read2.sequence)
    # merge them
    for m in match2
        if (m in match1) == false
            push!(match1, m)
        end
    end
    return match1
end

# simple fusion detection in a sequence
# this function is very fast and is not accurate
function detect_read(index::Index, seq::Sequence)
    counts, coords = stat(index, seq)
    matches = Array{Int16, 1}()
    for (k,v) in counts
        if v > THRESHOLD
            push!(matches, k)
        end
    end
    return matches
end

# stat the coordinations of kmers
function stat(index::Index, seq::Sequence)
    len = length(seq)
    coords = [unknown_coord() for i in 1:len-KMER+1]
    for i in 1:len-KMER+1
        kmer = seq[i:i+KMER-1]
        key = kmer2key(kmer)
        if key in keys(index.data)
            coords[i] = index.data[key]
        end
    end

    counts = Dict{Int16, Int}()
    for coord in coords
        if valid(coord)
            contig = coord.contig
            if contig in keys(counts)
                counts[contig] += 1
            else
                counts[contig] = 1
            end
        end
    end
    return counts, coords
end