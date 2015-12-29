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
                fusion_left, fusion_right = verify_fusion(index, bed, ref, pair)
                if fusion_left!=false && fusion_right!=false
                    name = ">"
                    name = name * bed[fusion_left.contig] * "-" * string(fusion_left.pos) * "-" 
                    name = name * bed[fusion_right.contig] * "-" * string(fusion_right.pos) * "-"
                    print(name, "R1\n",pair.read1.sequence.seq,"\n")
                    print(name, "R2\n",pair.read2.sequence.seq,"\n")
                end
            end
        end
    end
end

function verify_fusion(index, bed, ref, pair)
    offset, overlap_len, distance = overlap(pair)
    # this pair is overlapped, so merged it and segment the merged sequence
    if overlap_len>0 && distance<5
        seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len)
        seg, coords = segment(index, seq, ref)
        if length(seg) > 1
            return make_connected_fusion(index, ref, seg, coords)
        end
    else
        seg1, coords1 = segment(index, pair.read1.sequence, ref)
        if length(seg1) > 1
            return make_connected_fusion(index, ref, seg1, coords1)
        end
        seg2, coords2 = segment(index, pair.read2.sequence, ref)
        if length(seg2) > 1
            return make_connected_fusion(index, ref, seg2, coords2)
        end
        if length(seg1)==0 || length(seg2) == 0 || seg1[1][1] == 0 || seg2[1][1] == 0 
            return false, false
        end
        # fusion of different contigs on a pair
        if coords1[seg1[1][1]].contig != coords2[seg2[1][1]].contig
            l1 = seg1[1][1]
            l2 = seg2[1][1]
            len1 = length(coords1)
            len2= length(coords2)
            fusion_left = Coord(coords1[l1].contig, coords1[l1].pos + coords1[l1].strand * (len1 - l1))
            fusion_right = Coord(coords2[l2].contig, coords2[l2].pos + coords2[l2].strand * (len2 - l2))
            return fusion_left, fusion_right
        end
    end
    return false, false
end

function make_connected_fusion(index, ref, seg, coords)
    l1 = seg[1][1]
    r1 = seg[1][2]
    l2 = seg[2][1]
    r2 = seg[2][2]
    conjunct = div(r1+l2, 2)
    fusion_left = Coord(coords[l1].contig, coords[l1].pos + coords[l1].strand * (conjunct - l1))
    fusion_right = Coord(coords[l2].contig, coords[l2].pos + coords[l2].strand * (conjunct - l2))
    return fusion_left, fusion_right
end

# detect fusion and find the segmentations
function segment(index::Index, seq::Sequence, ref)
    counts, coords = stat(index, seq)
    clusters = cluster(coords)
    regions = Dict()
    for cluster in clusters
        # skip small clusters
        if length(cluster) < THRESHOLD
            continue
        end
        left, right = span(cluster, coords)
        regions[left]=right
    end
    # filter out those regions which are mostly inside other regions
    seg = filter_region(regions, coords)
    return seg, coords
end

# map to gene

# remove some inside regions
function filter_region(regions, coords)
    len = length(coords)
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
    # if some region is more or less a sub-region of the biggest region
    # then skip it
    for (l, r) in regions
        if r-right>10 || left-l>10
            clean_region[l]=r
        end
    end
    # sort the regions by the left
    sorted = sort(collect(clean_region))
    return sorted
end


# span the cluster on the ref to find the region
function span(cluster, coords)
    left = length(coords)
    right = 1
    for point in cluster
        left = min(left, point)
        right = max(right, point)
    end
    right = min(right+KMER-1, length(coords))
    return left, right
end

# clustering by distance
function cluster(coords::Array{Coord, 1})
    points = Set(1:length(coords))
    clusters = []
    while !isempty(points)
        # create a new cluster
        # get a valid seed
        seed = pop!(points)
        while(!valid(coords[seed]))
            seed = pop!(points)
        end
        cluster = Set(seed)
        # create a new working set
        working = Set(seed)
        while !isempty(working)
            current = pop!(working)
            for p in points
                # if it is not valid, just pop it out
                if !valid(coords[p])
                    pop!(points, p)
                elseif consistent(coords, current, p)
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
    coords = [unknown_coord() for i in 1:len]
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