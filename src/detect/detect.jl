const THRESHOLD = 30

# ref_path is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function detect(ref_path::AbstractString, bed_file::AbstractString, r1fq::AbstractString, r2fq::AbstractString="")
    index = index_bed(ref_path, bed_file)
    panel = index["panel"]
    panel_seq = index["seq"]
    panel_kmer_coord = index["kmer_coord"]
    # pair end sequencing
    if r2fq != ""
        io = fastq_open_pair(r1fq, r2fq)
        fusion_pairs = Dict()
        while (pair = fastq_read_pair(io)) != false
            matches = detect_pair(panel_kmer_coord, pair)
            if length(matches)>1
                fusion_left, fusion_right = verify_fusion_pair(panel_kmer_coord, panel, panel_seq, pair)
                if fusion_left!=false && fusion_right!=false && distance(fusion_left, fusion_right)>1000
                    add_to_fusion_pair(fusion_pairs, fusion_left, fusion_right, pair)
                end
            end
        end
        display_fusion_pair(fusion_pairs, panel_seq, panel)
    end
end

function display_fusion_pair(fusion_pairs, panel_seq, panel)
    for (fusion_key, fusion_reads) in fusion_pairs
        contig1, contig2 = fusion_key
        name1 = panel[contig1]["name"]
        name2 = panel[contig2]["name"]
        # give the fusion as a fasta comment line
        print("#Fusion:", name1, "-", name2, "\n")
        # display all reads support this fusion
        for reads in fusion_reads
            fusion_left, fusion_right, pair = reads
            name = ">"
            name = name * panel[fusion_left.contig]["name"] * "|" * strand_name(fusion_left.strand) * "|" * coord_to_chr(fusion_left, panel) * "_"
            name = name * panel[fusion_right.contig]["name"] * "|" * strand_name(fusion_right.strand)  * "|" * coord_to_chr(fusion_right, panel) * "/"
            print(name, "1\n",pair.read1.sequence.seq,"\n")
            print(name, "2\n",pair.read2.sequence.seq,"\n")
        end
    end
end

function strand_name(strand)
    return strand>0?"+":"-"
end

function coord_to_chr(coord, panel)
    from = panel[coord.contig]["from"]
    chr = panel[coord.contig]["chr"]
    pos = from + coord.pos 
    return chr * ":" * string(pos)
end

function add_to_fusion_pair(fusion_pairs, fusion_left, fusion_right, pair)
    key = (min(fusion_left.contig, fusion_right.contig), max(fusion_left.contig, fusion_right.contig))
    if (key in keys(fusion_pairs)) == false
        fusion_pairs[key]=[]
    end
    push!(fusion_pairs[key], (fusion_left, fusion_right, pair))
end

function verify_fusion_se(panel_kmer_coord, panel, panel_seq, sequence)
    seg1, coords1 = segment(panel_kmer_coord, sequence, panel_seq)
    if length(seg1) > 1
        return make_connected_fusion(panel_kmer_coord, panel_seq, seg1, coords1)
    end
end

function verify_fusion_pair(panel_kmer_coord, panel, panel_seq, pair)
    offset, overlap_len, distance = overlap(pair)
    # this pair is overlapped, so merged it and segment the merged sequence
    if overlap_len>0 && distance<5
        seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len)
        seg, coords = segment(panel_kmer_coord, seq, panel_seq)
        if length(seg) > 1
            return make_connected_fusion(panel_kmer_coord, panel_seq, seg, coords)
        end
    else
        seg1, coords1 = segment(panel_kmer_coord, pair.read1.sequence, panel_seq)
        if length(seg1) > 1
            return make_connected_fusion(panel_kmer_coord, panel_seq, seg1, coords1)
        end
        seg2, coords2 = segment(panel_kmer_coord, pair.read2.sequence, panel_seq)
        if length(seg2) > 1
            return make_connected_fusion(panel_kmer_coord, panel_seq, seg2, coords2)
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
            fusion_left = Coord(coords1[l1].contig, coords1[l1].pos + coords1[l1].strand * (len1 - l1), coords1[l1].strand)
            # reverse the strand of fusion right to align the pair in same direction
            fusion_right = Coord(coords2[l2].contig, coords2[l2].pos + coords2[l2].strand * (len2 - l2), coords2[l2].strand * -1)

            return fusion_left, fusion_right
        end
    end
    return false, false
end

function make_connected_fusion(panel_kmer_coord, panel_seq, seg, coords)
    l1 = seg[1][1]
    r1 = seg[1][2]
    l2 = seg[2][1]
    r2 = seg[2][2]
    conjunct = div(r1+l2, 2)
    #display_coords(coords)
    fusion_left = Coord(coords[l1].contig, coords[l1].pos + coords[l1].strand * (conjunct - l1), coords[l1].strand)
    fusion_right = Coord(coords[l2].contig, coords[l2].pos + coords[l2].strand * (conjunct - l2), coords[l2].strand)
    return fusion_left, fusion_right
end

# detect fusion and find the segmentations
function segment(panel_kmer_coord::KmerCoord, seq::Sequence, panel_seq)
    counts, coords = stat(panel_kmer_coord, seq)
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


# span the cluster on the panel_seq to find the region
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
        while(!valid(coords[seed]) && !isempty(points))
            seed = pop!(points)
        end
        # we don't find a valid seed, so break
        if !valid(coords[seed])
            break
        end
        # create a new cluster
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
function detect_pair(panel_kmer_coord::KmerCoord, pair::FastqPair)
    match1 = detect_read(panel_kmer_coord, pair.read1.sequence)
    match2 = detect_read(panel_kmer_coord, pair.read2.sequence)
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
function detect_read(panel_kmer_coord::KmerCoord, seq::Sequence)
    counts, coords = stat(panel_kmer_coord, seq)
    matches = Array{Int16, 1}()
    for (k,v) in counts
        if v > THRESHOLD
            push!(matches, k)
        end
    end
    return matches
end

# get the coordinations on the panel for given sequence
function stat(panel_kmer_coord::KmerCoord, seq::Sequence)
    len = length(seq)
    coords = [unknown_coord() for i in 1:len]
    for i in 1:len-KMER+1
        kmer = seq[i:i+KMER-1]
        key = kmer2key(kmer)
        if haskey(panel_kmer_coord, key)
            coords[i] = panel_kmer_coord[key]
        end
    end

    counts = Dict{Int16, Int}()
    for coord in coords
        if valid(coord)
            contig = coord.contig
            if haskey(counts, contig)
                counts[contig] += 1
            else
                counts[contig] = 1
            end
        end
    end
    return counts, coords
end

# get the coordinations on the whole reference for given sequence
function stat_ref(ref_kmer_coords::KmerCoordList, seq::Sequence)
    len = length(seq)
    coord_lists = Dict()
    for i in 1:len-KMER+1
        kmer = seq[i:i+KMER-1]
        key = kmer2key(kmer)
        if haskey(panel_kmer_coord, key)
            coord_lists[i] = panel_kmer_coord[key]
        end
    end
    return coord_lists
end

# clustering the kmer_coords
function cluster_ref(coord_lists)
    coords = Array{Coord, 1}()
    for (k, list) in coord_lists
        append!(coords, list)
    end
    return cluster(coords)
end