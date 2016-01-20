const THRESHOLD = 30

# ref_path is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function detect(ref_path::AbstractString, bed_file::AbstractString, r1fq::AbstractString, r2fq::AbstractString="")
    index = index_bed(ref_path, bed_file)
    panel = index["panel"]
    panel_seq = index["seq"]
    panel_kmer_coord = index["kmer_coord"]
    ref_kmer_coords = index["ref_kmer_coords"]
    # pair end sequencing
    if r2fq != ""
        io = fastq_open_pair(r1fq, r2fq)
        fusion_pairs = Dict()
        while (pair = fastq_read_pair(io)) != false
            matches = detect_pair(panel_kmer_coord, pair)
            if length(matches)>1
                fusion_left, fusion_right = verify_fusion_pair(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, pair)
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
            name = name * panel[fusion_left.contig]["name"] * "|" * strand_name(fusion_left) * "|" * coord_to_chr(fusion_left, panel) * "_"
            name = name * panel[fusion_right.contig]["name"] * "|" * strand_name(fusion_right)  * "|" * coord_to_chr(fusion_right, panel) * "/"
            print(name, "1\n",pair.read1.sequence.seq,"\n")
            print(name, "2\n",pair.read2.sequence.seq,"\n")
        end
    end
end

function strand_name(coord)
    return coord.pos>0?"+":"-"
end

function coord_to_chr(coord, panel)
    from = panel[coord.contig]["from"]
    chr = panel[coord.contig]["chr"]
    pos = from + abs(coord.pos)
    return chr * ":" * string(pos)
end

function add_to_fusion_pair(fusion_pairs, fusion_left, fusion_right, pair)
    key = (min(fusion_left.contig, fusion_right.contig), max(fusion_left.contig, fusion_right.contig))
    if (key in keys(fusion_pairs)) == false
        fusion_pairs[key]=[]
    end
    push!(fusion_pairs[key], (fusion_left, fusion_right, pair))
end

function verify_fusion_se(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, sequence)
    seg1, coords1 = segment(ref_kmer_coords, panel_kmer_coord, sequence, panel_seq)
    if length(seg1) > 1
        return make_connected_fusion(panel_kmer_coord, panel_seq, seg1, coords1)
    end
end

function verify_fusion_pair(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, pair)
    offset, overlap_len, distance = overlap(pair)
    # this pair is overlapped, so merged it and segment the merged sequence
    if overlap_len>0 && distance<5
        merged_seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len)
        seg_result, coords = segment(ref_kmer_coords, panel_kmer_coord, merged_seq, panel_seq)
        if length(seg_result) > 1
            if !is_seq_connected_on_ref(seg_result, ref_kmer_coords, merged_seq)
                return make_connected_fusion(panel_kmer_coord, panel_seq, seg_result, coords)
            end
        end
    else
        seg_result1, coords1 = segment(ref_kmer_coords, panel_kmer_coord, pair.read1.sequence, panel_seq)
        if length(seg_result1) > 1
            if !is_seq_connected_on_ref(seg_result1, ref_kmer_coords, pair.read1.sequence)
                return make_connected_fusion(panel_kmer_coord, panel_seq, seg_result1, coords1)
            end
        end
        seg_result2, coords2 = segment(ref_kmer_coords, panel_kmer_coord, pair.read2.sequence, panel_seq)
        if length(seg_result2) > 1
            if !is_seq_connected_on_ref(seg_result2, ref_kmer_coords, pair.read2.sequence)
                return make_connected_fusion(panel_kmer_coord, panel_seq, seg_result2, coords2)
            end
        end
        if length(seg_result1)==0 || length(seg_result2) == 0 || seg_result1[1][1] == 0 || seg_result2[1][1] == 0
            return false, false
        end
        # fusion of different contigs on a pair
        if coords1[seg_result1[1][1]].contig != coords2[seg_result2[1][1]].contig
            if is_pair_connected_on_ref(pair, ref_kmer_coords)
                return false, false
            end

            l1 = seg_result1[1][1]
            l2 = seg_result2[1][1]
            len1 = length(coords1)
            len2= length(coords2)
            fusion_left = Coord(coords1[l1].contig, coords1[l1].pos + (len1 - l1))
            # reverse the direction of fusion right to align the pair in same direction
            fusion_right = Coord(coords2[l2].contig, coords2[l2].pos + (len2 - l2), -1)

            return fusion_left, fusion_right
        end
    end
    return false, false
end

function is_seq_connected_on_ref(seg_result, ref_kmer_coords, seq)
    coords_list = stat_ref(ref_kmer_coords, seq)
    clusters_on_ref = clustering_ref(coords_list)
    for cluster in clusters_on_ref
        (left, right) = span_ref(cluster, length(seq))
        # nearly cover whold sequence
        if (right - left) > length(seq) - 30
            return true
        end
        # nearly covers the two regions
        if left < seg_result[1][1] + 10 && right > seg_result[2][2] - 10
            return true
        end
    end
    return false
end

function is_pair_connected_on_ref(pair, ref_kmer_coords)
    coords_list1 = stat_ref(ref_kmer_coords, pair.read1.sequence)
    clusters1 = clustering_ref(coords_list1)
    coords_list2 = stat_ref(ref_kmer_coords, pair.read2.sequence)
    clusters2 = clustering_ref(coords_list2)

    # try to find if there is any intersection between clusters_on_ref1 and clusters_on_ref2
    for c1 in clusters1
        (left, right) = span_ref(c1, length(seq))
        # for any cluster>30 in the read1
        if (right - left) > 30
            for c2 in clusters2
                # if we find a cluster from read2, which is >30 and near c1
                # it means they are connected
                (left2, right2) = span_ref(c1, length(seq))
                if (right2 - left2) > 30 && min_distance(c1, c2) < 500
                    return true
                end
            end
        end
    end

    return false
end

function make_connected_fusion(panel_kmer_coord, panel_seq, seg_result, coords)
    l1 = seg_result[1][1]
    r1 = seg_result[1][2]
    l2 = seg_result[2][1]
    r2 = seg_result[2][2]
    conjunct = div(r1+l2, 2)
    #display_coords(coords)
    fusion_left = Coord(coords[l1].contig, coords[l1].pos + (conjunct - l1))
    fusion_right = Coord(coords[l2].contig, coords[l2].pos + (conjunct - l2))
    return fusion_left, fusion_right
end

# detect fusion and find the segmentations
function segment(ref_kmer_coords, panel_kmer_coord::KmerCoord, seq::Sequence, panel_seq)
    counts, coords = stat(panel_kmer_coord, seq)
    clusters = clustering(coords)
    regions = Dict()
    for cluster in clusters
        # skip small clusters
        if length(cluster) < THRESHOLD
            continue
        end
        left, right = span(cluster, length(coords))
        regions[left]=right
    end
    # filter out those regions which are mostly inside other regions
    seg_result = filter_region(regions, coords)
    return seg_result, coords
end

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
function span(cluster, seqlen)
    left = seqlen
    right = 1
    for point in cluster
        left = min(left, point)
        right = max(right, point)
    end
    right = min(right+KMER-1, seqlen)
    return left, right
end

# span the cluster on the reference to find the region
function span_ref(cluster, seqlen)
    left = seqlen
    right = 1
    for point in cluster
        left = min(left, point.first)
        right = max(right, point.first)
    end
    right = min(right+KMER-1, seqlen)
    return left, right
end

# clustering by distance
function clustering(coords::Array{Coord, 1})
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
    if abs(p2-p1)<40 && abs((p2-p1) - dis * sign(coords[p1].pos)) < 4
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
        if haskey(ref_kmer_coords, key)
            coord_lists[i] = ref_kmer_coords[key]
        end
    end
    return coord_lists
end

# clustering of the kmer_coords
function clustering_ref(coord_lists)
    total = Set()
    # covert the coord list to a set of pos=>coord
    for (pos, arr) in coord_lists
        for coord in arr
            push!(total, pos=>coord)
        end
    end

    clusters = []
    while !isempty(total)
        # create a new cluster
        seed = pop!(total)
        # create a new cluster
        cluster = Set()
        push!(cluster, seed)
        # create a new working set
        working = Set()
        push!(working, seed)
        while !isempty(working)
            current = pop!(working)
            for p in total
                # if it is not valid, just pop it out
                if consistent_on_ref(current, p)
                    pop!(total, p)
                    push!(working, p)
                end
            end
            # current is done, push the pos to current cluster
            push!(cluster, current)
        end
        push!(clusters, cluster)
    end
    return clusters
end

function consistent_on_ref(p1, p2)
    dis = p2.second - p1.second
    if dis > 1000
        return false
    end

    # the line of these two points be near 45 degree
    if abs(p2.first-p1.first)<40 && abs((p2.first-p1.first) - dis * sign(p1.second.pos)) < 4
        return true
    end

    return false
end