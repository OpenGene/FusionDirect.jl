include("white_list.jl")

const THRESHOLD = 30
const MIN_READ_SUPPORT = 2
const MAX_CLUSTER_DISTANCE = 100
const MATCH_EDIT_DISTANCE_LIMIT = 5

const FUSION_ON_MERGED_READ = 0
const FUSION_ON_READ1 = 1
const FUSION_ON_READ2 = 2
const FUSION_ON_CROSS_READS = 3

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
                fusion_left, fusion_right, fusion_site, seg_result = verify_fusion_pair(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, pair)
                if fusion_left!=false && fusion_right!=false && distance(fusion_left, fusion_right)>1000
                    # drop the low quality reads
                    if pass_quality_filter(pair.read1.quality) && pass_quality_filter(pair.read2.quality)
                        add_to_fusion_pair(fusion_pairs, fusion_left, fusion_right, fusion_site, seg_result, pair)
                    end
                end
            end
        end
        printed = print_fusion_pair(fusion_pairs, panel_seq, panel)
        if printed == 0
            println("# no fusion detected")
        end
    else
        error("FusionDirect only supports pair-end sequencing data")
    end
end

export pass_quality_filter
function pass_quality_filter(quality)
    # '0' = q15 as threshold
    const QUAL_T = '0'
    lowqual = 0
    for q in quality.qual
        if q<QUAL_T
            lowqual += 1
        end
    end
    # allow up to 8% low quality bases
    if lowqual/length(quality) > 0.08
        return false
    else
        return true
    end
end

function print_fusion_pair(fusion_pairs, panel_seq, panel)
    printed = 0
    gencode = gencode_load("GRCh37")
    for (fusion_key, all_fusion_reads) in fusion_pairs
        contig1, contig2 = fusion_key
        name1 = panel[contig1]["name"]
        name2 = panel[contig2]["name"]

        # clustering of fusion pairs
        clustered_fusion_reads = fusion_clustering(all_fusion_reads, panel)

        for fusion_reads in clustered_fusion_reads
            unique_fusion_reads, read_support = get_unique_fusion_pairs(fusion_reads)
            total_num = length(fusion_reads)
            unique_num = length(unique_fusion_reads)

            # filter out those fusion with very few reads support (like only 1 unique reads) or have no consistent fusion readsd
            if unique_num< MIN_READ_SUPPORT || !have_consistency(unique_fusion_reads)
                # check if this fusion is in the white list
                if !is_important_fusion(unique_fusion_reads, panel, gencode)
                    continue
                end
            end

            # give the fusion as a fasta comment line
            println("#Fusion:$name1-$name2 (total: $total_num, unique: $unique_num)")
            # display all reads support this fusion
            for read in unique_fusion_reads
                fusion_left, fusion_right, fusion_site, conjunct, pair = read
                support = read_support[read]
                name = ">$support\_"
                name = name * get_fusion_site_string(conjunct, fusion_site, pair) * "_"
                name = name * get_gene_coord_string(fusion_left, panel, gencode) * "_"
                name = name * get_gene_coord_string(fusion_right, panel, gencode) * "/"
                print(name, "1\n",pair.read1.sequence.seq,"\n")
                # println("# ", pair.read1.quality.qual)
                print(name, "2\n",pair.read2.sequence.seq,"\n")
                # println("# ", pair.read2.quality.qual)

                if fusion_site == FUSION_ON_MERGED_READ
                    offset, overlap_len, distance = overlap(pair)
                    merged_seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len, offset)
                    print(name, "merged\n",merged_seq.seq,"\n")
                end

                printed += 1
            end
        end
    end
    return printed
end

function is_important_fusion(fusion_reads, panel, gencode)
    if length(IMPORTANT_FUSIONS) == 0
        return false
    end
    for read in fusion_reads
        fusion_left, fusion_right, fusion_site, conjunct, pair = read

        find1 = get_coord_gene(panel, fusion_left, gencode)
        find2 = get_coord_gene(panel, fusion_right, gencode)

        if find1[1] == false || find2[1] == false
            continue
        end

        for (fusion1, fusion2) in IMPORTANT_FUSIONS
            if is_fusion_identical(find1, fusion1) && is_fusion_identical(find2, fusion2)
                return true
            elseif is_fusion_identical(find1, fusion2) && is_fusion_identical(find2, fusion1)
                return true
            end
        end
    end
    return false
end

function is_fusion_identical(find, fusion)
    if find[1] != fusion[1]
        return false
    elseif fusion[2] != "*" && find[2] != fusion[2]
        return false
    elseif fusion[3] != "*" && find[3] != parse(Int, fusion[3])
        return false
    else
        return true
    end
end

function get_coord_gene(panel, coord, gencode)
    genename = panel[coord.contig]["name"]
    from = panel[coord.contig]["from"]
    chr = panel[coord.contig]["chr"]
    pos = from + abs(coord.pos)

    str = genename
    genelocs = gencode_locate(gencode, chr, pos)
    if length(genelocs) > 0
        for loc in genelocs
            if genename == loc["gene"]
                kind = loc["type"]
                number = loc["number"]
                return genename, kind, number
            end
        end
    end
    return false, false, false
end

function have_consistency(unique_fusion_reads)
    for i in 1:length(unique_fusion_reads) - 1
        for j in i+1:length(unique_fusion_reads)
            fusion_left1, fusion_right1, fusion_site1, conjunct1, pair1 = unique_fusion_reads[i]
            fusion_left2, fusion_right2, fusion_site2, conjunct2, pair2 = unique_fusion_reads[j]
            if fusion_site1 == FUSION_ON_CROSS_READS || fusion_site2 == FUSION_ON_CROSS_READS
                if is_bad_pair(pair1, pair2)
                    continue
                else
                    return true
                end
            else
                r11, r12 = get_fusion_seqs(unique_fusion_reads[i])
                r21, r22 = get_fusion_seqs(unique_fusion_reads[j])
                if is_consistent(r11, r21) && is_consistent(r12, r22)
                    return true
                elseif is_consistent(r11, r22) && is_consistent(r12, r21)
                    return true
                end
            end
        end
    end
    return false
end

function is_bad_pair(pair1, pair2)
    # filter the case that one read is found as dup, but the other read has lots of wrong bases
    dup, ed = is_dup(pair1.read1.sequence.seq, pair2.read1.sequence.seq)
    if dup
        if hamming(pair1.read2.sequence.seq[1:20], pair2.read2.sequence.seq[1:20]) <= 5
            return true
        end
    end
    dup, ed = is_dup(pair1.read1.sequence.seq, pair2.read2.sequence.seq)
    if dup
        if hamming(pair1.read2.sequence.seq[1:20], pair2.read1.sequence.seq[1:20]) <= 5
            return true
        end
    end
    return false
end

function is_consistent(r1, r2)
    const CONSISTENT_T = 5
    overlap_len = min(length(r1), length(r2))
    overlap_r1 = r1[1:overlap_len]
    overlap_r2 = r2[1:overlap_len]
    ed = edit_distance(overlap_r1.seq, overlap_r2.seq)
    return ed <= CONSISTENT_T
end

function get_fusion_seqs(fusion_read)
    fusion_left, fusion_right, fusion_site, conjunct, pair = fusion_read
    seq = nothing
    if fusion_site == FUSION_ON_READ1
        seq = pair.read1.sequence
    elseif fusion_site == FUSION_ON_READ2
        seq = pair.read2.sequence
    else
        offset, overlap_len, distance = overlap(pair)
        seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len, offset)
    end
    len = length(seq)
    return ~(seq[1:conjunct]), seq[conjunct:len]
end

function get_gene_coord_string(coord, panel, gencode)
    genename = panel[coord.contig]["name"]
    from = panel[coord.contig]["from"]
    chr = panel[coord.contig]["chr"]
    pos = from + abs(coord.pos)

    str = genename
    genelocs = gencode_locate(gencode, chr, pos)
    if length(genelocs) > 0
        for loc in genelocs
            if genename == loc["gene"]
                kind = loc["type"]
                number = loc["number"]
                str *= ":$kind"
                str *= ":$number"
                break
            end
        end
    end
    chrloc = chr * ":" * string(pos)
    str *= "|" * strand_name(coord) * chrloc
end

function fusion_clustering(all_fusion_reads, panel)
    const FUSION_THRESHOLD = 100
    fusion_clusters = []
    for r in all_fusion_reads
        found = false
        for cluster in fusion_clusters
            for read in cluster
                left = read[1]
                right = read[2]
                left1 = r[1]
                right1 = r[2]
                if (abs_distance(left,left1)< MAX_CLUSTER_DISTANCE && abs_distance(right, right1)<MAX_CLUSTER_DISTANCE) || (abs_distance(left, right1)< MAX_CLUSTER_DISTANCE && abs_distance(left1, right)<MAX_CLUSTER_DISTANCE)
                    found = true
                    push!(cluster, r)
                    break
                end
            end
            if found == true
                break
            end
        end
        if found == false
            new_cluster = [r]
            push!(fusion_clusters, new_cluster)
        end
    end
    return fusion_clusters
end

function get_fusion_site_string(conjunct, fusion_site, pair)
    if fusion_site == FUSION_ON_MERGED_READ
        return "merged_$conjunct"
    elseif fusion_site == FUSION_ON_READ1
        return "read1_$conjunct"
    elseif fusion_site == FUSION_ON_READ2
        return "read2_$conjunct"
    else
        return "crosspair_0"
    end
end

function get_unique_fusion_pairs(fusion_reads)
    unique_fusion_reads=[]
    read_support = Dict()
    # display all reads support this fusion
    for read in fusion_reads
        fusion_left, fusion_right, fusion_site, conjunct, pair = read
        unique = true
        for uread in unique_fusion_reads
            uleft, uright, usite, uconjunct, upair = uread
            if is_dup_pair(upair, pair)
                unique = false
                read_support[uread] += 1
                break
            end
        end
        if unique
            push!(unique_fusion_reads, read)
            read_support[read] = 1
        end
    end
    return unique_fusion_reads, read_support
end

function is_dup_pair(pair1::FastqPair, pair2::FastqPair)
    # edit distance threshold
    const ED_T = length(pair1.read1) * 0.1

    seq11 = pair1.read1.sequence.seq
    seq12 = pair1.read2.sequence.seq
    seq21 = pair2.read1.sequence.seq
    seq22 = pair2.read2.sequence.seq

    dup1, ed1 = is_dup(seq11, seq21)
    dup2, ed2 = is_dup(seq12, seq22)
    if dup1 == true && (dup2 == true || ed2 <= ED_T || share_start_bases(seq12, seq22, 5))
        return true
    elseif dup2 == true && (dup1 == true || ed1 <= ED_T || share_start_bases(seq11, seq21, 5))
        return true
    end
    dup1, ed1 = is_dup(seq11, seq22)
    dup2, ed2 = is_dup(seq12, seq21)
    if dup1 == true && (dup2 == true || ed2 <= ED_T || share_start_bases(seq12, seq21, 5))
        return true
    elseif dup2 == true && (dup1 == true || ed1 <= ED_T || share_start_bases(seq11, seq22, 5))
        return true
    end

    # a work around for bad sequence quality
    if share_start_bases(seq11, seq21, 10) && share_start_bases(seq11, seq21, 10)
        return true
    elseif share_start_bases(seq12, seq21, 10) && share_start_bases(seq11, seq22, 10)
        return true
    end

    return false
end

function share_start_bases(s1::ASCIIString, s2::ASCIIString, len = 5)
    return hamming(s1[1:len], s2[1:len]) <= 1
end

function is_dup(s1::ASCIIString, s2::ASCIIString)
    ed = edit_distance(s1, s2)
    hm = hamming(s1, s2)
    if ed <= 3
        return true, ed
    elseif hm <= length(s1) * 0.1
        return true, ed
    elseif edit_distance(s1, s2) <=5 && ( hamming(s1[1:5], s2[1:5])<=1 || hamming(s1[length(s1)-4:length(s1)], s2[length(s2)-4:length(s2)])<=1 )
        return true, ed
    else
        return false, ed
    end
end

# hamming distance
function hamming(s1::ASCIIString, s2::ASCIIString)
    d = 0
    len = min(length(s1), length(s2))
    for i in 1:len
        if s1[i] != s2[i]
            d += 1
        end
    end
    return d
end

function strand_name(coord)
    return coord.pos>0?"+":"-"
end

function add_to_fusion_pair(fusion_pairs, fusion_left, fusion_right, fusion_site, seg_result, pair)
    key = (min(fusion_left.contig, fusion_right.contig), max(fusion_left.contig, fusion_right.contig))
    if (key in keys(fusion_pairs)) == false
        fusion_pairs[key]=[]
    end
    push!(fusion_pairs[key], (fusion_left, fusion_right, fusion_site, seg_result, pair))
end

function verify_fusion_se(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, sequence)
    seg1, coords1 = segment(ref_kmer_coords, panel_kmer_coord, sequence, panel_seq)
    if length(seg1) > 1
        return make_connected_fusion(seg1, coords1)
    end
end

function found_seq_near_coord(panel_seq, seq, coord, pos_in_seq, threshold = 10)
    search_pos = abs(coord.pos -  pos_in_seq)
    seqlen = length(seq)
    # in reverse strand
    if coord.pos < 0
        seq = ~seq
        search_pos -= seqlen
    end
    const SEARCH_WINDOW = 10

    geneseq = panel_seq[coord.contig]
    genelen = length(geneseq)

    for start = max(search_pos-SEARCH_WINDOW, 1) : min(search_pos+SEARCH_WINDOW, genelen-seqlen-SEARCH_WINDOW)
        genepart = uppercase(geneseq.seq[start : start+seqlen-1])
        ed = edit_distance(genepart, seq.seq)
        if ed < threshold
            return true
        end
    end
    return false
end

# check if the splitted two sequences are consistent with ref
function splitted_consistent_with_ref(panel_seq, seq, fusion_left, fusion_right, left, right, conjunct)
    const CONJUNCT_WINDOW = 5
    seqlen = length(seq)

    #println("splitted_consistent_with_ref", " ", conjunct)
    #println(seq)
    for conj = -CONJUNCT_WINDOW:CONJUNCT_WINDOW
        break_point = conjunct + conj
        left_seq = seq[left:break_point]
        right_seq = seq[break_point+1:right]
        # coord of the starts of these two seqs
        left_coord = Coord(fusion_left.contig, fusion_left.pos - (conjunct-left))
        right_coord = Coord(fusion_right.contig, fusion_right.pos + conj)
        if found_seq_near_coord(panel_seq, left_seq, left_coord, 1, MATCH_EDIT_DISTANCE_LIMIT) && found_seq_near_coord(panel_seq, right_seq, right_coord, 1, MATCH_EDIT_DISTANCE_LIMIT)
            return true
        end
    end
    return false

end

function verify_fusion_pair(ref_kmer_coords, panel_kmer_coord, panel, panel_seq, pair)
    offset, overlap_len, distance = overlap(pair)
    # this pair is overlapped, so merged it and segment the merged sequence
    if overlap_len>0 && distance<5
        merged_seq = simple_merge(pair.read1.sequence, pair.read2.sequence, overlap_len, offset)
        seg_result, coords = segment(ref_kmer_coords, panel_kmer_coord, merged_seq, panel_seq)
        if length(seg_result) > 1
            if !is_seq_connected_on_ref(seg_result, ref_kmer_coords, merged_seq)
                fusion_left, fusion_right, conjunct = make_connected_fusion(seg_result, coords)
                l = seg_result[1][1]
                r = seg_result[2][2]
                if !found_seq_near_coord(panel_seq, merged_seq[l:r], fusion_left, conjunct-l+1) && !found_seq_near_coord(panel_seq, merged_seq[l:r], fusion_right, conjunct-l+1)
                    if splitted_consistent_with_ref(panel_seq, merged_seq, fusion_left, fusion_right, l, r, conjunct)
                        return (fusion_left, fusion_right, FUSION_ON_MERGED_READ, conjunct)
                    end
                end
            end
        end
    else
        seg_result1, coords1 = segment(ref_kmer_coords, panel_kmer_coord, pair.read1.sequence, panel_seq)
        if length(seg_result1) > 1
            if !is_seq_connected_on_ref(seg_result1, ref_kmer_coords, pair.read1.sequence)
                fusion_left, fusion_right, conjunct = make_connected_fusion(seg_result1, coords1)
                l = seg_result1[1][1]
                r = seg_result1[2][2]
                if !found_seq_near_coord(panel_seq, pair.read1.sequence[l:r], fusion_left, conjunct-l+1) && !found_seq_near_coord(panel_seq, pair.read1.sequence, fusion_right, conjunct-l+1)
                    if splitted_consistent_with_ref(panel_seq, pair.read1.sequence, fusion_left, fusion_right, l, r, conjunct)
                        return (fusion_left, fusion_right, FUSION_ON_READ1, conjunct)
                    end
                end
            end
        end
        seg_result2, coords2 = segment(ref_kmer_coords, panel_kmer_coord, pair.read2.sequence, panel_seq)
        if length(seg_result2) > 1
            if !is_seq_connected_on_ref(seg_result2, ref_kmer_coords, pair.read2.sequence)
                fusion_left, fusion_right, conjunct = make_connected_fusion(seg_result2, coords2)
                l = seg_result2[1][1]
                r = seg_result2[2][2]
                if !found_seq_near_coord(panel_seq, pair.read2.sequence[l:r], fusion_left, conjunct-l+1) && !found_seq_near_coord(panel_seq, pair.read2.sequence, fusion_right, conjunct-l+1)
                    if splitted_consistent_with_ref(panel_seq, pair.read2.sequence, fusion_left, fusion_right, l, r, conjunct)
                        return (fusion_left, fusion_right, FUSION_ON_READ2, conjunct)
                    end
                end
            end
        end
        if length(seg_result1)==0 || length(seg_result2) == 0 || seg_result1[1][1] == 0 || seg_result2[1][1] == 0
            return false, false, false, 0
        end
        # fusion of different contigs on a pair
        if coords1[seg_result1[1][1]].contig != coords2[seg_result2[1][1]].contig
            if is_pair_connected_on_ref(pair, ref_kmer_coords)
                return false, false, false, 0
            end

            l1 = seg_result1[1][1]
            l2 = seg_result2[1][1]
            len1 = length(coords1)
            len2= length(coords2)
            fusion_left = Coord(coords1[l1].contig, coords1[l1].pos + (len1 - l1))
            # reverse the direction of fusion right to align the pair in same direction
            fusion_right = Coord(coords2[l2].contig, coords2[l2].pos + (len2 - l2), -1)

            # check if consistent with ref
            left_coord = Coord(coords1[l1].contig, coords1[l1].pos - l1)
            right_coord = Coord(coords2[l2].contig, coords2[l2].pos - l2)
            if found_seq_near_coord(panel_seq, pair.read1.sequence, left_coord, 1, MATCH_EDIT_DISTANCE_LIMIT) && found_seq_near_coord(panel_seq, pair.read2.sequence, right_coord, 1, MATCH_EDIT_DISTANCE_LIMIT)
                return (fusion_left, fusion_right, FUSION_ON_CROSS_READS, 0)
            end
        end
    end
    return false, false, false, 0
end

function is_seq_connected_on_ref(seg_result, ref_kmer_coords, original_seq)
    l = seg_result[1][1]
    r = seg_result[2][2]
    seq = original_seq[l:r]
    coords_list = stat_ref(ref_kmer_coords, seq)
    clusters_on_ref = clustering_ref(coords_list)
    if length(clusters_on_ref) == 0
        # the cluster set is too big
        # so it is considered in a heavy repeating region
        # treat it as not fusion directly
        # println("seq not clustered")
        return true
    end
    for cluster in clusters_on_ref
        (left, right) = span_ref(cluster, length(seq))
        # nearly cover whole sequence
        if (right - left) > length(seq) - 20
            # println("whole seq on ref")
            return true
        end
        # nearly covers the two regions
        if left <  15 && right > length(seq) - 15
            # println("two regions on ref")
            return true
        end
    end
    # println("seq not on ref")
    return false
end

function is_pair_connected_on_ref(pair, ref_kmer_coords)
    coords_list1 = stat_ref(ref_kmer_coords, pair.read1.sequence)
    clusters1 = clustering_ref(coords_list1)
    coords_list2 = stat_ref(ref_kmer_coords, pair.read2.sequence)
    clusters2 = clustering_ref(coords_list2)
    # TODO
    # one or two reads are not clustered because the cluster set is too big
    # so they are considered in a heavy repeating region
    # treat it as not fusion directly
    if length(clusters1) == 0 || length(clusters2) == 0
        # println("pair not clustered")
        return true
    end

    # try to find if there is any intersection between clusters_on_ref1 and clusters_on_ref2
    for c1 in clusters1
        (left, right) = span_ref(c1, length(pair.read1.sequence))
        # for any cluster>30 in the read1
        if (right - left) > 30
            for c2 in clusters2
                # if we find a cluster from read2, which is >30 and near c1
                # it means they are connected
                (left2, right2) = span_ref(c2, length(pair.read2.sequence))
                if (right2 - left2) > 30 && min_distance(c1, c2) < 500
                    # two reads are close
                    # println("\npair too close, distance:", min_distance(c1, c2))
                    return true
                end
            end
        end
    end

    return false
end

function min_distance(c1, c2)
    dis = Inf
    for (pos1, coord1) in c1
        for (pos2, coord2) in c2
            # the pair should be in different strand
            d = abs(distance(coord1, Coord(coord2.contig, -coord2.pos)))
            if d < dis
                dis = d
            end
        end
    end
    return dis
end

function make_connected_fusion(seg_result, coords)
    l1 = seg_result[1][1]
    r1 = seg_result[1][2]
    l2 = seg_result[2][1]
    r2 = seg_result[2][2]
    conjunct = div(r1+l2, 2)
    #display_coords(coords)
    fusion_left = Coord(coords[l1].contig, coords[l1].pos + (conjunct - l1))
    fusion_right = Coord(coords[l2].contig, coords[l2].pos + (conjunct - l2))
    return fusion_left, fusion_right, conjunct
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
        # if two many entries on this kmer
        # we consider this kmer_coords as less informative
        # for performance considering, we skip it
        if length(arr) < 300
            for coord in arr
                push!(total, pos=>coord)
            end
        end
    end

    clusters = []
    # println(length(total))
    # TODO: handle big cluster set
    # the cluster set is too big
    # so it is considered in a heavy repeating region
    # treat it as not fusion directly for performance considering
    if(length(total) > 10000)
        return clusters
    end
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