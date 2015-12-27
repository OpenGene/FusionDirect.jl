
function detect_fusion(index::Index, seq::Sequence)
    len = length(seq)
    coords = [unknown_coord() for i in 1:len-KMER+1]
    for i in 1:len-KMER+1
        kmer = seq[i:i+KMER-1]
        key = kmer2key(kmer)
        if key in keys(index.data)
            coords[i] = index.data[key]
        end
    end

    # do segmentation based on coords
end