const THRESHOLD = 30

# ref_folder is a folder contains fasta files by chromosomes
# like chr1.fa, chr2.fa ...
function detect(ref_folder::AbstractString, bed_file::AbstractString, r1fq::AbstractString, r2fq::AbstractString)
    index, bed = index_bed(ref_folder, bed_file)
    if r2fq != ""
        io = fastq_open_pair(r1fq, r2fq)
        i = 0
        while (pair = fastq_read_pair(io)) != false
            i += 1
            matches = detect_pair(index, pair)
            if length(matches)>1
                print("\n")
                println(pair.read1.sequence)
                println(pair.read2.sequence)
                for m in matches
                    print(bed[m] * " ")
                end
            elseif i%10000 == 0
                print("$i ")
            end
        end
    end
end

# simple fusion detection
function detect_read(index::Index, seq::Sequence)
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

    matches = Array{Int16, 1}()
    for (k,v) in counts
        if v > THRESHOLD
            push!(matches, k)
        end
    end
    return matches
end

function detect_pair(index::Index, pair::FastqPair)
    match1 = detect_read(index, pair.read1.sequence)
    match2 = detect_read(index, pair.read2.sequence)
    for m in match2
        if (m in match1) == false
            push!(match1, m)
        end
    end
    return match1
end