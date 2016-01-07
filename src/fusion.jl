# this standalone script can be put anywhere
#
# run from shell/command-line:
# julia fusion.jl -f <REF> -b <BED> -l <READ1> -r <READ2>
#
# required arguments:
#   -f, --ref REF      the reference folder, which contains chr1.fa,chr2fa..., 
#                      can be downloaded from
#                      http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
#   -b, --bed BED      a gene list with their coordination intervals
#   -l, --read1 READ1  read1 input fastq file, can be gzipped
#   -r, --read2 READ2  read2 input fastq file, can be gzipped

using FusionDirect
using ArgParse

function main(args)

    s = ArgParseSettings(description = "FusionDirect: detect gene fusion directly from fastq, without alignment required")

    @add_arg_table s begin
        "--ref",  "-f"
            help = "the reference folder, which contains chr1.fa, chr2fa..., can be downloaded from http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz"
            required = true
        "--bed",  "-b"
            help = "a gene list with their coordination intervals, see the example bed files in data folder"
            required = true
        "--read1",   "-l"
            help = "read1 input fastq file, can be gzipped"
            required = true
        "--read2",   "-r"
            help = "read2 input fastq file, can be gzipped"
            required = true
    end

    options = parse_args(s) # the result is a Dict{String,Any}
    detect(options["ref"], options["bed"], options["read1"],options["read2"])
end

# REPL is not supported
if !isinteractive()
    main(ARGS)
end