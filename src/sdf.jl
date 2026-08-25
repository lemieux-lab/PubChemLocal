# Low-level SDF parsing shared by the Compound and Substance ingestion paths.
#
# PubChem bulk dumps (Compound and Substance) are both plain multi-record SDF:
# each record is a CTAB mol block followed by zero or more
#   > <TAG>
#   value line(s)
# blocks, blank-line separated, with records themselves separated by a line
# containing only "$$$$".

using CodecZlib

"""
    read_gz_text(fn) -> String

Decompress a `.gz` file fully into memory as text.
"""
function read_gz_text(fn::AbstractString)
    open(fn, "r") do fio
        io = GzipDecompressorStream(fio)
        read(io, String)
    end
end

"""
    split_records(text::AbstractString) -> Vector{SubString}

Split raw (decompressed) SDF text on the "\$\$\$\$" record delimiter, dropping
empty fragments (trailing newline, stray blank records).
"""
split_records(text::AbstractString) = filter(!isempty, split(text, "\$\$\$\$\n"))

const TAG_BLOCK_RE = r"> +<(.+?)>\s*\n(.*)"s  # 's' = dotall, so '.' spans lines

# Not anchored at the block start: the mol block's "M  END" and the *first*
# property tag share a blank-line-delimited chunk (blank lines only appear
# *between* tag blocks, not before the first one), so the first chunk after
# split(record, "\n\n") looks like "...M  END\n> <TAG>\nvalue" — the tag
# pattern has to be found within the chunk, not required to start it.

"""
    parse_tags(record::AbstractString) -> Dict{String, Vector{String}}

Parse the `> <TAG>\\nvalue...` property blocks of one SDF record into a
`Dict` mapping tag name to *all* of its values.

A tag's value can itself span several lines (PubChem uses this for
multi-valued fields such as `PUBCHEM_SUBSTANCE_SYNONYM`, one synonym per
line) and, rarely, a tag block can repeat within a record; both cases are
collected rather than silently dropped to the first line/first occurrence.

Note: an earlier version of this parser (ported from the original
`PubChem.jl` script) only kept the first line of each block, which is
immaterial for the largely single-valued Compound tags but would have
silently discarded most synonyms once applied to Substance records.
"""
function parse_tags(record::AbstractString)
    tags = Dict{String,Vector{String}}()
    for block in split(record, "\n\n")
        m = match(TAG_BLOCK_RE, block)
        m === nothing && continue
        tag = m[1]
        values = filter(!isempty, strip.(split(rstrip(m[2]), '\n')))
        isempty(values) && continue
        append!(get!(() -> String[], tags, tag), values)
    end
    return tags
end

"""
    chunk(v::AbstractVector, blk_size::Integer) -> Vector{<:AbstractVector}

Split `v` into contiguous chunks of at most `blk_size` elements, for
parallel processing over `Threads.@spawn`/`Threads.foreach`.
"""
function chunk(v::AbstractVector, blk_size::Integer)
    n = length(v)
    nblk = cld(n, blk_size)
    return [v[((b-1)*blk_size+1):min(n, b*blk_size)] for b in 1:nblk]
end
