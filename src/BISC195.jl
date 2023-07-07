module BISC195

export normalizeDNA
export composition
export gc_content
export complement
export reverse_complement
export parse_fasta

# # uncomment the following line if you intend to use BioSequences types
#using BioSequences

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a String.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end


# Your code here.
# Don't forget to export your functions!

"""
    composition(sequence)

Counts the number of each type of base
in a DNA sequence and returns a dict of those counts
in the order A, C, G, T

Examples  
≡≡≡≡≡≡≡≡≡≡

julia> composition("ACCGGGTTTTN")
   Dict{Char,Int64} with 5 entries:
     'A' => 1
     'G' => 3
     'T' => 4
     'N' => 1
     'C' => 2
    
julia> composition("AAX")
   ERROR: Invalid base, X  
"""
function composition(sequence)
    sequence = normalizeDNA(sequence) # make uppercase string, check invalid bases
    a = c = g = t = n = 0 # sets all 4 variables to `0`

    for base in sequence
        ## add 1 to each base as it occurs
        if base == 'A'
            a = a + 1
        elseif base == 'C'
            c = c + 1
        elseif base == 'T'
            t = t + 1
        elseif base == 'G'
            g = g + 1
        elseif base == 'N'
            n = n + 1
        end
    end
    
    #create dictionary to return
    ans = Dict()
    ans['A'] = a
    ans['G'] = g
    ans['T'] = t
    ans['N'] = n
    ans['C'] = c
    
    return ans
end


"""
    gc_content(sequence)

Calculates the GC ratio of a DNA sequence.
The GC ratio is the total number of G and C bases divided by the total length of the sequence.

Examples  
≡≡≡≡≡≡≡≡≡≡

julia> gc_content("ATNG")
   0.25
   
julia> gc_content("ccccggggn")
   0.8888888888888888
"""
function gc_content(sequence)
    ans = composition(sequence)
    a = ans['A']
    g = ans['G']
    t = ans['T']
    n = ans['N']
    c = ans['C']
    seqlength = length(sequence)
    return (c+g) / seqlength
end


"""
    complement(base)

Get the DNA complement of the provided base:

    A <-> T
    G <-> C

Accepts uppercase or lowercase `String` or `Char`,
but always returns an uppercase `Char`.
If a valid base is not provided, the function throws an error.

Examples  
≡≡≡≡≡≡≡≡≡≡
julia> complement("ATTN")
"TAAN"

julia> complement("ATTAGC")
"TAATCG"
"""
function complement(sequence)
    complements = Dict('A' => "T",
                       'T' => "A",
                       'G' => "C",
                       'C' => "G",
                       'N' => "N")
    
    ans = ""
    sequence = uppercase(string(sequence))
    for base in sequence
        if (base ∉ "ATCGN")
            error("Invalid base $base")
        else
          ans *= complements[base]
        end
    end
    return ans
end


"""
    reverse_complement(sequence)

Takes a DNA sequence and returns the reverse complement
of that sequence.

Takes lowercase or uppercase sequences,
but always returns uppercase.

Examples
≡≡≡≡≡≡≡≡≡≡
    julia> reverse_complement("AAATTT")
    "AAATTT"

    julia> reverse_complement("GCAT")
    "ATGC"

    julia> reverse_complement("ATTAGC")
    "GCTAAT"

    julia> reverse_complement("ATN")
    "NAT"

    julia> rc = reverse_complement("TTGGG");

    julia> println(rc)
    CCCAA
"""
function reverse_complement(sequence)
    ## your code here
    sequence = string(uppercase(sequence))
    rev = reverse(sequence)
    ans = ""
    for base in rev
        ans*=complement(base)
    end
    return ans
end


"""
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing headers,
the other containing the entire sequence as a `String`.

Note: function does validate DNA sequences for correctness.

Example
≡≡≡≡≡≡≡≡≡
julia> ex1 = parse_fasta("data/ex1.fasta");

julia> ex1[1]
2-element Array{String,1}:
   "ex1.1 | easy"
   "ex1.2 | multiline"

julia> ex1[2]
2-element Array{String,1}:
   "AATTATAGC"
   "CGCCCCCCAGTCGGATT"

julia> ex2 = parse_fasta("data/ex2.fasta");
ERROR: invalid base H
"""
function parse_fasta(path)
    headers = []
    sequences = []
    sequence = ""
    for line in eachline(path)
        if length(line) >=1 && line[1] == '>'
            push!(headers, line[2:end])
            push!(sequences, sequence)
            sequence = ""
        else
            bases = ""
            for base in line
                if base ∉ "ATCGN"
                    error("invalid base $base")
                else
                    bases *= base
                end
            end
            sequence *= bases
        end
    end
    push!(sequences, sequence)
    
    if length(sequences) >= 1
        popfirst!(sequences)
    end
    
    return (headers, sequences)
end

end # module BISC195.jl
