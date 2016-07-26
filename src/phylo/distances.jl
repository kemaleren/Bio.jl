# phylo/distances.jl
# ==================
#
# Types and methods for computing evolutionary distances.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

abstract MutationModel

immutable Raw <: MutationModel end

immutable JukesCantor69 <: MutationModel end

macro checkambig(expr)
    quote
        if isambiguous(a[i]) || isambiguous(b[i])
            break
        else
            l += 1
            $expr
        end
    end
end

macro countTsTv()

end

@inline function distance(a::BioSequence, b::BioSequence, model::Raw)
    @assert length(a) == length(b)
    return mismatches(a, b) / length(a)
end

@inline function distance_pairdel(a::BioSequence, b::BioSequence, model::Raw)
    @assert length(a) == length(b)
    d = 0
    l = 0
    for i in eachindex(a)
        @checkambig begin
            d += a[i] != b[i] ? 1 : 0
        end
    end
    return d / l
end
