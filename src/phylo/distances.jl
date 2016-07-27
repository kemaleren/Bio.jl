# phylo/distances.jl
# ==================
#
# Types and methods for computing evolutionary and genetic distances.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Distance calculation infrastructure/api
# ---------------------------------------

# A type inheriting from EvoDist must implement the following methods.

# computes_variance(::Type{})
# return_type(::Type{})

abstract EvoDist


function count_differences(a::BioSequence, b::BioSequence, ::Type{EvoDist})
    @assert length(a) == length(b)
    return mismatches(a, b), length(a)
end

function count_differences_pairdel(a::BioSequence, b::BioSequence, ::Type{EvoDist})
    @assert length(a) == length(b)
    d = 0
    l = 0
    for i in eachindex(a)
        if isambiguous(a[i]) || isambiguous(b[i])
            continue
        else
            l += 1
            if a[i] != b[i]
                d += 1
            end
        end
    end
    return d, l
end


macro pairwise(expr)
    quote
        r = 1
        for i = 1:endof(seqs)
            for j = i+1:endof(seqs)
                res[r] = $(expr)(seqs[i], seqs[j], model)
                r += 1
            end
        end
    end
end

macro pairdel(expr)
    quote

    end
end







# Distance methods
# ----------------

# Indicate inline on these functions as they are called by other functions repeatedly in a loop.

## Raw distances.

immutable Raw <: EvoDist end

computes_variance(::Type{Raw}) = false
gamma_correction(::Type{Raw}) = false

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
            if a[i] != b[i]
                d += 1
            end
        end
    end
    return d / l
end

#function distance_pairdel(seqs::Vector{BioSequence}, res::Vector{Float64}, model::Raw)
#    @pairwise distance_pairdel
#end
#
#function distance(seqs::Vector{BioSequence}, res::Vector{Float64}, model::Raw)
#    @pairwise distance
#end

immutable JC69 <: EvoDist end

computes_variance(::Type{JC69}) = true
gamma_correction(::Type{JC69}) = true

@inline function distance(a::BioSequence, b::BioSequence, model::JC69)
    p = distance(a, b, Raw())
    return -0.75 * log(1 - 4 * p / 3)
end

@inline function distance(a::BioSequence, b::BioSequence, model::JC69, alpha::Float64)
    p = distance(a, b, Raw())
end

if (*gamma)\
  d[target] = 0.75 * *alpha*(pow(1 - 4*p/3, -1/ *alpha) - 1);\
else d[target] = ;\
if (*variance) {\
    if (*gamma) var[target] = p*(1 - p)/(pow(1 - 4*p/3, -2/(*alpha + 1)) * L);\
else var[target] = p*(1 - p)/(pow(1 - 4*p/3, 2)*L);\
}
