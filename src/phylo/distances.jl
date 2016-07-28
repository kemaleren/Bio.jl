# phylo/distances.jl
# ==================
#
# Types and methods for computing evolutionary and genetic distances.
#
# Part of the Bio.Phylo module.
#
# This file is a part of BioJulia. License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md


# Types
# -----

abstract EvoDist
abstract UncorrectedDist <: EvoDist
abstract CorrectedDist <: EvoDist
abstract TsTv <: CorrectedDist

immutable Raw <: UncorrectedDist end
immutable JC69 <: CorrectedDist end
immutable K80 <: TsTv end
immutable F81 end
immutable K81 end
immutable F84 end
immutable T92 end
immutable TN93 end

# Macros
# ------
macro pairwiseDel(expr)
    :(if !isambiguous(a[i]) && !isambiguous(b[i])
        nSites += 1
        $expr
    end)
end

# Methods
# -------

# Count the different types of difference.

function count_differences(a::BioSequence, b::BioSequence, ::Type{EvoDist})
    @assert length(a) == length(b)
    return mismatches(a, b), length(a)
end

function count_differences_pairdel(a::BioSequence, b::BioSequence, ::Type{EvoDist})
    @assert length(a) == length(b)
    nDifferences = 0
    nSites = 0
    for i in eachindex(a)
        @pairwiseDel if a[i] != b[i]
            nDifferences += 1
        end
    end
    return nDifferences, nSites
end

function countTsTv(a::BioSequence, b::BioSequence, model::TsTv)
    @assert length(a) == length(b)
    nDifferences = 0
    nTransitions = 0
    nSites = length(a)
    for i in eachindex(a)
        if a[i] != b[i]
            nd += 1
            if (isPurine(a[i]) && isPurine(b[i])) || isPyrimidine(a[i]) && isPyrimidine(a[i])
                nTransitions += 1
            end
        end
    end
    nTransversions = nDifferences - nTransitions
    return nDifferences, nTransitions, nTransversions, nSites
end

function countTsTv_pairdel(a::BioSequence, b::BioSequence, model::TsTv)
    @assert length(a) == length(b)
    nDifferences = 0
    nTransitions = 0
    nSites = 0
    for i in eachindex(a)
        @pairwiseDel if a[i] != b[i]
            nDifferences += 1
            if (isPurine(a[i]) && isPurine(b[i])) || isPyrimidine(a[i]) && isPyrimidine(a[i])
                nTransitions += 1
            end
        end
    end
    nTransversions = nDifferences - nTransitions
    return nDifferences, nTransitions, nTransversions, nSites
end







# Distance methods
# ----------------

# Indicate inline on these functions as they are called by other functions repeatedly in a loop.

## Raw distances.



@inline function distance(a::BioSequence, b::BioSequence, model::Raw)
    d, l = count_differences(a, b, model)
    return d / l
end

@inline function distance_pairdel(a::BioSequence, b::BioSequence, model::Raw)
    d, l = count_differences_pairdel(a, b, model)
    return d / l
end


# JC69 Distance computation.

function model_correction(x::Float64, model::JC69)
    return -0.75 * log(1 - 4 * x / 3)
end

function model_correction(x::Float64, gamma::Float64, model::JC69)
    return 0.75 * alpha * ( (1 - 4 * p / 3) ^ (-1 / alpha) - 1)
end

function variance(x::Float64, L::Int, model::JC69)
    return x * (1 - x) / (((1 - 4 * p / 3) ^ 2) * L)
end

function variance(x::Float64, L::Int, gamma::Float64, model::JC69)
    return x * (1 - x)/(((1 - 4 * x / 3) ^ (-2 / (alpha + 1))) * L)
end

function distance(a::BioSequence, b::BioSequence, model::JC69)
    d, l = count_differences(a, b, model)
    p = d / l
    D = model_correction(p, model)
    V = variance(p, l, model)
    return D, V
end

function distance(a::BioSequence, b::BioSequence, model::JC69, gamma::Float64)
    d, l = count_differences(a, b, model)
    p = d / l
    D = model_correction(p, gamma, model)
    V = variance(p, l, gamma, model)
    return D, V
end

function distance_pairdel(a::BioSequence, b::BioSequence, model::JC69)
    d, l = count_differences_pairdel(a, b, model)
    p = d / l
    D = model_correction(p, model)
    V = variance(p, l, model)
    return D, V
end

function distance_pairdel(a::BioSequence, b::BioSequence, model::JC69, gamma::Float64)
    d, l = count_differences_pairdel(a, b, model)
    p = d / l
    D = model_correction(p, gamma, model)
    V = variance(p, l, gamma, model)
    return D, V
end



# K80 Distance computation.

function distance(a::BioSequence, b::BioSequence, model::K80)
    nd, ns, nv, l = countTsTv(a::BioSequence, b::BioSequence, model)

end
