module State

# Declare-top
abstract type SatelliteState{T} <: Array{<:T, 1} end

struct ECIState{T} <: SatelliteState{T}
    x::Array{T, 1}
end

struct ECEFState{T} <: SatelliteState{T}
    x::Array{T, 1}
end

struct GeodeticState{T} <: SatelliteState{T}
    x::Array{T, 1}
end

struct GeocentricState{T} <: SatelliteState{T}
    x::Array{T, 1}
end

end