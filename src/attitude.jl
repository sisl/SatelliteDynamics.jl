__precompile__(true)
module Attitude

using LinearAlgebra

#############
# Rotations #
#############

export Rx
"""
Rotation matrix, for a rotation about the x-axis.

Arguments:
- `angle::Real`: Counter-clockwise angle of rotation as viewed looking back along the postive direction of the rotation axis.
- `use_degrees:Bool`: If `true` interpret input as being in degrees.

Returns:
- `r::Array{<:Real, 2}`: Rotation matrix

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.27.
"""
function Rx(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +1.0  0.0  0.0;
              0.0 +c   +s;
              0.0 -s   +c]
end

export Ry
"""
Rotation matrix, for a rotation about the y-axis.

Arguments:
- `angle::Real`: Counter-clockwise angle of rotation as viewed looking back along the postive direction of the rotation axis.
- `use_degrees:Bool`: If `true` interpret input as being in degrees.

Returns:
- `r::Array{<:Real, 2}`: Rotation matrix

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.27.
"""
function Ry(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +c    0.0 -s;
              0.0 +1.0  0.0;
             +s    0.0 +c]
end

export Rz
"""
Rotation matrix, for a rotation about the z-axis.

Arguments:
- `angle::Real`: Counter-clockwise angle of rotation as viewed looking back along the postive direction of the rotation axis.
- `use_degrees:Bool`: If `true` interpret input as being in degrees.

Returns:
- `r::Array{<:Real, 2}`: Rotation matrix

References:
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and Applications_, 2012, p.27.
"""
function Rz(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +c    +s   0.0;
             -s    +c   0.0;
              0.0  0.0 +1.0]
end

###################
# Data Structures #
###################


export Quaternion
"""
The `Quaternion` type defines a _scalar-first_ quaternion for representation
of attitude.

Data members:
- `q0::Float64`: Scalar part of quaternion
- `q1::Float64`: First vector component of quaternion
- `q2::Float64`: Second vector component of quaternion
- `q3::Float64`: Third vector component of quaternion

References:
1. J. Diebel, _Representing attitude: Euler angles, unit quaternions, and rotation vectors._ Matrix 58(15-16) (2006).
"""
mutable struct Quaternion
    q0::Float64
    q1::Float64
    q2::Float64
    q3::Float64

    function Quaternion(q0::Real, q1::Real, q2::Real, q3::Real)
        n = sqrt(q0^2 + q1^2 + q2^2 + q3^2)
        new(q0/n, q1/n, q2/n, q3/n)
    end
end

# export EulerAngle
# """
# The `EulerAngle` type provides a represenation of EulerAngles for storing attitude
# information.

# Valid sequences are: `:121, :123, :131, :132, :212, :213, :231, :232, :312, :313, :321, :323`.

# Data members:
# - `seq::symbol`: Order of application of angles with respect to body axis.
# - `phi::Float64`: First Euler angle
# - `theta::Float64`: Second Euler angle
# - `psi::Float64`: Third Euler angle
#
# References:
# 1. J. Diebel, _Representing attitude: Euler angles, unit quaternions, and rotation vectors._ Matrix 58(15-16) (2006).
# """
# mutable struct EulerAngle
#     seq::symbol
#     phi::Float64
#     theta::Float64
#     psi::Float64
# end

# export EulerAxis
# """
# The `EulerAxis` type provides a representation of the Euler angle-and-axis attitude
# representation.

# Data members:
# - `theta::Float64`: Angle of rotation
# - `v::Array{Float64, 1}`: Axis of rotation
# """
# mutable struct EulerAxis
#     theta::Float64
#     v::Array{Float64, 1}
#
# References:
# 1. J. Diebel, _Representing attitude: Euler angles, unit quaternions, and rotation vectors._ Matrix 58(15-16) (2006).
# end

##############
# Quaternion #
##############

# Quaternion Constructors 
function Quaternion(vec::Array{<:Real, 1})
    if length(vec) != 4
        throw(ArgumentError("Invalid array for Quaternion initialization. Input length: $(length(vec)), Required length: 4"))
    end

    Quaternion(vec...)
end

function Quaternion(mat::Array{<:Real, 2})
    if size(mat) == (1,4)
        # Actually vector initialization. so it and return early
        return Quaternion(mat...)
    elseif size(mat) != (3,3)
        throw(ArgumentError("Invalid array for Quaternion initialization. Input size: $(size(mat)), Required size: (3,3)"))
    end

    temp = zeros(Float64, 4)
    temp[1] = 1 + mat[1, 1] + mat[2, 2] + mat[3, 3]
    temp[2] = 1 + mat[1, 1] - mat[2, 2] - mat[3, 3]
    temp[3] = 1 - mat[1, 1] + mat[2, 2] - mat[3, 3]
    temp[4] = 1 - mat[1, 1] - mat[2, 2] + mat[3, 3]

    # Get the maximum value and its index
    den, ind = findmax(temp)
    den      = sqrt(den) # Short-cut equivalence to reuse test information

    # Select optimal inverse mapping
    q0, q1, q2, q3 = 0.0, 0.0, 0.0, 0.0
    if ind == 1
        q0 = 0.5 * den
        q1 = 0.5 * (mat[2, 3] - mat[3, 2]) / den
        q2 = 0.5 * (mat[3, 1] - mat[1, 3]) / den
        q3 = 0.5 * (mat[1, 2] - mat[2, 1]) / den
    elseif ind == 2
        q0 = 0.5 * (mat[2, 3] - mat[3, 2]) / den
        q1 = 0.5 * den
        q2 = 0.5 * (mat[1, 2] + mat[2, 1]) / den
        q3 = 0.5 * (mat[3, 1] + mat[1, 3]) / den
    elseif ind == 3
        q0 = 0.5 * (mat[3, 1] - mat[1, 3]) / den
        q1 = 0.5 * (mat[1, 2] + mat[2, 1]) / den
        q2 = 0.5 * den
        q3 = 0.5 * (mat[2, 3] + mat[3, 2]) / den
    elseif ind == 4
        q0 = 0.5 * (mat[1, 2] - mat[2, 1]) / den
        q1 = 0.5 * (mat[3, 1] + mat[1, 3]) / den
        q2 = 0.5 * (mat[2, 3] + mat[3, 2]) / den
        q3 = 0.5 * den
    end

    return Quaternion(q0, q1, q2, q3)
end

# Return quaternion as a matrix
# function as_matrix(q::Quaternion)

# end

# Quaternion Operators
function norm(q::Quaternion)
    return sqrt(q.q0^2 + q.q1^2 + q.q2^2 + q.q3^2)
end

# Type conversions


##############
# EulerAngle #
##############

#############
# EulerAxis #
#############

end # End module Coordinates