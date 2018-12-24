__precompile__(true)
module Attitude

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
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012, p.27.
"""
function Rx(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +1  0  0 ;
              0 +c +s ;
              0 -s +c]
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
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012, p.27.
"""
function Ry(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +c  0 -s ;
              0 +1  0 ;
             +s  0 +c]
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
1. O. Montenbruck, and E. Gill, _Satellite Orbits: Models, Methods and 
Applications_, 2012, p.27.
"""
function Rz(angle::Real ; use_degrees=false::Bool)
    
    if use_degrees
        angle *= pi/180.0
    end

    c = cos(angle)
    s = sin(angle)

    return [ +c +s  0 ;
             -s +c  0 ;
              0  0 +1]
end

##############
# Quaternion #
##############

##############
# EulerAngle #
##############

#############
# EulerAxis #
#############

end # End module Coordinates