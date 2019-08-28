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
function Rx(angle::Real ; use_degrees::Bool=false)
    
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
function Ry(angle::Real ; use_degrees::Bool=false)
    
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
function Rz(angle::Real ; use_degrees::Bool=false)
    
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
        
        return new(q0/n, q1/n, q2/n, q3/n)
    end
end

export EulerAngle
"""
The `EulerAngle` type provides a represenation of EulerAngles for storing attitude
information.

Valid sequences are: `121, 123, 131, 132, 212, 213, 231, 232, 312, 313, 321, 323`.

Data members:
- `seq::Integer`: Order of application of angles with respect to body axis.
- `phi::Float64`: First Euler angle
- `theta::Float64`: Second Euler angle
- `psi::Float64`: Third Euler angle

References:
1. J. Diebel, _Representing attitude: Euler angles, unit quaternions, and rotation vectors._ Matrix 58(15-16) (2006).
"""
mutable struct EulerAngle
    seq::Int
    phi::Float64
    theta::Float64
    psi::Float64

    function EulerAngle(seq::Integer, phi::Real, theta::Real, psi::Real)
        if !(seq in [121, 123, 131, 132, 212, 213, 231, 232, 312, 313, 321, 323])
            throw(ArgumentError("Invalid EulerAngle sequence: $seq"))
        end

        return new(seq, phi, theta, psi)
    end
end

export EulerAxis
"""
The `EulerAxis` type provides a representation of the Euler angle-and-axis attitude
representation.

Data members:
- `theta::Float64`: Angle of rotation
- `vec::Array{Float64, 1}`: Axis of rotation

References:
1. J. Diebel, _Representing attitude: Euler angles, unit quaternions, and rotation vectors._ Matrix 58(15-16) (2006).
"""
mutable struct EulerAxis
    angle::Float64
    axis::Array{Float64, 1}

    function EulerAxis(angle::Real, axis::Array{<:Real, 1})
        if length(axis) != 3
            throw(ArgumentError("Invalid array for EulerAxis initialization. Input size: $(size(axis)), Required size: (3,)"))
        end

        return new(angle, axis)
    end
end

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

function Quaternion(e::EulerAngle)
    # Extract Quaternion components
    q0, q1, q2, q3 = 0.0, 0.0, 0.0, 0.0

    # Compute sine and cosine values for angles

    # Get Quaternion components depending on Euler Angle sequence
    # reduce number of trig calls
    c1 = cos(e.phi/2.0)
    c2 = cos(e.theta/2.0)
    c3 = cos(e.psi/2.0)
    s1 = sin(e.phi/2.0)
    s2 = sin(e.theta/2.0)
    s3 = sin(e.psi/2.0)

    # Populate the quaternion
    if e.seq == 121
        q0 = c1*c2*c3 - s1*c2*s3
        q1 = c1*c2*s3 + c2*c3*s1
        q2 = c1*c3*s2 + s1*s2*s3
        q3 = c1*s2*s3 - s1*c3*s2

    elseif e.seq == 123
        q0 =  c1*c2*c3 + s1*s2*s3
        q1 = -c1*s2*s3 + c2*c3*s1
        q2 =  c1*c3*s2 + s1*c2*s3
        q3 =  c1*c2*s3 - s1*c3*s2

    elseif e.seq == 131
        q0 =  c1*c2*c3 - s1*c2*s3
        q1 =  c1*c2*s3 + c2*c3*s1
        q2 = -c1*s2*s3 + s1*c3*s2
        q3 =  c1*c3*s2 + s1*s2*s3

    elseif e.seq == 132
        q0 =  c1*c2*c3 - s1*s2*s3
        q1 =  c1*s2*s3 + c2*c3*s1
        q2 =  c1*c2*s3 + s1*c3*s2
        q3 =  c1*c3*s2 - s1*c2*s3

    elseif e.seq == 212
        q0 =  c1*c2*c3 - s1*c2*s3
        q1 =  c1*c3*s2 + s1*s2*s3
        q2 =  c1*c2*s3 + c2*c3*s1
        q3 = -c1*s2*s3 + s1*c3*s2

    elseif e.seq == 213
        q0 =  c1*c2*c3 - s1*s2*s3
        q1 =  c1*c3*s2 - s1*c2*s3
        q2 =  c1*s2*s3 + c2*c3*s1
        q3 =  c1*c2*s3 + s1*c3*s2

    elseif e.seq == 231
        q0 =  c1*c2*c3 + s1*s2*s3
        q1 =  c1*c2*s3 - s1*c3*s2
        q2 = -c1*s2*s3 + c2*c3*s1
        q3 =  c1*c3*s2 + s1*c2*s3

    elseif e.seq == 232
        q0 =  c1*c2*c3 - s1*c2*s3
        q1 =  c1*s2*s3 - s1*c3*s2
        q2 =  c1*c2*s3 + c2*c3*s1
        q3 =  c1*c3*s2 + s1*s2*s3

    elseif e.seq == 312
        q0 =  c1*c2*c3 + s1*s2*s3
        q1 =  c1*c3*s2 + s1*c2*s3
        q2 =  c1*c2*s3 - s1*c3*s2
        q3 = -c1*s2*s3 + c2*c3*s1

    elseif e.seq == 313
        q0 =  c1*c2*c3 - s1*c2*s3
        q1 =  c1*c3*s2 + s1*s2*s3
        q2 =  c1*s2*s3 - s1*c3*s2
        q3 =  c1*c2*s3 + c2*c3*s1

    elseif e.seq == 321
        q0 =  c1*c2*c3 - s1*s2*s3
        q1 =  c1*c2*s3 + s1*c3*s2
        q2 =  c1*c3*s2 - s1*c2*s3
        q3 =  c1*s2*s3 + c2*c3*s1

    elseif e.seq == 323
        q0 =  c1*c2*c3 - s1*c2*s3
        q1 = -c1*s2*s3 + s1*c3*s2
        q2 =  c1*c3*s2 + s1*s2*s3
        q3 =  c1*c2*s3 + c2*c3*s1
    else
        # Should get an invalid sequence, but it is possible if a user
        # Directly sets the sequence number
        throw(ArgumentError("Invalid EulerAngle sequence: $e.seq"))
    end

    return Quaternion(q0, q1, q2, q3)
end

function Quaternion(e::EulerAxis)
    # Extract Quaternion components
    q0 = cos(e.angle/2.0)
    q1 = e.axis[1]*sin(e.angle/2.0)
    q2 = e.axis[2]*sin(e.angle/2.0)
    q3 = e.axis[3]*sin(e.angle/2.0)

    return Quaternion(q0, q1, q2, q3)
end

#########################
# Quaternion Operations #
#########################

function Base.getindex(q::Quaternion, I::UnitRange{<:Integer})
    # Allocate vector once
    vec = as_vector(q)

    # Return selected index or range
    return [vec[i] for i in I]
end

function Base.getindex(q::Quaternion, I::Integer)
    if I == 1
        return q.q0
    elseif I == 2
        return q.q1
    elseif I == 3
        return q.q2
    elseif I == 4
        return q.q3
    else
        throw(BoundsError())
    end
end

Base.getindex(q::Quaternion, ::Colon) = [q.q0, q.q1, q.q2, q.q3]

# Return quaternion as a vector
export as_vector
"""
Return quaternion as a vector. 

Equivalent to q[:]

Arguments:
- `q::Quaternion`: Quaternion

Returns:
- `vec::Array{Float64, 1}`: Quaternion as a (4,) vector
"""
function as_vector(q::Quaternion)
    return q[:]
end

# Return quaternion as a matrix
export as_matrix
"""
Return the rotation matrix representation of a Quaternion.

Arguments:
- `q::Quaternion`: Quaternion

Returns:
- `mat::Array{Float64, 2}`: Rotation Matrix on SO(3).
"""
function as_matrix(q::Quaternion)
    # initialize Empty Matrix
    mat = zeros(Float64, 3, 3)

    # Construct matrix from Quaternion
    mat[1, 1] = q.q0*q.q0 + q.q1*q.q1 - q.q2*q.q2 - q.q3*q.q3
    mat[1, 2] = 2*q.q1*q.q2 + 2*q.q0*q.q3
    mat[1, 3] = 2*q.q1*q.q3 - 2*q.q0*q.q2
    mat[2, 1] = 2*q.q1*q.q2 - 2*q.q0*q.q3
    mat[2, 2] = q.q0*q.q0 - q.q1*q.q1 + q.q2*q.q2 - q.q3*q.q3
    mat[2, 3] = 2*q.q2*q.q3 + 2*q.q0*q.q1
    mat[3, 1] = 2*q.q1*q.q3 + 2*q.q0*q.q2
    mat[3, 2] = 2*q.q2*q.q3 - 2*q.q0*q.q1
    mat[3, 3] = q.q0*q.q0 - q.q1*q.q1 - q.q2*q.q2 + q.q3*q.q3

    return mat
end

function Base.copy(q::Quaternion)
    return Quaternion(q.q0, q.q1, q.q2, q.q3)
end

function Base.deepcopy(q::Quaternion)
    return Quaternion(q.q0, q.q1, q.q2, q.q3)
end

"""
Compute the norm of a Quaternion.

Equivalent to `sqrt(q0^2 + q1^2 + q2^2 + q3^2)`

Arguments:
- `q::Quaternion`: Quaternion

Returns:
- `q_norm::Float64`: Norm of quaternion.
"""
function LinearAlgebra.norm(q::Quaternion)
    return sqrt(q.q0^2 + q.q1^2 + q.q2^2 + q.q3^2)
end

"""
Normalize a Quaternion in-place.

Equivalent to q=q/norm(q)

Arguments:
- `q::Quaternion`: Quaternion

Returns:
- `q_norm::Float64`: Norm of quaternion.
"""
function LinearAlgebra.normalize(q::Quaternion)
    # Get Quaternion norm
    q_norm = norm(q)

    # Normalize q in-place
    q.q0 = q.q0/q_norm
    q.q1 = q.q1/q_norm
    q.q2 = q.q2/q_norm
    q.q3 = q.q3/q_norm

    # Ensure return value is nothing
    nothing
end

"""
Get conjugate Quaternion.

Arguments:
- `q::Quaternion`: Input Quaternion

Returns:
- `q_conj::Quaternion`: Conjugate Quaternion of input
"""
function Base.conj(q::Quaternion)
    return Quaternion(q.q0, -q.q1, -q.q2, -q.q3)
end

"""
Get Quaternion inverse.

Arguments:
- `q::Quaternion`: Input Quaternion

Returns:
- `q_inv::Quaternion`: Inverse Quaternion of input
"""
function Base.inv(q::Quaternion)
    # Same as Quaternion conjugate since all quaternions are normalized to have
    # unit norm on construction
    return conj(q)
end

function Base.:-(q::Quaternion)
    return Quaternion(-q.q0, -q.q1, -q.q2, -q.q3)
end

function Base.:-(qa::Quaternion, qb::Quaternion)
    return [qa.q0 - qb.q0
            qa.q1 - qb.q1
            qa.q2 - qb.q2
            qa.q3 - qb.q3]
end

function Base.:+(qa::Quaternion, qb::Quaternion)
    return [qa.q0 + qb.q0
            qa.q1 + qb.q1
            qa.q2 + qb.q2
            qa.q3 + qb.q3]
end

function Base.:+(q::Quaternion, n::Real)
    return [q.q0 + n
            q.q1 + n
            q.q2 + n
            q.q3 + n]
end

function Base.:+(n::Real, q::Quaternion)
    return q+n
end

function Base.:*(qa::Quaternion, qb::Quaternion)
    # # Quaternion Multiplication
    # qcos = self.data[0]*other.data[0] - np.dot(self.data[1:4], other.data[1:4])
    # qvec = self.data[0]*other.data[1:4] + other.data[0]*self.data[1:4] + np.cross(self.data[1:4], other.data[1:4])

    # Quaternion Multiplication
    qcos = qa.q0*qb.q0 - dot(qa[2:4], qb[2:4])
    qvec = qa.q0*qb[2:4] + qb.q0*qa[2:4] + cross(qa[2:4], qb[2:4])

    return Quaternion(qcos, qvec...)
end

function Base.:*(q::Quaternion, n::Real)
    return q[:]*n
end

function Base.:*(n::Real, q::Quaternion)
    return q*n
end

export slerp
"""
Perform spherical linear interpolation (SLERP) on two quaternions. Interpolatles 
from quaternion, `q1`, to quaternion, `q2`, at normalized interpolation time, `t`.

Interpolation time must be in the range `[0, 1]` a value of `0` will return `q1`,
while a value of `1` will return `q2`.

Arguments:
- `q1::Quaternion`: Starting Quaternion
- `q2::Quaternion`: Ending Quaternion
- `t::Real`: Normalized interpolation time. [0, 1]

Returns:
- `q:Quaternion`: Quaternion attitude interpolation from q1 toward q2 at time t.
"""
function slerp(q0::Quaternion, q1::Quaternion, t::Real)
    # Check Range on t
    if t < 0.0 || t > 1.0
        throw(ArgumentError("Invalid interpolation time $t. t must be in the range [0, 1]."))
    end

    # Extract vectors and normalize
    q0 = copy(q0)[:]
    q1 = copy(q1)[:]

    # Compute cosine of the angle between the two vectors
    dp = dot(q0, q1)

    # If the dot product is negative, the quaternions have opposite handed-ness 
    # and slerp won't take the shortest path. Fix by reversing one quaternion.
    if dp < 0.0
        q1  = -q1
        dp = -dp
    end

    # If the inputs are too close we use linear interpolation instead
    if dp > 0.9995
        return Quaternion(q0 + (q1 - q0)*t)
    end

    theta0 = acos(dp) # Angle between input vectors
    theta  = theta0*t  # Angle between q0 and result quaternion

    s0 = cos(theta) - dp*sin(theta)/sin(theta0)
    s1 = sin(theta) / sin(theta0)

    return Quaternion((s0 * q0) + (s1 * q1))
end

##############
# EulerAngle #
##############

function EulerAngle(seq::Integer, vec::Array{<:Real, 1})
    if length(vec) != 3
        throw(ArgumentError("Invalid array for EulerAngle initialization. Input length: $(length(vec)), Required length: 3"))
    end

    EulerAngle(seq, vec...)
end

function EulerAngle(seq::Integer, mat::Array{<:Real, 2})
    if size(mat) != (3,3)
        throw(ArgumentError("Invalid array for Quaternion initialization. Input size: $(size(mat)), Required size: (3,3)"))
    end

    # Extract elements out of rotation matrix
    r11 = mat[1, 1]
    r12 = mat[1, 2]
    r13 = mat[1, 3]
    r21 = mat[2, 1]
    r22 = mat[2, 2]
    r23 = mat[2, 3]
    r31 = mat[3, 1]
    r32 = mat[3, 2]
    r33 = mat[3, 3]

    # Select euler angle sequence
    phi, theta, psi = 0.0, 0.0, 0.0
    if seq == 121
        phi   = atan(r21, r31)
        theta = acos(r11)
        psi   = atan(r12, -r13)
    elseif seq == 123
        phi   = atan(r23, r33)
        theta = -asin(r13)
        psi   = atan(r12, r11)
    elseif seq == 131
        phi   = atan(r31, -r21)
        theta = acos(r11)
        psi   = atan(r13, r12)
    elseif seq == 132
        phi   = atan(-r32, r22)
        theta = asin(r12)
        psi   = atan(-r13, r11)
    elseif seq == 212
        phi   = atan(r12, -r32)
        theta = acos(r22)
        psi   = atan(r21, r23)
    elseif seq == 213
        phi   = atan(-r13, r33)
        theta = asin(r23)
        psi   = atan(-r21, r22)
    elseif seq == 231
        phi   = atan(r31, r11)
        theta = -asin(r21)
        psi   = atan(r23, r22)
    elseif seq == 232
        phi   = atan(r32, r12)
        theta = acos(r22)
        psi   = atan(r23, -r21)
    elseif seq == 312
        phi   = atan(r12, r22)
        theta = -asin(r32)
        psi   = atan(r31, r33)
    elseif seq == 313
        phi   = atan(r13, r23)
        theta = acos(r33)
        psi   = atan(r31, -r32)
    elseif seq == 321
        phi   = atan(-r21, r11)
        theta = asin(r31)
        psi   = atan(-r32, r33)
    elseif seq == 323
        phi   = atan(r23, -r13)
        theta = acos(r33)
        psi   = atan(r32, r31)
    else
        throw(ArgumentError("Invalid EulerAngle sequence: $seq"))
    end

    return EulerAngle(seq, phi, theta, psi)
end

function EulerAngle(seq::Integer, q::Quaternion)
    # Construct angle from Quaternion by going through a rotation matrix
    return EulerAngle(seq::Integer, as_matrix(q))
end

function EulerAngle(seq::Integer, ea::EulerAxis)
    # Construct angle from EulerAxis by going through a rotation matrix
    return EulerAngle(seq::Integer, as_matrix(ea))
end

# Access Operators
function Base.getindex(e::EulerAngle, I::UnitRange{<:Integer})
    # Allocate vector once
    vec = as_vector(e)

    # Return selected index or range
    return [vec[i] for i in I]
end

function Base.getindex(e::EulerAngle, I::Integer)
    if I == 1
        return e.phi
    elseif I == 2
        return e.theta
    elseif I == 3
        return e.psi
    else
        throw(BoundsError())
    end
end

Base.getindex(e::EulerAngle, ::Colon) = [e.phi, e.theta, e.psi]

"""
Return Euler angles as a vector.

Equivalent to: `[e.phi, e.theta, e.psi]` for `EulerAngle` `e`

Arguments:
- `e::EulerAngle` Euler Angle

Returns:
- `evec::Array{Float64, 1}` Euler angles components in vector form.
"""
function as_vector(e::EulerAngle)
    return [e.phi, e.theta, e.psi]
end

function as_matrix(e::EulerAngle)
    # Get EulerAngle as matrix by going through Quaternions
    return as_matrix(Quaternion(e))
end

function Base.copy(e::EulerAngle)
    return EulerAngle(e.seq, e.phi, e.theta, e.psi)
end

function Base.deepcopy(e::EulerAngle)
    return EulerAngle(e.seq, e.phi, e.theta, e.psi)
end

#############
# EulerAxis #
#############

function EulerAxis(angle::Real, v1::Real, v2::Real, v3::Real)
    return EulerAxis(angle, [v1, v2, v3])
end

function EulerAxis(vec::Array{<:Real, 1})
    if length(vec) != 4
        throw(ArgumentError("Invalid array for EulerAxis initialization. Input size: $(size(vec)), Required size: (4,)"))
    end

    return EulerAxis(vec[1], [vec[2], vec[3], vec[4]])
end

function EulerAxis(mat::Array{<:Real, 2})
    if size(mat) != (3,3)
        throw(ArgumentError("Invalid array for EulerAxis initialization. Input size: $(size(mat)), Required size: (3,3)"))
    end

    return EulerAxis(Quaternion(mat))
end

function EulerAxis(q::Quaternion)
    # Extract quaternion vector and normalize
    qv = as_vector(q)
    q  = qv/norm(qv)
            
    # Ensure first element is positive
    if q[1] < 0
        q = -q
    end
    
    # Compute Euler Angle
    angle     = 2*acos(q[1])
    qvec_norm = norm(q[2:4])
    vec       = [0.0, 0.0, 0.0]

    if qvec_norm > 1e-15
        vec = q[2:4]/qvec_norm
    end

    return EulerAxis(angle, vec)
end

function EulerAxis(e::EulerAngle)
    # If input is an EulerAngle first compute a quaternion then get the
    # EulerAxis form
    EulerAxis(Quaternion(e))
end

# Access operators

function Base.getindex(e::EulerAxis, I::UnitRange{<:Integer})
    # Allocate vector once
    vec = as_vector(e)

    # Return selected index or range
    return [vec[i] for i in I]
end

function Base.getindex(e::EulerAxis, I::Integer)
    if I == 1
        return e.angle
    elseif I == 2
        return e.axis[1]
    elseif I == 3
        return e.axis[2]
    elseif I == 4
        return e.axis[3]
    else
        throw(BoundsError())
    end
end

Base.getindex(e::EulerAxis, ::Colon) = [e.angle, e.axis[1], e.axis[2], e.axis[3]]

function as_vector(e::EulerAxis)
    return e[:]
end

function as_matrix(e::EulerAxis)
    # Get matrix form from Quaternion for ease
    return as_matrix(Quaternion(e))
end

function Base.copy(e::EulerAxis)
    return EulerAxis(e.angle, e.axis)
end

function Base.deepcopy(e::EulerAxis)
    return EulerAxis(e.angle, e.axis)
end

####################
# Type Conversions #
####################