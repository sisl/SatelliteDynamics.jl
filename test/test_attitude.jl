let
    deg = 45.0
    rad = deg*pi/180
    r = Rx(deg, use_degrees=true)

    tol = 1e-8
    @test isapprox(r[1, 1], 1.0,       atol=tol)
    @test isapprox(r[1, 2], 0.0,       atol=tol)
    @test isapprox(r[1, 3], 0.0,       atol=tol)

    @test isapprox(r[2, 1], 0.0,       atol=tol)
    @test isapprox(r[2, 2], +cos(rad), atol=tol)
    @test isapprox(r[2, 3], +sin(rad), atol=tol)

    @test isapprox(r[3, 1], 0.0,       atol=tol)
    @test isapprox(r[3, 2], -sin(rad), atol=tol)
    @test isapprox(r[3, 3], +cos(rad), atol=tol)
end

let
    deg = 45.0
    rad = deg*pi/180
    r = Ry(deg, use_degrees=true)

    tol = 1e-8
    @test isapprox(r[1, 1], +cos(rad), atol=tol)
    @test isapprox(r[1, 2], 0.0,       atol=tol)
    @test isapprox(r[1, 3], -sin(rad), atol=tol)

    @test isapprox(r[2, 1], 0.0,       atol=tol)
    @test isapprox(r[2, 2], 1.0,       atol=tol)
    @test isapprox(r[2, 3], 0.0,       atol=tol)

    @test isapprox(r[3, 1], +sin(rad), atol=tol)
    @test isapprox(r[3, 2], 0.0,       atol=tol)
    @test isapprox(r[3, 3], +cos(rad), atol=tol)
end

let
    deg = 45.0
    rad = deg*pi/180
    r = Rz(deg, use_degrees=true)

    tol = 1e-8
    @test isapprox(r[1, 1], +cos(rad), atol=tol)
    @test isapprox(r[1, 2], +sin(rad), atol=tol)
    @test isapprox(r[1, 3], 0.0,       atol=tol)

    @test isapprox(r[2, 1], -sin(rad), atol=tol)
    @test isapprox(r[2, 2], +cos(rad), atol=tol)
    @test isapprox(r[2, 3], 0.0,       atol=tol)

    @test isapprox(r[3, 1], 0.0,       atol=tol)
    @test isapprox(r[3, 2], 0.0,       atol=tol)
    @test isapprox(r[3, 3], 1.0,       atol=tol)
end

###########################
# Quaternion Constructors #
###########################

let
    # Vector initialization
    q = Quaternion([1.1, 0.0, 0.0, 0.0])

    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # Test vector initialization errors
    @test_throws ArgumentError Quaternion([1 2 3 4 5])
end

let
    # Matrix initialization Quaternion Fomulation 1
    mat = [1.0 0.0 0.0;
           0.0 1.0 0.0;
           0.0 0.0 1.0]

    q = Quaternion(mat)

    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # Matrix initialization Quaternion Fomulation 2
    mat = [1.0 0.0 0.0;
           0.0 -1.0 0.0;
           0.0 0.0 -1.0]

    q = Quaternion(mat)

    tol = 1e-12
    @test isapprox(q.q0, 0.0, atol=tol)
    @test isapprox(q.q1, 1.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # Matrix initialization Quaternion Fomulation 3
    mat = [-1.0 0.0 0.0;
           0.0 1.0 0.0;
           0.0 0.0 -1.0]

    q = Quaternion(mat)

    tol = 1e-12
    @test isapprox(q.q0, 0.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 1.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # Matrix initialization Quaternion Fomulation 4
    mat = [-1.0 0.0 0.0;
           0.0 -1.0 0.0;
           0.0 0.0 1.0]

    q = Quaternion(mat)

    tol = 1e-12
    @test isapprox(q.q0, 0.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 1.0, atol=tol)
end

let
    # EulerAngle initialization - 121
    e = EulerAngle(121, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 123
    e = EulerAngle(123, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)
    
    # EulerAngle initialization - 131
    e = EulerAngle(131, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 132
    e = EulerAngle(132, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 212
    e = EulerAngle(212, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 213
    e = EulerAngle(213, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 231
    e = EulerAngle(231, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 232
    e = EulerAngle(232, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 312
    e = EulerAngle(312, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 313
    e = EulerAngle(313, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 321
    e = EulerAngle(321, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)

    # EulerAngle initialization - 323
    e = EulerAngle(323, 0.0, 0.0, 0.0)

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)
end

let 
    # EulerAngle initialization - 323
    e = EulerAngle(323, 0.0, 0.0, 0.0)

    # Change to 999 to confirm error
    e.seq = 999

    @test typeof(e) == EulerAngle
    @test_throws ArgumentError Quaternion(e)
end

let
    # EulerAxis initialization
    e = EulerAxis(0.0, [0.0, 0.0, 0.0])

    q = Quaternion(e)
    tol = 1e-12
    @test isapprox(q.q0, 1.0, atol=tol)
    @test isapprox(q.q1, 0.0, atol=tol)
    @test isapprox(q.q2, 0.0, atol=tol)
    @test isapprox(q.q3, 0.0, atol=tol)
end


#########################
# Quaternion Operations #
#########################

let
    q = Quaternion(1.0, 1.0, 1.0, 1.0)

    v = q[:]

    tol = 1e-12
    array_isapprox(v, 0.5, atol=tol)

    # Test Access
    @test isapprox(q[1], 0.5, atol=tol)
    @test isapprox(q[2], 0.5, atol=tol)
    @test isapprox(q[3], 0.5, atol=tol)
    @test isapprox(q[4], 0.5, atol=tol)

    @test length(q[2:3]) == 2

    @test_throws BoundsError q5 = q[5]
end

let
    q = Quaternion(1.0, 1.0, 1.0, 1.0)

    v = as_vector(q)

    tol = 1e-12
    array_isapprox(v, 0.5, atol=tol)
end

let
    q = Quaternion(1.0, 0.0, 0.0, 0.0)

    mat = as_matrix(q)

    tol = 1.0e-12
    array_isapprox(mat, one(mat), atol=tol)
end

let
    q1 = Quaternion(randn(4))
    q2 = copy(q1)

    @test pointer_from_objref(q1) != pointer_from_objref(q2)
end

let
    q1 = Quaternion(randn(4))
    q2 = deepcopy(q1)

    @test pointer_from_objref(q1) != pointer_from_objref(q2)
end

let
    q = Quaternion([1.0 0.0 0.0 0.0])

    # Circumvent normalization during Quaternion initialization
    q.q0 = 1.0
    q.q1 = 1.0
    q.q2 = 1.0
    q.q3 = 1.0

    # Norm
    @test norm(q) == 2.0

    # Normalize Quaternion
    normalize(q)
    
    tol = 1e-12
    array_isapprox(q[:], 0.5, atol=tol)

    @test normalize(q) == nothing

    @test_throws ArgumentError Quaternion([1 2 3 4 5])
end

let
    q = Quaternion(randn(4))

    qc = conj(q)

    tol = 1e-12
    @test isapprox(qc.q0,  q.q0, atol=tol)
    @test isapprox(qc.q1, -q.q1, atol=tol)
    @test isapprox(qc.q2, -q.q2, atol=tol)
    @test isapprox(qc.q3, -q.q3, atol=tol)
end

let
    q = Quaternion(randn(4))

    qi = inv(q)

    tol = 1e-12
    @test isapprox(qi.q0,  q.q0, atol=tol)
    @test isapprox(qi.q1, -q.q1, atol=tol)
    @test isapprox(qi.q2, -q.q2, atol=tol)
    @test isapprox(qi.q3, -q.q3, atol=tol)
end

let
    q  = Quaternion(randn(4))
    q2 = -q

    tol = 1e-12
    @test isapprox(q2.q0, -q.q0, atol=tol)
    @test isapprox(q2.q1, -q.q1, atol=tol)
    @test isapprox(q2.q2, -q.q2, atol=tol)
    @test isapprox(q2.q3, -q.q3, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])
    q2 = Quaternion([1.0 1.0 1.0 1.0])

    q = q2 - q1

    tol = 1e-12
    array_isapprox(q, 0.0, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])
    q2 = Quaternion([1.0 1.0 1.0 1.0])

    q = q2 + q1

    tol = 1e-12
    array_isapprox(q, 1.0, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])

    q = q1 + 1.5

    tol = 1e-12
    array_isapprox(q, 2.0, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])

    q = 1.5 + q1

    tol = 1e-12
    array_isapprox(q, 2.0, atol=tol)
end

let
    q1 = Quaternion([1.0 0.0 0.0 0.0])
    q2 = Quaternion([1.0 0.0 0.0 0.0])

    qp = q2*q1

    tol = 1e-12
    isapprox(qp.q0, 1.0, atol=tol)
    isapprox(qp.q1, 0.0, atol=tol)
    isapprox(qp.q2, 0.0, atol=tol)
    isapprox(qp.q3, 0.0, atol=tol)

    qp = q1*q2

    tol = 1e-12
    isapprox(qp.q0, 1.0, atol=tol)
    isapprox(qp.q1, 0.0, atol=tol)
    isapprox(qp.q2, 0.0, atol=tol)
    isapprox(qp.q3, 0.0, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])

    q = q1*4

    tol = 1e-12
    array_isapprox(q, 2.0, atol=tol)
end

let
    q1 = Quaternion([1.0 1.0 1.0 1.0])

    q = 4*q1

    tol = 1e-12
    array_isapprox(q, 2.0, atol=tol)
end

let
    q1  = Quaternion([1.0 0.0 0.0 0.0])
    q2  = Quaternion([0.0 0.0 0.0 1.0])
    tol = 1e-12

    q = slerp(q1, q2, 0.0)
    array_isapprox(q[:], [1.0 0.0 0.0 0.0], atol=tol)

    q = slerp(q1, q2, 1.0)
    array_isapprox(q[:], [0.0 0.0 0.0 1.0], atol=tol)

    q = slerp(q1, q2, 0.5)
    array_isapprox(q[:], [sqrt(2)/2.0 0.0 0.0 sqrt(2)/2.0], atol=tol)

    # Edge Cases
    q1  = Quaternion([1.0 0.0 0.0 0.0])
    q2  = Quaternion([-1.0 0.0 0.0 0.0])

    q = slerp(q1, q2, 0.0)
    array_isapprox(q[:], [1.0 0.0 0.0 0.0], atol=tol)
end

###########################
# EulerAngle Constructors #
###########################

let
    # Test invalid constructors
    @test_throws ArgumentError EulerAngle(999, 0.1, 0.2, 0.3)
    @test_throws ArgumentError EulerAngle(123, [0.1, 0.2, 0.3, 0.4])
    @test_throws ArgumentError EulerAngle(123, zeros(Float64, 4, 4))
end

let
    # Test valid constructors
    e = EulerAngle(123, 0.1, 0.2, 0.3)

    tol = 1.0e-12
    @test isapprox(e.phi, 0.1, atol=tol)
    @test isapprox(e.theta, 0.2, atol=tol)
    @test isapprox(e.psi, 0.3, atol=tol)

    # Test valid constructors
    e = EulerAngle(123, [0.4, 0.5, 0.6])

    tol = 1.0e-12
    @test isapprox(e.phi, 0.4, atol=tol)
    @test isapprox(e.theta, 0.5, atol=tol)
    @test isapprox(e.psi, 0.6, atol=tol)
end

# Test alternate constructors
let
    e = EulerAngle(123, Quaternion(-1.0, 0.0, 0.0, 0.0))
    tol = 1.0e-12
    array_isapprox((e[:]), 0.0, atol=tol)

    e = EulerAngle(123, EulerAxis(0.0, 0.1, 0.2, 0.3))
    tol = 1.0e-12
    array_isapprox((e[:]), 0.0, atol=tol)
end

let
    # Test error cases
    mat = zeros(Float64, 4, 4)
    @test_throws ArgumentError EulerAngle(123, mat)

    mat = [1.0 0.0 0.0;
           0.0 1.0 0.0;
           0.0 0.0 1.0]
    tol = 1.0e-12
    
    @test_throws ArgumentError EulerAngle(999, mat)

    # Test angle construction    
    e = EulerAngle(121, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, pi, atol=tol)

    e = EulerAngle(123, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(131, mat)
    @test isapprox(e.phi, pi, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(132, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(212, mat)
    @test isapprox(e.phi, pi, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(213, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(231, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(232, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, pi, atol=tol)

    e = EulerAngle(312, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(313, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, pi, atol=tol)

    e = EulerAngle(321, mat)
    @test isapprox(e.phi, 0.0, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)

    e = EulerAngle(323, mat)
    @test isapprox(e.phi, pi, atol=tol)
    @test isapprox(e.theta, 0.0, atol=tol)
    @test isapprox(e.psi, 0.0, atol=tol)
end

# Test Operators
let
    e = EulerAngle(123, 0.1, 0.2, 0.3)

    v = e[:]

    tol = 1e-12
    array_isapprox(e[:], [0.1, 0.2, 0.3], atol=tol)

    # Test Access
    @test isapprox(e[1], 0.1, atol=tol)
    @test isapprox(e[2], 0.2, atol=tol)
    @test isapprox(e[3], 0.3, atol=tol)

    @test length(e[2:3]) == 2

    @test_throws BoundsError e[4]
end

let
    q = EulerAngle(123, 0.5, 0.5, 0.5)

    v = as_vector(q)

    tol = 1.0e-12
    array_isapprox(v, 0.5, atol=tol)
end

let
    mat = as_matrix(EulerAngle(123, 0.0, 0.0, 0.0))

    tol = 1.0e-12
    array_isapprox(mat, one(mat), atol=tol)
end

let
    e1 = EulerAngle(123, randn(3))
    e2 = copy(e1)

    @test pointer_from_objref(e1) != pointer_from_objref(e2)
end

let
    e1 = EulerAngle(123, randn(3))
    e2 = deepcopy(e1)

    @test pointer_from_objref(e1) != pointer_from_objref(e2)
end


##########################
# EulerAxis Constructors #
##########################

# Constructor Exceptions
let
    @test_throws ArgumentError EulerAxis([0.1, 0.2, 0.3])
    @test_throws ArgumentError EulerAxis([0.1, 0.2, 0.3])
    @test_throws ArgumentError EulerAxis(zeros(Float64, 4, 4))
end

# Vector Constructor
let
    e = EulerAxis([0.1, 0.2, 0.3, 0.4])

    tol = 1.0e-12
    array_isapprox(e[:], [0.1, 0.2, 0.3, 0.4], atol=tol)

    e = EulerAxis(0.1, [0.2, 0.3, 0.4])

    tol = 1.0e-12
    array_isapprox(e[:], [0.1, 0.2, 0.3, 0.4], atol=tol)
end

# Euler Angle constructor
let
    e = EulerAxis(EulerAngle(123, 0.0, 0.0, 0.0))

    tol = 1.0e-12
    array_isapprox(e[:], 0.0, atol=tol)
end

# Quaternion constructor
let
    mat = [1.0 0.0 0.0;
           0.0 1.0 0.0;
           0.0 0.0 1.0]
    e = EulerAxis(mat)

    tol = 1.0e-12
    array_isapprox(e[:], 0.0, atol=tol)
end

# Matrix constructor
let
    q = Quaternion(1.0, 0.0, 0.0, 0.0)

    e = EulerAxis(q)

    tol = 1.0e-12
    array_isapprox(e[:], 0.0, atol=tol)
end

# Test Operators
let
    e = EulerAxis(0.4, 0.1, 0.2, 0.3)

    v = e[:]

    tol = 1e-12
    array_isapprox((e[:]), [0.4, 0.1, 0.2, 0.3], atol=tol)

    # Test Access
    @test isapprox(e[1], 0.4, atol=tol)
    @test isapprox(e[2], 0.1, atol=tol)
    @test isapprox(e[3], 0.2, atol=tol)
    @test isapprox(e[4], 0.3, atol=tol)

    @test length(e[2:3]) == 2

    @test_throws BoundsError e[5]
end

let
    e = EulerAxis(0.5, 0.5, 0.5, 0.5)

    v = as_vector(e)

    tol = 1.0e-12
    array_isapprox(v, 0.5, atol=tol)
end

let
    e = EulerAxis(0.0, 0.0, 0.0, 0.0)

    mat = as_matrix(e)

    tol = 1.0e-12
    array_isapprox(mat, one(mat), atol=tol)
end

let
    e1 = EulerAxis(randn(4))
    e2 = copy(e1)

    @test pointer_from_objref(e1) != pointer_from_objref(e2)
end

let
    e1 = EulerAxis(randn(4))
    e2 = deepcopy(e1)

    @test pointer_from_objref(e1) != pointer_from_objref(e2)
end
