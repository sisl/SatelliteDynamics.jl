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
    @test isapprox(v[1], 0.5, atol=tol)
    @test isapprox(v[2], 0.5, atol=tol)
    @test isapprox(v[3], 0.5, atol=tol)
    @test isapprox(v[4], 0.5, atol=tol)

    # Test Access
    # @test isapprox(q[1], 0.5, atol=tol)
    # @test isapprox(q[2], 0.5, atol=tol)
    # @test isapprox(q[3], 0.5, atol=tol)
    # @test isapprox(q[4], 0.5, atol=tol)
end

let
    q = Quaternion(1.0, 1.0, 1.0, 1.0)

    v = as_vector(q)

    tol = 1e-12
    @test isapprox(v[1], 0.5, atol=tol)
    @test isapprox(v[2], 0.5, atol=tol)
    @test isapprox(v[3], 0.5, atol=tol)
    @test isapprox(v[4], 0.5, atol=tol)
end

let
    q = Quaternion(1.0, 0.0, 0.0, 0.0)

    mat = as_matrix(q)

    tol = 1.0e-12
    @test isapprox(mat[1, 1], 1.0, atol=tol)
    @test isapprox(mat[1, 2], 0.0, atol=tol)
    @test isapprox(mat[1, 3], 0.0, atol=tol)
    @test isapprox(mat[2, 1], 0.0, atol=tol)
    @test isapprox(mat[2, 2], 1.0, atol=tol)
    @test isapprox(mat[2, 3], 0.0, atol=tol)
    @test isapprox(mat[3, 1], 0.0, atol=tol)
    @test isapprox(mat[3, 2], 0.0, atol=tol)
    @test isapprox(mat[3, 3], 1.0, atol=tol)
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
    @test isapprox(q.q0, 0.5, atol=tol)
    @test isapprox(q.q1, 0.5, atol=tol)
    @test isapprox(q.q2, 0.5, atol=tol)
    @test isapprox(q.q3, 0.5, atol=tol)

    @test normalize(q) == nothing
end

let
    q = Quaternion(randn(4))

    qc = conj(q)

    tol = 1e-12
    @test isapprox(qc.q0, q.q0, atol=tol)
    @test isapprox(qc.q1, -q.q1, atol=tol)
    @test isapprox(qc.q2, -q.q2, atol=tol)
    @test isapprox(qc.q3, -q.q3, atol=tol)
end

let
    q = Quaternion(randn(4))

    qi = inv(q)

    tol = 1e-12
    @test isapprox(qi.q0, q.q0, atol=tol)
    @test isapprox(qi.q1, -q.q1, atol=tol)
    @test isapprox(qi.q2, -q.q2, atol=tol)
    @test isapprox(qi.q3, -q.q3, atol=tol)
end

###########################
# EulerAngle Constructors #
###########################

##########################
# EulerAxis Constructors #
##########################