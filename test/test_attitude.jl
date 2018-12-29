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

####################
# Quaternion Tests #
####################

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
    q = Quaternion([1.0 2.0 3.0 4.0])
    
end