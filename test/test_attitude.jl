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