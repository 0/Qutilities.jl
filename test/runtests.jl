using Qutilities
using Base.Test

@testset "Qutilities" begin

@testset "sigma_x, sigma_y, sigma_z" begin
    @test sigma_x^2 == eye(2)
    @test sigma_y^2 == eye(2)
    @test sigma_z^2 == eye(2)
    @test -im * sigma_x * sigma_y * sigma_z == eye(2)
end

@testset "ptrace, ptranspose" begin
    let
        A = eye(Int, 2)
        B = reshape(1:16, 4, 4)
        C = [[0, 0, 1] [0, 1, 0] [1, 0, 0]]

        ABC = kron(A, B, C)
        dims = 2, 4, 3

        @test ptrace(ABC, dims, 1) == trace(A) * kron(B, C)
        @test ptranspose(ABC, dims, 1) == kron(A', B, C)
        @test ptrace(ABC, dims, 2) == trace(B) * kron(A, C)
        @test ptranspose(ABC, dims, 2) == kron(A, B', C)
        @test ptrace(ABC, dims, 3) == trace(C) * kron(A, B)
        @test ptranspose(ABC, dims, 3) == kron(A, B, C')
    end

    let M = reshape(1.:16., 4, 4)
        MT1 = [[1., 2., 9., 10.] [5., 6., 13., 14.] [3., 4., 11., 12.] [7., 8., 15., 16.]]
        MT2 = [[1., 5., 3., 7.] [2., 6., 4., 8.] [9., 13., 11., 15.] [10., 14., 12., 16.]]

        @test ptrace(M, (1, 4), 1) == M
        @test ptranspose(M, (1, 4), 1) == M
        @test ptrace(M, (1, 4), 2) == fill(trace(M), 1, 1)
        @test ptranspose(M, (1, 4), 2) == transpose(M)
        @test ptrace(M, (2, 2), 1) == [[12., 14.] [20., 22.]]
        @test ptranspose(M, (2, 2), 1) == MT1
        @test ptrace(M, (2, 2), 2) == [[7., 11.] [23., 27.]]
        @test ptranspose(M, (2, 2), 2) == MT2
        @test ptrace(M, (4, 1), 1) == fill(trace(M), 1, 1)
        @test ptranspose(M, (4, 1), 1) == transpose(M)
        @test ptrace(M, (4, 1), 2) == M
        @test ptranspose(M, (4, 1), 2) == M

        @test ptrace(M, 1) == ptrace(M, (2, 2), 1)
        @test ptrace(M, 2) == ptrace(M, (2, 2), 2)

        @test ptrace(M) == ptrace(M, (2, 2), 2)
        @test ptranspose(M) == ptranspose(M, (2, 2), 2)
    end

    let M = [[1., im] [im, 1.]]
        @test ptrace(M, (1, 2), 1) == M
        @test ptranspose(M, (1, 2), 1) == M
        @test ptrace(M, (1, 2), 2) == fill(trace(M), 1, 1)
        @test ptranspose(M, (1, 2), 2) == transpose(M)
    end
end

@testset "binent" begin
    @test binent(0.) == 0.
    @test binent(0.5) == 1.
    @test binent(1) == 0.
end

@testset "purity, S_vn, S_renyi" begin
    let
        rho1 = eye(4) / 4.
        rho2 = [[2., im] [-im, 2.]] / 2.
        eigs2 = [1., 3.] / 2.

        @test purity(rho1) == 0.25
        @test purity(rho2) == sum(eigs2.^2)

        @test S_renyi(rho1, 0) == 2.
        @test S_renyi(rho2, 0) == 1.

        @test S_vn(rho1) == 2.
        @test S_vn(rho2) == -sum(eigs2 .* log2(eigs2))

        @test S_renyi(rho1) == 2.
        @test S_renyi(rho2) == -log2(sum(eigs2.^2))

        @test S_renyi(rho1, Inf) == 2.
        @test S_renyi(rho2, Inf) == -log2(maximum(eigs2))
    end
end

@testset "mutinf" begin
    let rho = diagm([3, 2, 1, 2]) / 8.
        @test isapprox(mutinf(rho), (1.5 - 5. * log2(5.) / 8.))
        @test isapprox(mutinf(rho, S_renyi), (1. + 2. * log2(3.) - log2(17.)))
    end
end


@testset "spinflip, concurrence, concurrence_lb, formation, negativity" begin
    let rho = reshape(1.:16., 4, 4)
        rho_f = [[16., -15., -14., 13.] [-12., 11., 10., -9.] [-8., 7., 6., -5.] [4., -3., -2., 1.]]

        @test spinflip(rho) == rho_f
    end

    let rho = zeros(4, 4)
        rho[1, 1] = 0.5
        rho[4, 4] = 0.5
        C = concurrence(rho)

        @test spinflip(rho) == rho
        @test C == 0.
        @test concurrence_lb(rho) == 0.
        @test formation(C) == 0.
        @test negativity(rho) == 0.
    end

    let rho = zeros(4, 4)
        rho[1, 1] = 0.125
        for i=2:3, j=2:3
            rho[i, j] = 0.375
        end
        rho[4, 4] = 0.125
        C = concurrence(rho)
        a = (2 + sqrt(3)) / 4.
        b = (2 - sqrt(3)) / 4.

        @test spinflip(rho) == rho
        @test C == 0.5
        @test concurrence_lb(rho) == sqrt(3.)/4.
        @test formation(C) == -a * log2(a) - b * log2(b)
        @test negativity(rho) == log2(1.5)
    end

    let rho = zeros(4, 4)
        for i=2:3, j=2:3
            rho[i, j] = 0.5
        end
        C = concurrence(rho)

        @test spinflip(rho) == rho
        @test C == 1.
        @test concurrence_lb(rho) == 1.
        @test formation(C) == 1.
        @test negativity(rho) == 1.
    end

    let rho = Array(Complex128, 4, 4)
        for i=1:4
            for j=1:4
                if i < j
                    rho[i, j] = 0.0625im
                elseif i == j
                    rho[i, j] = 0.25
                else
                    rho[i, j] = -0.0625im
                end
            end
        end
        rho_f = copy(transpose(rho))
        rho_f[1, 4] *= -1
        rho_f[2, 3] *= -1
        rho_f[3, 2] *= -1
        rho_f[4, 1] *= -1
        C = concurrence(rho)

        @test spinflip(rho) == rho_f
        @test C == 0.
        @test concurrence_lb(rho) == 0.
        @test formation(C) == 0.
        @test isapprox(negativity(rho), 0., atol=1e-15)
    end
end

end
