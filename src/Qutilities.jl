module Qutilities

export
    binent,
    concurrence,
    concurrence_lb,
    formation,
    mutinf,
    negativity,
    ptrace,
    ptranspose,
    purity,
    sigma_x,
    sigma_y,
    sigma_z,
    spinflip,
    S_renyi,
    S_vn


# All logarithms are base 2.
const LOG = log2

"""
Make Hermitian, but carefully.
"""
function hermitize(rho::AbstractMatrix)
    # Only square matrices are supported.
    size(rho, 1) == size(rho, 2) || throw(DomainError())

    err = maximum(abs(rho - rho'))
    if err > 1e-12
        warn("Matrix is not Hermitian: $(err)")
    end

    # Make the diagonal strictly real.
    rho_H = copy(rho)
    for i=1:size(rho_H, 1)
        rho_H[i, i] = real(rho_H[i, i])
    end

    Hermitian(rho_H)
end

"""
Change negative values to zero.
"""
nonneg(A::Real) = A < 0 ? zero(A) : A
nonneg(A::AbstractArray) = map(nonneg, A)

"""
Shannon entropy.
"""
shannon(xs) = -sum([x * LOG(x) for x in xs[xs .> 0]])


# Single-qubit Pauli matrices.
const sigma_x = [[0, 1] [1, 0]]
const sigma_y = [[0, im] [-im, 0]]
const sigma_z = [[1, 0] [0, -1]]

"""
Partial trace.

If no details are provided, it splits rho into two halves and traces over the
second half.
"""
function ptrace{T}(rho::AbstractMatrix{T}, dims, which::Int)
    # Only square matrices are supported.
    size(rho) == (prod(dims), prod(dims)) || throw(DomainError())

    size_before = prod(dims[1:which-1])
    size_at = dims[which]
    size_after = prod(dims[which+1:end])

    result = zeros(T, size_before*size_after, size_before*size_after)
    for i1=1:size_before
        for j1=1:size_before
            for k=1:size_at
                for i2=1:size_after
                    for j2=1:size_after
                        row1 = size_after * (i1 - 1) + i2
                        col1 = size_after * (j1 - 1) + j2
                        row2 = size_at * size_after * (i1 - 1) + size_after * (k - 1) + i2
                        col2 = size_at * size_after * (j1 - 1) + size_after * (k - 1) + j2
                        result[row1, col1] += rho[row2, col2]
                    end
                end
            end
        end
    end
    result
end

function ptrace(rho::AbstractMatrix, which::Int=2)
    size(rho, 1) % 2 == 0 || throw(DomainError())

    s = div(size(rho, 1), 2)
    ptrace(rho, (s, s), which)
end

"""
Partial transpose.

If no details are provided, it splits rho into two halves and transposes the
second half.
"""
function ptranspose(rho::AbstractMatrix, dims, which::Int)
    # Only square matrices are supported.
    size(rho) == (prod(dims), prod(dims)) || throw(DomainError())

    size_before = prod(dims[1:which-1])
    size_at = dims[which]
    size_after = prod(dims[which+1:end])

    result = similar(rho)
    for i1=1:size_before
        for j1=1:size_before
            for i2=1:size_at
                for j2=1:size_at
                    for i3=1:size_after
                        for j3=1:size_after
                            row1 = size_at * size_after * (i1 - 1) + size_after * (i2 - 1) + i3
                            col1 = size_at * size_after * (j1 - 1) + size_after * (j2 - 1) + j3
                            row2 = size_at * size_after * (i1 - 1) + size_after * (j2 - 1) + i3
                            col2 = size_at * size_after * (j1 - 1) + size_after * (i2 - 1) + j3
                            result[row1, col1] = rho[row2, col2]
                        end
                    end
                end
            end
        end
    end
    result
end

function ptranspose(rho::AbstractMatrix)
    size(rho, 1) % 2 == 0 || throw(DomainError())

    s = div(size(rho, 1), 2)
    ptranspose(rho, (s, s), 2)
end

"""
Binary entropy.
"""
binent(x::Real) = shannon([x, one(x) - x])

"""
Purity.
"""
purity(rho::AbstractMatrix) = rho^2 |> trace |> real

"""
Von Neumann entropy.
"""
S_vn(rho::AbstractMatrix) = rho |> hermitize |> eigvals |> shannon

"""
Rényi entropy.

The alpha parameter may have any value on [0, Inf] except 1. It defaults to 2.
"""
function S_renyi(rho::AbstractMatrix, alpha::Real=2)
    E = rho |> hermitize |> eigvals
    alpha == Inf && return E |> maximum |> LOG |> -
    LOG(sum(E.^alpha)) / (1 - alpha)
end

"""
Mutual information.
"""
function mutinf(rho::AbstractMatrix, S::Function=S_vn)
    S(ptrace(rho, 1)) + S(ptrace(rho, 2)) - S(rho)
end

"""
Wootters spin-flip operation for two qubits in the standard basis.

Ref: Wootters, W. K. (1998). Entanglement of formation of an arbitrary state of
two qubits. Physical Review Letters, 80(10), 2245.
"""
function spinflip(rho::AbstractMatrix)
    size(rho) == (4, 4) || throw(DomainError())

    Y = kron(sigma_y, sigma_y)
    Y * conj(rho) * Y
end

"""
Concurrence of a mixed state for two qubits in the standard basis.

Ref: Wootters, W. K. (1998). Entanglement of formation of an arbitrary state of
two qubits. Physical Review Letters, 80(10), 2245.
"""
function concurrence(rho::AbstractMatrix)
    size(rho) == (4, 4) || throw(DomainError())

    E = rho * spinflip(rho) |> eigvals
    if any(imag(E) .> 1e-15)
        warn("Complex eigenvalues: $(maximum(imag(E)))")
    end
    if any(real(E) .< -1e-12)
        warn("Negative eigenvalues: $(minimum(real(E)))")
    end
    F = E |> real |> nonneg |> sqrt |> sort
    nonneg(F[end] - sum(F[1:end-1]))
end

"""
Lower bound on the concurrence for two qubits in the standard basis.

Ref: Mintert, F., & Buchleitner, A. (2007). Observable entanglement measure for
mixed quantum states. Physical Review Letters, 98(14), 140505.
"""
function concurrence_lb(rho::AbstractMatrix)
    size(rho) == (4, 4) || throw(DomainError())

    2. * (purity(rho) - purity(ptrace(rho))) |> nonneg |> sqrt
end

"""
Entanglement of formation for two qubits (given the concurrence).

Ref: Wootters, W. K. (1998). Entanglement of formation of an arbitrary state of
two qubits. Physical Review Letters, 80(10), 2245.
"""
formation(C::Real) = binent(0.5 * (1. + sqrt(1. - C^2)))
formation(rho::AbstractMatrix) = rho |> concurrence |> formation

"""
Logarithmic negativity for a symmetric bipartition.

Ref: Plenio, M. B. (2005). Logarithmic negativity: A full entanglement monotone
that is not convex. Physical Review Letters, 95(9), 090503.
"""
negativity(rho::AbstractMatrix) = rho |> ptranspose |> svdvals |> sum |> LOG

end
