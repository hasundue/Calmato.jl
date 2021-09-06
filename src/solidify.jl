struct SolidifyResult
    sys::System
    eqs::Vector{EquilibResult}
end

function Base.display(res::SolidifyResult)
    N = length(res.eqs)
    N == 0 && return
    Tl = res.eqs[1].T
    Tu = res.eqs[end].T
    println("Solidification:")
    println("\t$Tl ≤ T ≤ $Tu")
    println("\t$N calculations")
end

function Base.print(res::SolidifyResult)
    println("Solidification Result:\n")
    for eq in res.eqs
        @printf "T = %.2f\n" eq.T
        print(eq)
    end
end

function solidify(sys::System, X::Vector{<:Real} = equiatom(sys),
                  Trange::AbstractRange = defaulttemprange(sys))
    eqs = Vector{Calmato.EquilibResult}()
    for T in Trange
        eq = equilib(sys, X, T)
        push!(eqs, eq)
    end
    return SolidifyResult(sys, eqs)
end

function defaulttemprange(sys::System)
    Tl, Tu = sys.temp
    N::Int = div(Tu - Tl, 100)
    return LinRange(Tl, Tl + 100N, N)
end

@recipe function plot(res::SolidifyResult)
    sys = res.sys
    eqs = res.eqs
    K = sys.nphas
    N = length(eqs)

    names = reshape([ sys.phass[k].name for k in 1:K ], 1,K)
    Ts = [ eqs[n].T for n in 1:N ]
    xss = [ [ eqs[n].Y[k] / sum(eqs[n].Y) for n in 1:N ] for k in 1:K ]

    xguide --> "T (K)"
    yguide --> "mole fraction"
    label --> names

    (Ts, xss)
end