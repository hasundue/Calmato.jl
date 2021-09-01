struct System
    elem::Vector{Element}
    phas::Vector{Phase}
end

function init_system(db::Database, elem::Vector{Element}, phas::Vector{<:Phase})
    for fn in db.func
        eval(Meta.parse(fn.func))
    end
    for ph in phas
        for pr in ph.para
            eval(Meta.parse(pr.func))
        end
    end
    return System(elem, phas)
end

function init_system(db::Database)
    elem = Vector{Element}()
    phas = Vector{Phase}()

    for ph in db.phas
        if 'O' in ph.model
            @warn """
            Order/disorder transition is not supported currently.
            Phase "$(ph.name)" is not included in the system.
            """
            continue
        end
        push!(phas, ph)
    end

    for el in db.elem
        el.name in ["/-", "VA"] && continue
        push!(elem, el)
    end

    return init_system(db, elem, phas)
end

function equiatom(sys::System)
    N = length(sys.elem)
    return fill(1//N, N)
end

function Base.display(sys::System)
    println("System:")
    nelem = length(sys.elem)
    println("\tElements: $nelem")
    nphas = length(sys.phas)
    println("\tPhases: $nphas")
end