## Plots.jl - plot partition
using RecipesBase
@recipe function f(rs::Vector{<:AbstractRegion}; z_id=0, fix_ids = zeros(0), free_ids = zeros(0), 
        fix_ids = zeros(0), fix_vals=zeros(0), color_mapping=r->r.iter)
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)
    if haskey(plotattributes, :CR_attr)
        z_id,free_ids,fix_ids,fix_vals = pop!(plotattributes,:CR_attr)
    end
    if isempty(free_ids)
        ids = isempty(fix_ids) ? collect(3:nth) : fix_ids
        values = isempty(fix_vals) ? zeros(nth-2) : fix_vals
        free_ids = setdiff(1:nth,ids)
    elseif length(free_ids) != 2 
        error("The number of parameters to plot needs to be 2, not $(length(free_ids))")
    else
        ids = setdiff(1:nth,free_ids) 
        values = isempty(fix_vals) ? zeros(nth-2) : fix_vals
    end

    plotattributes[:CR_attr] = (z_id,free_ids,fix_ids,fix_vals)

    xlims --> (-1,1)
    ylims --> (-1,1)
    xlabel --> "\\theta [$(free_ids[1])]"
    ylabel --> "\\theta [$(free_ids[2])]"
    colorbar --> true
    title --> "Number of iterations" 
    fill_z --> [color_mapping(r) for r in rs] 

    ps = [Polyhedron(slice(r.Ath,r.bth,ids;values)...) for r in rs]
    return ps
end

## PGFPlotsX - plot partition
function pplot(rs::Vector{<:AbstractRegion};key=nothing, fix_ids = nothing, fix_vals=nothing,opts=Dict{Symbol,Any}(), clabel=nothing)
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)
    ids = isnothing(fix_ids) ? collect(3:nth) : fix_ids
    values = isnothing(fix_vals) ? zeros(nth-2) : fix_vals
    free_ids = setdiff(1:nth,ids)

    if(isnothing(key))
        mapping = x->getfield(x,:iter)
        if(isnothing(clabel))
            clabel = "\\# of iterations"
        end
    else
        mapping = key
    end

    ps = PolyDAQP.Polyhedron[]
    cs = typeof(mapping(rs[1]))[]
    for r in rs
        p = Polyhedron(slice(r.Ath,r.bth,ids;values)...)
        isempty(p) && continue
        push!(ps,minrep(p))
        push!(cs,mapping(r))
    end
    lopts = Dict(
                 :xlabel=>"\\large\$\\theta_"*string(free_ids[1])*"\$",
                 :ylabel=>"\\large\$\\theta_"*string(free_ids[2])*"\$",
                 :xlabel_style=> "{yshift={15pt}}",
                 :ylabel_style=> "{yshift={-20pt}}",
                 :xtick=>"{-1,1}",
                 :ytick=>"{-1,1}",
                 :xmin=>-1, :xmax=>1,
                 :ymin=>-1, :ymax=>1,
                )

    if clabel isa String
        push!(lopts, :colorbar_style => "{xlabel ="*clabel*"}")
    end
    opts = merge(lopts,opts)
    PolyDAQP.pplot(ps;cs,opts)
end
## Print ASs in a readable way
function print_ASs(ASs::BitMatrix)
    # TODO: will not work if AS0 ≂̸ ∅
    inds = collect(1:size(ASs,1))
    AS = inds[ASs[:,1]]
    printstyled("$(lpad("ini",3," ")) "; color = :yellow)
    for i in 1:size(ASs,2)-1
        println("$AS ")
        sumi, sumip1= sum(ASs[:,i]), sum(ASs[:,i+1])

        if(sumi < sumip1) # Addition
            add_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            push!(AS,add_ind[1])
            printstyled("+$(lpad(add_ind[1],2," ")) "; color = :green)
        elseif (sumi > sumip1)# Removal
            rm_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            deleteat!(AS,findfirst(AS.==rm_ind[1]))
            printstyled("-$(lpad(rm_ind[1],2," ")) "; color = :red)
        else # Exchange
            diff_inds = .!(ASs[:,i].⊻ .!ASs[:,i+1]);
            rm_ind = inds[ASs[:,i].&& diff_inds];
            add_ind = inds[ASs[:,i+1].&& diff_inds];
            printstyled("$(lpad(rm_ind[1],2," ")) "; color = :red)
            print("-> ")
            printstyled("$(lpad(add_ind[1],2," ")) "; color = :green)
            deleteat!(AS,findfirst(AS.==rm_ind[1]))
            push!(AS,add_ind[1])
        end
    end
    println("$AS")
end
## Display
function Base.:display(r::Region)
    println("=====================================")
    println("Iter : $(r.iter)")
    println("State: $(r.state)")
    println("AS   : $(r.AS)")
    if(!isnan(last(r.chebyball)))
        println("Size : $(round(last(r.chebyball),sigdigits=1))")
    end
    println("-------------------------------------")
    println("History:")
    print_ASs(r.ASs)
    println("=====================================")
end
function Base.:display(rs::Vector{Region})
    N = length(rs)
    max_iter, max_id = findmax(r.iter for r in rs)
    n_optimal = sum(r.state==OPTIMAL for r in rs)
    println("=====================================")
    println("# regions: $N")
    println("Max iter : $(max_iter)")
    println("Optimal  : $(Int(round(100*n_optimal/N)))%")
    println("-------------------------------------")
    println("Maximum sequence:")
    print_ASs(rs[max_id].ASs)
    println("=====================================")
end
