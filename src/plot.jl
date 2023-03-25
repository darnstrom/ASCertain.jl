## Print ASs in a readable way
function print_ASs(ASs::BitMatrix)
    AS = Int64[];
    inds = collect(1:size(ASs,1))
    printstyled("$(lpad("ini",3," ")) "; color = :yellow)
    for i in 1:size(ASs,2)-1
        println("$AS ")
        if(sum(ASs[:,i])<sum(ASs[:,i+1])) # Addition	
            add_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            push!(AS,add_ind[1])
            printstyled("+$(lpad(add_ind[1],2," ")) "; color = :green)
        else # Removal
            rm_ind = inds[.!(ASs[:,i].⊻ .!ASs[:,i+1])]
            deleteat!(AS,findfirst(AS.==rm_ind[1]))
            printstyled("-$(lpad(rm_ind[1],2," ")) "; color = :red)
        end
    end
    println("$AS")
end
## Plot partition
function pplot(rs::Vector{<:AbstractRegion};key=:iter, fix_ids = nothing, fix_vals=nothing,opts=Dict{Symbol,Any}(), clabel=nothing)
    isempty(rs) && error("Cannot plot empty collection")
    nth = size(rs[1].Ath,1)
    ids = isnothing(fix_ids) ? collect(3:nth) : fix_ids
    values = isnothing(fix_vals) ? zeros(nth-2) : fix_vals
    free_ids = setdiff(1:nth,ids)

    ps = PolyDAQP.Polyhedron[]
    cs = typeof(getfield(rs[1],key))[]
    for r in rs
        p = Polyhedron(slice(r.Ath,r.bth,ids;values)...)
        isempty(p) && continue
        push!(ps,minrep(p))
        push!(cs,getfield(r,key))
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
    else
        push!(lopts, :colorbar_style => "{xlabel = \\# of "*String(key)*"}")
    end
    opts = merge(lopts,opts)
    PolyDAQP.pplot(ps;cs,opts)
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
    println("=====================================")
    println("History:")
    print_ASs(r.ASs)
end
