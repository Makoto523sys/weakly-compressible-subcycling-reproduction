# Galerkin aggregation for an SPD pressure graph with a Dirichlet outlet.
# No geometry or fine-grid pressure coefficient is changed by this preconditioner.
using SparseArrays
struct GraphMGLevel
    a::SparseMatrixCSC{Float64,Int}
    prolongation::SparseMatrixCSC{Float64,Int}
    diagonal::Vector{Float64}
    x::Vector{Float64}
    b::Vector{Float64}
    residual::Vector{Float64}
end
function graph_matrix(g::PressureGraph)
    n=length(g.diagonal);ii=collect(1:n);jj=copy(ii);vv=copy(g.diagonal)
    for k in eachindex(g.left)
        a,b=g.left[k],g.right[k];b==0 && continue;w=g.conductance[k]
        push!(ii,a);push!(jj,b);push!(vv,-w)
        push!(ii,b);push!(jj,a);push!(vv,-w)
    end
    sparse(ii,jj,vv,n,n)
end
function graph_mg(a::PressureGraph,cv,g)
    matrix=graph_matrix(a);xx=copy(cv.representative_x);yy=copy(cv.representative_y)
    spacingx=g.dx;spacingy=g.dy;levels=GraphMGLevel[]
    while size(matrix,1)>128
        spacingx*=2;spacingy*=2;groups=Dict{Tuple{Int,Int},Int}();map=Int[];cx=Float64[];cy=Float64[]
        for (x,y) in zip(xx,yy)
            key=(floor(Int,x/spacingx),floor(Int,y/spacingy))
            if !haskey(groups,key)
                groups[key]=length(groups)+1;push!(cx,(key[1]+.5)*spacingx);push!(cy,(key[2]+.5)*spacingy)
            end
            push!(map,groups[key])
        end
        n=size(matrix,1);nc=length(groups);nc<n || error("Pressure aggregation did not coarsen")
        prolong=sparse(1:n,map,ones(n),n,nc)
        push!(levels,GraphMGLevel(matrix,prolong,diag(matrix),zeros(n),zeros(n),zeros(n)))
        matrix=transpose(prolong)*matrix*prolong;xx=cx;yy=cy
    end
    n=size(matrix,1)
    push!(levels,GraphMGLevel(matrix,spzeros(n,0),diag(matrix),zeros(n),zeros(n),zeros(n)))
    levels,cholesky(Symmetric(Matrix(matrix)))
end
function graph_cycle!(levels,coarse,k=1)
    l=levels[k]
    if k==length(levels)
        l.x.=coarse\l.b;return l.x
    end
    @. l.x=(2/3)*l.b/l.diagonal
    mul!(l.residual,l.a,l.x);l.residual.=l.b.-l.residual
    child=levels[k+1];mul!(child.b,transpose(l.prolongation),l.residual)
    graph_cycle!(levels,coarse,k+1)
    mul!(l.residual,l.prolongation,child.x);l.x.+=l.residual
    mul!(l.residual,l.a,l.x)
    @. l.x+=(2/3)*(l.b-l.residual)/l.diagonal
    l.x
end
function graph_precondition!(z,r,mg)
    levels,coarse=mg;levels[1].b.=r;z.=graph_cycle!(levels,coarse);z
end
