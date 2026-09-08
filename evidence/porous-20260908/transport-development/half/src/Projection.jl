# Conservative agglomerated control volumes for tiny cut cells. This is an
# explicit implementation choice, not a claim to duplicate Meyer's mixing.
struct ControlVolumes
    owner::Matrix{Int}
    volumes::Vector{Float64}
    representative_x::Vector{Float64}
    representative_y::Vector{Float64}
end
function control_volumes(g::CutGeometry;minimum_fraction=.5)
    0<minimum_fraction<=1 || error("Invalid agglomeration threshold")
    nx,ny=g.nx,g.ny;n=nx*ny;parent=collect(1:n)
    function root(k)
        while parent[k]!=k
            parent[k]=parent[parent[k]];k=parent[k]
        end
        k
    end
    small=sort(findall(v->0<v<minimum_fraction,vec(g.volume));by=k->g.volume[k])
    for k in small
        i=mod1(k,nx);j=(k-1)÷nx+1;best=0;score=-Inf
        for (ii,jj,ap) in ((i-1,j,g.aperture_x[i,j]),(i+1,j,g.aperture_x[i+1,j]),
                           (i,j-1,g.aperture_y[i,j]),(i,j+1,g.aperture_y[i,j+1]))
            1<=ii<=nx && 1<=jj<=ny || continue
            kk=ii+(jj-1)*nx
            if ap>1e-14 && g.volume[kk]>0 && root(kk)!=root(k)
                value=g.volume[kk]*ap
                if value>score;score=value;best=kk;end
            end
        end
        best>0 || error("Isolated tiny fluid control volume at ($i,$j), fraction=$(g.volume[i,j])")
        parent[root(k)]=root(best)
    end
    owner=zeros(Int,nx,ny);ids=Dict{Int,Int}();vol=Float64[];xx=Float64[];yy=Float64[]
    for j in 1:ny,i in 1:nx
        g.volume[i,j]==0 && continue
        k=root(i+(j-1)*nx)
        if !haskey(ids,k)
            ids[k]=length(vol)+1;push!(vol,0);push!(xx,0);push!(yy,0)
        end
        a=ids[k];owner[i,j]=a;v=g.volume[i,j]*g.dx*g.dy
        vol[a]+=v;xx[a]+=(i-.5)*g.dx*v;yy[a]+=(j-.5)*g.dy*v
    end
    xx./=vol;yy./=vol
    minimum(vol)>=minimum_fraction*g.dx*g.dy*(1-1e-10) || error("Agglomeration left a tiny control volume")
    ControlVolumes(owner,vol,xx,yy)
end

struct PressureGraph
    left::Vector{Int}
    right::Vector{Int}
    conductance::Vector{Float64}
    face_axis::Vector{Int}
    face_i::Vector{Int}
    face_j::Vector{Int}
    mobility::Vector{Float64}
    diagonal::Vector{Float64}
end
function pressure_graph(g,cv,rho)
    size(rho)==size(g.volume) && minimum(rho)>0 || error("Invalid fluid density")
    left=Int[];right=Int[];coef=Float64[];axis=Int[];ii=Int[];jj=Int[];mob=Float64[]
    function edge(a,b,ap,rf,spacing,length,ax,i,j)
        (a==0 && b==0 || a==b || ap==0) && return
        # A zero owner on an internal face is solid, never an implicit outlet.
        a>0 || error("Unexpected solid-left open face")
        push!(left,a);push!(right,b);push!(coef,ap*length/(rf*spacing))
        push!(axis,ax);push!(ii,i);push!(jj,j);push!(mob,1/(rf*spacing))
    end
    for j in 1:g.ny,i in 2:g.nx
        a,b=cv.owner[i-1,j],cv.owner[i,j]
        (a==0 || b==0) && continue
        edge(a,b,g.aperture_x[i,j],(rho[i-1,j]+rho[i,j])/2,g.dx,g.dy,1,i,j)
    end
    for j in 2:g.ny,i in 1:g.nx
        a,b=cv.owner[i,j-1],cv.owner[i,j]
        (a==0 || b==0) && continue
        edge(a,b,g.aperture_y[i,j],(rho[i,j-1]+rho[i,j])/2,g.dy,g.dx,2,i,j)
    end
    for j in 1:g.ny
        a=cv.owner[end,j];a==0 && continue
        edge(a,0,g.aperture_x[end,j],rho[end,j],g.dx/2,g.dy,1,g.nx+1,j)
    end
    diag=zeros(length(cv.volumes))
    for (a,b,w) in zip(left,right,coef)
        diag[a]+=w;b>0 && (diag[b]+=w)
    end
    minimum(diag)>0 || error("Unconnected pressure control volume")
    PressureGraph(left,right,coef,axis,ii,jj,mob,diag)
end
function graph_apply!(out,p,a::PressureGraph)
    fill!(out,0)
    @inbounds for k in eachindex(a.left)
        i,j=a.left[k],a.right[k];v=a.conductance[k]*(p[i]-(j==0 ? 0.0 : p[j]))
        out[i]+=v;j>0 && (out[j]-=v)
    end
    out
end
function graph_pcg!(p,b,a;rtol=1e-10,maxiter=20000,mg=nothing)
    ap=similar(p);graph_apply!(ap,p,a);r=b-ap;tol=max(1e-14,rtol*norm(b))
    norm(r)<=tol && return 0,norm(r)
    z=r./a.diagonal;mg!==nothing && graph_precondition!(z,r,mg);d=copy(z);rz=dot(r,z)
    for k in 1:maxiter
        graph_apply!(ap,d,a);den=dot(d,ap);den>0 || error("Pressure graph lost positive definiteness")
        alpha=rz/den;p.+=alpha.*d;r.-=alpha.*ap
        if norm(r)<=tol || k%50==0
            graph_apply!(ap,p,a);r.=b.-ap
            norm(r)<=tol && return k,norm(r)
        end
        z.=r./a.diagonal;mg!==nothing && graph_precondition!(z,r,mg);next=dot(r,z);d.=z.+(next/rz).*d;rz=next
    end
    error("Cut-cell pressure solve did not converge: residual=$(norm(r)), tolerance=$tol")
end
function integrated_divergence(uf,vf,g,cv)
    d=zeros(length(cv.volumes))
    for j in 1:g.ny,i in 1:g.nx
        a=cv.owner[i,j];a==0 && continue
        d[a]+=(g.aperture_x[i+1,j]*uf[i+1,j]-g.aperture_x[i,j]*uf[i,j])*g.dy+
              (g.aperture_y[i,j+1]*vf[i,j+1]-g.aperture_y[i,j]*vf[i,j])*g.dx
    end
    d
end
function project_cut!(uf,vf,p,dt,g,cv,rho;inlet_velocity=.003,rtol=1e-10,backend=:jacobi)
    dt>0 || error("Positive projection timestep required")
    uf[1,:].=inlet_velocity;vf[:,1].=0;vf[:,end].=0
    uf[g.aperture_x.==0].=0;vf[g.aperture_y.==0].=0
    a=pressure_graph(g,cv,rho);b=-integrated_divergence(uf,vf,g,cv)./dt
    backend in (:jacobi,:multigrid) || error("Unknown pressure backend")
    mg=backend==:multigrid ? graph_mg(a,cv,g) : nothing
    its,res=graph_pcg!(p,b,a;rtol=rtol,mg=mg)
    for k in eachindex(a.left)
        i,j=a.left[k],a.right[k];dp=(j==0 ? 0.0 : p[j])-p[i]
        if a.face_axis[k]==1
            uf[a.face_i[k],a.face_j[k]]-=dt*a.mobility[k]*dp
        else
            vf[a.face_i[k],a.face_j[k]]-=dt*a.mobility[k]*dp
        end
    end
    its,res
end
# Mix a conservative cell quantity within each agglomerated volume. Both sum(Vq)
# and constants are preserved; no threshold clipping is applied to q itself.
function agglomerate!(q,g,cv)
    mass=zeros(length(cv.volumes))
    for k in eachindex(q)
        a=cv.owner[k];a==0 && continue
        mass[a]+=g.volume[k]*g.dx*g.dy*q[k]
    end
    mass./=cv.volumes
    for k in eachindex(q)
        a=cv.owner[k];a==0 && continue;q[k]=mass[a]
    end
    q
end
