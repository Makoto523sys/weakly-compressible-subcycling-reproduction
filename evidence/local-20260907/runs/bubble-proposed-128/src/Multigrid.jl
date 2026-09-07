# Symmetric aggregation V-cycle used only as a PCG preconditioner. The solved
# fine-grid operator and true-residual tolerance remain unchanged.
struct MGLevel
    ax::Matrix{Float64}
    ay::Matrix{Float64}
    diag::Matrix{Float64}
    x::Matrix{Float64}
    b::Matrix{Float64}
    tmp::Matrix{Float64}
end
function mg_level(ax,ay)
    nx=size(ay,1);ny=size(ax,2)
    d=[ax[i,j]+ax[i+1,j]+ay[i,j]+ay[i,j+1] for i in 1:nx,j in 1:ny]
    MGLevel(ax,ay,d,zeros(nx,ny),zeros(nx,ny),zeros(nx,ny))
end
function mg_apply!(out,x,l)
    nx,ny=size(x)
    @inbounds for j in 1:ny,i in 1:nx
        v=l.diag[i,j]*x[i,j]
        i>1 && (v-=l.ax[i,j]*x[i-1,j])
        i<nx && (v-=l.ax[i+1,j]*x[i+1,j])
        j>1 && (v-=l.ay[i,j]*x[i,j-1])
        j<ny && (v-=l.ay[i,j+1]*x[i,j+1])
        out[i,j]=v
    end
    out
end
function mg_setup(rx,ry,c)
    ax=1.0./(rx.*dx(c)^2);ay=1.0./(ry.*dy(c)^2)
    ax[1,:].=0;ax[end,:].=0;ay[:,1].=0;ay[:,end].=0
    levels=[mg_level(ax,ay)]
    nx,ny=c.nx,c.ny
    while iseven(nx) && iseven(ny) && min(nx,ny)>4
        nx÷=2;ny÷=2
        fine=levels[end]
        ax=[fine.ax[2i-1,2j-1]+fine.ax[2i-1,2j] for i in 1:nx+1,j in 1:ny]
        ay=[fine.ay[2i-1,2j-1]+fine.ay[2i,2j-1] for i in 1:nx,j in 1:ny+1]
        push!(levels,mg_level(ax,ay))
    end
    l=levels[end];n=length(l.x);a=zeros(n,n)
    for k in 1:n
        fill!(l.x,0);l.x[k]=1
        mg_apply!(l.tmp,l.x,l);a[:,k].=vec(l.tmp)
    end
    # Remove the constant nullspace with a rank-one term at the coarse scale.
    scale=maximum(l.diag)
    a .+= scale/n
    levels,cholesky(Symmetric(a))
end
function mg_cycle!(levels,coarse,k=1)
    l=levels[k];fill!(l.x,0)
    if k==length(levels)
        l.x .= reshape(coarse\vec(l.b),size(l.x))
        zero_mean!(l.x)
        return l.x
    end
    for _ in 1:1
        mg_apply!(l.tmp,l.x,l)
        @. l.x += (2/3)*(l.b-l.tmp)/l.diag
    end
    mg_apply!(l.tmp,l.x,l)
    child=levels[k+1];nx,ny=size(child.b)
    @inbounds for j in 1:ny,i in 1:nx
        v=0.0
        for jj in 2j-1:2j,ii in 2i-1:2i
            v+=l.b[ii,jj]-l.tmp[ii,jj]
        end
        child.b[i,j]=v
    end
    zero_mean!(child.b);mg_cycle!(levels,coarse,k+1)
    @inbounds for j in 1:ny,i in 1:nx
        for jj in 2j-1:2j,ii in 2i-1:2i
            l.x[ii,jj]+=child.x[i,j]
        end
    end
    for _ in 1:1
        mg_apply!(l.tmp,l.x,l)
        @. l.x += (2/3)*(l.b-l.tmp)/l.diag
    end
    zero_mean!(l.x)
end
function mg_precondition!(z,r,mg)
    levels,coarse=mg
    levels[1].b .= r
    z .= mg_cycle!(levels,coarse)
    z
end
