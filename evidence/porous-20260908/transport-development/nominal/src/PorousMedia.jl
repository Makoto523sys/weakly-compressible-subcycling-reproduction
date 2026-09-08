module PorousMedia
using LinearAlgebra
export Circle, CutGeometry, paper_circles, cut_geometry, cut_divergence,
       extend_contact!, impose_ghost_velocity!, liquid_contact_angle

struct Circle
    x::Float64
    y::Float64
    radius::Float64
end
struct CutGeometry
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
    volume::Matrix{Float64}
    aperture_x::Matrix{Float64}
    aperture_y::Matrix{Float64}
    distance::Matrix{Float64}
    normal_x::Matrix{Float64}
    normal_y::Matrix{Float64}
end
function paper_circles()
    shifts=Dict((0,2)=>.00002,(1,1)=>.00002,(2,2)=>.00002,(3,1)=>-.00002,(4,1)=>-.00002)
    [Circle(.0006+sqrt(3)/2*.001*k,
            .001*j+(isodd(k) ? .0005 : 0)+get(shifts,(k,j),0.0),.0004)
     for k in 0:4 for j in (iseven(k) ? (0:3) : (0:2))]
end
cross2(a,b)=a[1]*b[2]-a[2]*b[1]
dot2(a,b)=a[1]*b[1]+a[2]*b[2]
# Exact line/sector integration of a polygon edge intersected with a disk.
function edge_disk(a,b,r)
    d=(b[1]-a[1],b[2]-a[2]); aa=dot2(d,d);bb=2dot2(a,d);cc=dot2(a,a)-r*r
    times=[zero(r),one(r)];disc=bb*bb-4aa*cc
    if disc>0
        for t in ((-bb-sqrt(disc))/(2aa),(-bb+sqrt(disc))/(2aa))
            0<t<1 && push!(times,t)
        end
    end
    sort!(times);area=zero(r)
    for k in 1:length(times)-1
        t0,t1=times[k],times[k+1]
        p=(a[1]+t0*d[1],a[2]+t0*d[2]);q=(a[1]+t1*d[1],a[2]+t1*d[2])
        mid=((p[1]+q[1])/2,(p[2]+q[2])/2)
        area+=dot2(mid,mid)<r*r ? cross2(p,q)/2 : r*r*atan(cross2(p,q),dot2(p,q))/2
    end
    area
end
function precise_fluid_fraction(x0,x1,y0,y1,circles)
    # Near a disk corner tangency, subtracting two O(cell area) Float64
    # quantities can create a disconnected phantom fluid sliver. Recompute
    # that cancellation at 256-bit precision, rather than raising the cutoff.
    setprecision(BigFloat,256) do
        a,b,c,d=BigFloat.((x0,x1,y0,y1));area=(b-a)*(d-c);covered=zero(area)
        for circle in circles
            x,y,r=BigFloat.((circle.x,circle.y,circle.radius))
            hypot(clamp(x,a,b)-x,clamp(y,c,d)-y)>=r && continue
            p=((a-x,c-y),(b-x,c-y),(b-x,d-y),(a-x,d-y))
            covered+=sum(edge_disk(p[k],p[mod1(k+1,4)],r) for k in 1:4)
        end
        Float64(clamp((area-covered)/area,0,1))
    end
end
function disk_rectangle(x0,x1,y0,y1,c::Circle)
    # Fast exclusion/inclusion also avoids cancellation far from the circle.
    near=hypot(clamp(c.x,x0,x1)-c.x,clamp(c.y,y0,y1)-c.y)
    near>=c.radius && return 0.0
    far=hypot(max(abs(x0-c.x),abs(x1-c.x)),max(abs(y0-c.y),abs(y1-c.y)))
    far<=c.radius && return (x1-x0)*(y1-y0)
    p=((x0-c.x,y0-c.y),(x1-c.x,y0-c.y),(x1-c.x,y1-c.y),(x0-c.x,y1-c.y))
    clamp(sum(edge_disk(p[k],p[mod1(k+1,4)],c.radius) for k in 1:4),0.0,(x1-x0)*(y1-y0))
end
function blocked_length(fixed,lo,hi,c::Circle,axis)
    center=axis==1 ? c.x : c.y;along=axis==1 ? c.y : c.x
    q=c.radius^2-(fixed-center)^2;q<=0 && return 0.0
    a=sqrt(q);max(0.0,min(hi,along+a)-max(lo,along-a))
end
function cut_geometry(nx,ny,lx,ly,circles)
    nx>=4 && ny>=4 && min(lx,ly)>0 || error("Invalid grid")
    for a in 1:length(circles),b in a+1:length(circles)
        c,d=circles[a],circles[b]
        hypot(c.x-d.x,c.y-d.y)>=c.radius+d.radius || error("Overlapping circles require a union-area algorithm")
    end
    dx=lx/nx;dy=ly/ny;vol=ones(nx,ny);ax=ones(nx+1,ny);ay=ones(nx,ny+1)
    dist=fill(Inf,nx,ny);normalx=zeros(nx,ny);normaly=zeros(nx,ny)
    for j in 1:ny,i in 1:nx
        x=(i-.5)*dx;y=(j-.5)*dy
        for c in circles
            r=hypot(x-c.x,y-c.y);d=r-c.radius
            if d<dist[i,j]
                dist[i,j]=d
                normalx[i,j]=r==0 ? 1.0 : (x-c.x)/r
                normaly[i,j]=r==0 ? 0.0 : (y-c.y)/r
            end
            # Normalize with the same endpoint differences used in the intersection.
            # Using dx*dy instead can leave phantom fluid in a fully solid cell
            # when i*dx-(i-1)*dx rounds differently from dx.
            cellarea=(i*dx-(i-1)*dx)*(j*dy-(j-1)*dy)
            vol[i,j]-=disk_rectangle((i-1)*dx,i*dx,(j-1)*dy,j*dy,c)/cellarea
        end
        vol[i,j]=clamp(vol[i,j],0.0,1.0)
        if 0<vol[i,j]<1e-10
            vol[i,j]=precise_fluid_fraction((i-1)*dx,i*dx,(j-1)*dy,j*dy,circles)
        end
        vol[i,j]<64eps(Float64) && (vol[i,j]=0.0)
    end
    for j in 1:ny,i in 1:nx+1
        ax[i,j]=clamp(1-sum(blocked_length((i-1)*dx,(j-1)*dy,j*dy,c,1) for c in circles;init=0.0)/(j*dy-(j-1)*dy),0.0,1.0)
    end
    for j in 1:ny+1,i in 1:nx
        ay[i,j]=clamp(1-sum(blocked_length((j-1)*dy,(i-1)*dx,i*dx,c,2) for c in circles;init=0.0)/(i*dx-(i-1)*dx),0.0,1.0)
    end
    CutGeometry(nx,ny,dx,dy,vol,ax,ay,dist,normalx,normaly)
end
# Zero normal flux on the circular solid surfaces is implicit. Outer face fluxes
# remain explicit; the caller chooses inlet/outlet/wall boundary conditions.
function cut_divergence(fx,fy,g::CutGeometry)
    d=zeros(g.nx,g.ny)
    for j in 1:g.ny,i in 1:g.nx
        v=g.volume[i,j];v==0 && continue
        d[i,j]=((g.aperture_x[i+1,j]*fx[i+1,j]-g.aperture_x[i,j]*fx[i,j])/g.dx+
                (g.aperture_y[i,j+1]*fy[i,j+1]-g.aperture_y[i,j]*fy[i,j])/g.dy)/v
    end
    d
end

# Bilinear cell-centre interpolation with constant extension at domain edges.
function stencil(x,y,nx,ny,dx,dy)
    xx=clamp(x/dx+.5,1.0,Float64(nx));yy=clamp(y/dy+.5,1.0,Float64(ny))
    i=min(floor(Int,xx),nx-1);j=min(floor(Int,yy),ny-1);a=xx-i;b=yy-j
    ((i,j,(1-a)*(1-b)),(i+1,j,a*(1-b)),(i,j+1,(1-a)*b),(i+1,j+1,a*b))
end
function interpolate(q,x,y,dx,dy)
    points=stencil(x,y,size(q)...,dx,dy)
    i,j,_=points[1];base=q[i,j]
    # Difference form preserves constants exactly, including inside repeated
    # contact extensions where tiny errors otherwise acquire a unit normal.
    base+sum(w*(q[ii,jj]-base) for (ii,jj,w) in points)
end

"""Extrapolate only solid values along contact-angle characteristics.
Normals point from solid into fluid. phi=1 is liquid. This is a boundary
component, not a validated moving-contact-line CFD solver.
"""
function extend_contact!(phi,distance,nxwall,nywall,dx,dy,theta;iterations=80,band=3max(dx,dy))
    0<theta<180 || error("Liquid angle must be strictly between 0 and 180 degrees")
    nx,ny=size(phi);cotangent=cotd(theta);dtau=.5min(dx,dy);old=similar(phi)
    for _ in 1:iterations
        old.=phi
        for j in 1:ny,i in 1:nx
            -band<distance[i,j]<0 || continue
            gx=(old[min(i+1,nx),j]-old[max(i-1,1),j])/(2dx)
            gy=(old[i,min(j+1,ny)]-old[i,max(j-1,1)])/(2dy)
            wx,wy=nxwall[i,j],nywall[i,j]
            dn=gx*wx+gy*wy;tx=gx-dn*wx;ty=gy-dn*wy;tn=hypot(tx,ty)
            tx=tn>1e-30 ? tx/tn : 0.0;ty=tn>1e-30 ? ty/tn : 0.0
            # u_ext is tangent to the prescribed interface and points into solid.
            ux=-wx-cotangent*tx;uy=-wy-cotangent*ty;un=hypot(ux,uy)
            x=(i-.5)*dx-dtau*ux/un;y=(j-.5)*dy-dtau*uy/un
            phi[i,j]=interpolate(old,x,y,dx,dy)
        end
    end
    phi
end
function liquid_contact_angle(gx,gy,wx,wy)
    r=hypot(gx,gy);r>0 || error("Undefined phase normal")
    acosd(clamp(-(gx*wx+gy*wy)/r,-1,1))
end

"""Stationary sharp-wall image/ghost velocity relation, solved iteratively.
The image point is the reflection across the nearest circular surface. Its own
interpolation weight is eliminated algebraically; adjacent ghosts are iterated.
"""
function impose_ghost_velocity!(u,v,g::CutGeometry;iterations=100,tolerance=1e-12,band=3max(g.dx,g.dy))
    ghosts=[(i,j) for j in 1:g.ny for i in 1:g.nx if -band<g.distance[i,j]<0]
    probes=[stencil((i-.5)*g.dx-2g.distance[i,j]*g.normal_x[i,j],
                    (j-.5)*g.dy-2g.distance[i,j]*g.normal_y[i,j],g.nx,g.ny,g.dx,g.dy) for (i,j) in ghosts]
    for iteration in 1:iterations
        change=0.0
        for ((i,j),p) in zip(ghosts,probes)
            own=sum(w for (ii,jj,w) in p if ii==i && jj==j;init=0.0)
            a=-sum(w*u[ii,jj] for (ii,jj,w) in p if ii!=i || jj!=j;init=0.0)/(1+own)
            b=-sum(w*v[ii,jj] for (ii,jj,w) in p if ii!=i || jj!=j;init=0.0)/(1+own)
            change=max(change,abs(a-u[i,j]),abs(b-v[i,j]));u[i,j]=a;v[i,j]=b
        end
        change<=tolerance && return iteration,change
    end
    error("Ghost velocity iteration failed to converge")
end
include("Projection.jl")
include("GraphMultigrid.jl")
include("WallGeometry.jl")
include("BoundaryCurvature.jl")
include("Transport.jl")

end
