struct WallArc
    i::Int
    j::Int
    center_x::Float64
    center_y::Float64
    radius::Float64
    angle0::Float64
    angle1::Float64
end
function wall_arcs(g::CutGeometry,circles)
    arcs=WallArc[]
    for j in 1:g.ny,i in 1:g.nx
        0<g.volume[i,j]<1 || continue
        x0,x1=(i-1)*g.dx,i*g.dx;y0,y1=(j-1)*g.dy,j*g.dy
        for c in circles
            hypot(clamp(c.x,x0,x1)-c.x,clamp(c.y,y0,y1)-c.y)>=c.radius && continue
            angles=[0.,2pi]
            for x in (x0,x1)
                z=(x-c.x)/c.radius
                if abs(z)<1
                    a=acos(z);push!(angles,a,2pi-a)
                end
            end
            for y in (y0,y1)
                z=(y-c.y)/c.radius
                if abs(z)<1
                    a=asin(z);push!(angles,mod(a,2pi),mod(pi-a,2pi))
                end
            end
            sort!(angles)
            for k in 1:length(angles)-1
                a,b=angles[k],angles[k+1];t=(a+b)/2
                x=c.x+c.radius*cos(t);y=c.y+c.radius*sin(t)
                if x0<=x<=x1 && y0<=y<=y1 && b>a
                    push!(arcs,WallArc(i,j,c.x,c.y,c.radius,a,b))
                end
            end
        end
    end
    arcs
end
# Integral of fluid-outward normal over the solid arc (opposite circle normal).
wall_normal_integral(a::WallArc)=(-a.radius*(sin(a.angle1)-sin(a.angle0)),
                                 a.radius*(cos(a.angle1)-cos(a.angle0)))
function wall_viscous_force(u,v,mu,g,arcs)
    fx=zeros(size(u));fy=zeros(size(v));probe=1.5max(g.dx,g.dy)
    for a in arcs
        theta=(a.angle0+a.angle1)/2;nx,ny=cos(theta),sin(theta)
        x=a.center_x+a.radius*nx;y=a.center_y+a.radius*ny
        # Two normal probes provide a second-order one-sided wall derivative.
        u1=interpolate(u,x+probe*nx,y+probe*ny,g.dx,g.dy)
        v1=interpolate(v,x+probe*nx,y+probe*ny,g.dx,g.dy)
        u2=interpolate(u,x+2probe*nx,y+2probe*ny,g.dx,g.dy)
        v2=interpolate(v,x+2probe*nx,y+2probe*ny,g.dx,g.dy)
        du=(4u1-u2)/(2probe);dv=(4v1-v2)/(2probe)
        muv=interpolate(mu,x+probe*nx,y+probe*ny,g.dx,g.dy)
        length=a.radius*(a.angle1-a.angle0);normal=du*nx+dv*ny
        fx[a.i,a.j]-=muv*length*(du+normal*nx)
        fy[a.i,a.j]-=muv*length*(dv+normal*ny)
    end
    fx,fy
end
