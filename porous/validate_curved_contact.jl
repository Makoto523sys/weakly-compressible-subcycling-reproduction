using LinearAlgebra
include("src/PorousMedia.jl")
using .PorousMedia
const P=PorousMedia
function main()
    theta=150.;rs=.23;rl=.12;sep=sqrt(rs^2+rl^2-2rs*rl*cosd(theta))
    sx=.4;sy=.5;cx=sx+sep;cy=sy
    xcontact=(rs^2-rl^2+sep^2)/(2sep)
    yc=sqrt(rs^2-xcontact^2)
    println("nx,contact_angle_degrees,central_curvature,fluid_side_curvature,expected_curvature")
    for n in (100,200,400)
        h=1/n;g=cut_geometry(n,n,1.,1.,[Circle(sx,sy,rs)])
        phi=[.5*(1+tanh((rl-hypot((i-.5)*h-cx,(j-.5)*h-cy))/(2h))) for i in 1:n,j in 1:n]
        # Remove analytic solid extension, so the boundary algorithm must supply it.
        for j in 1:n,i in 1:n
            if g.distance[i,j]<0
                x=(i-.5)*h-g.distance[i,j]*g.normal_x[i,j]
                y=(j-.5)*h-g.distance[i,j]*g.normal_y[i,j]
                phi[i,j]=.5*(1+tanh((rl-hypot(x-cx,y-cy))/(2h)))
            end
        end
        extend_contact!(phi,g.distance,g.normal_x,g.normal_y,h,h,theta;iterations=100,band=5h)
        psi=@. h*log((clamp(phi,0,1)+1e-100)/(1-clamp(phi,0,1)+1e-100))
        gx=zeros(n,n);gy=zeros(n,n)
        for j in 2:n-1,i in 2:n-1
            gx[i,j]=(psi[i+1,j]-psi[i-1,j])/(2h);gy[i,j]=(psi[i,j+1]-psi[i,j-1])/(2h)
        end
        normg=hypot.(gx,gy).+1e-100;nx=gx./normg;ny=gy./normg;k=zeros(n,n)
        for j in 2:n-1,i in 2:n-1
            k[i,j]=-(nx[i+1,j]-nx[i-1,j]+ny[i,j+1]-ny[i,j-1])/(2h)
        end
        x=sx+xcontact;y=sy+yc;wx=xcontact/rs;wy=yc/rs
        angle=liquid_contact_angle(P.interpolate(gx,x,y,h,h),P.interpolate(gy,x,y,h,h),wx,wy)
        curvature=P.interpolate(k,x,y,h,h)
        P.boundary_curvature!(k,psi,g,P.boundary_fits(g))
        fitted=P.interpolate(k,x,y,h,h)
        println(join((n,angle,curvature,fitted,1/rl),','));flush(stdout)
    end
end
main()
