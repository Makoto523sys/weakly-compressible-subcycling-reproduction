using LinearAlgebra
include("src/PorousMedia.jl")
using .PorousMedia
const P=PorousMedia
println("nx,torque,expected_torque,relative_error")
for n in (40,80,160)
    c=Circle(.5,.5,.21);g=cut_geometry(n,n,1.,1.,[c]);arcs=P.wall_arcs(g,[c]);h=1/n
    u=zeros(n,n);v=zeros(n,n)
    for j in 1:n,i in 1:n
        x=(i-.5)*h-c.x;y=(j-.5)*h-c.y;r2=x*x+y*y
        if r2>c.radius^2
            u[i,j]=-(1-c.radius^2/r2)*y;v[i,j]=(1-c.radius^2/r2)*x
        end
    end
    fx,fy=P.wall_viscous_force(u,v,ones(n,n),g,arcs)
    torque=sum(((i-.5)*h-c.x)*fy[i,j]-((j-.5)*h-c.y)*fx[i,j] for i in 1:n,j in 1:n)
    expected=-4pi*c.radius^2
    println(join((n,torque,expected,abs(torque-expected)/abs(expected)),','))
end
