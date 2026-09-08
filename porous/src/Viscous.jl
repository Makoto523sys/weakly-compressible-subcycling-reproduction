"""Integrated Newtonian viscous force on each raw fluid cut volume.
Includes open-face stresses and the solid-circle traction. Exterior top/bottom
are no slip, inlet velocity prescribed, outlet zero velocity gradient. Callers
must subsequently accumulate raw forces into merged volumes before advancing.
"""
function viscous_force(phi,u,v,g,arcs,c::TransportParameters)
    uu=copy(u);vv=copy(v);impose_ghost_velocity!(uu,vv,g)
    mu=@. c.mu_g+(c.mu_l-c.mu_g)*phi
    minimum(mu[g.volume.>0])>=0 || error("Negative fluid viscosity")
    # Extend material values to wall ghosts without changing conserved phase.
    q=copy(phi);extend_contact!(q,g.distance,g.normal_x,g.normal_y,g.dx,g.dy,c.theta)
    mu=@. c.mu_g+(c.mu_l-c.mu_g)*q
    ms(i,j)=mu[clamp(i,1,g.nx),clamp(j,1,g.ny)]
    us(i,j)=transport_sample(uu,i,j,:u,c);vs(i,j)=transport_sample(vv,i,j,:v,c)
    ux=zeros(g.nx+1,g.ny);vx=similar(ux);uy=zeros(g.nx,g.ny+1);vy=similar(uy)
    for j in 1:g.ny,i in 1:g.nx+1
        muf=(ms(i-1,j)+ms(i,j))/2
        ux[i,j]=2muf*(us(i,j)-us(i-1,j))/g.dx
        vx[i,j]=muf*((vs(i,j)-vs(i-1,j))/g.dx+
            (us(i,j+1)+us(i-1,j+1)-us(i,j-1)-us(i-1,j-1))/(4g.dy))
    end
    for j in 1:g.ny+1,i in 1:g.nx
        muf=(ms(i,j-1)+ms(i,j))/2
        vy[i,j]=2muf*(vs(i,j)-vs(i,j-1))/g.dy
        uy[i,j]=muf*((us(i,j)-us(i,j-1))/g.dy+
            (vs(i+1,j)+vs(i+1,j-1)-vs(i-1,j)-vs(i-1,j-1))/(4g.dx))
    end
    wallx,wally=wall_viscous_force(uu,vv,mu,g,arcs)
    fx=copy(wallx);fy=copy(wally)
    for j in 1:g.ny,i in 1:g.nx
        g.volume[i,j]>0 || continue
        fx[i,j]+=(g.aperture_x[i+1,j]*ux[i+1,j]-g.aperture_x[i,j]*ux[i,j])*g.dy+
            (g.aperture_y[i,j+1]*uy[i,j+1]-g.aperture_y[i,j]*uy[i,j])*g.dx
        fy[i,j]+=(g.aperture_x[i+1,j]*vx[i+1,j]-g.aperture_x[i,j]*vx[i,j])*g.dy+
            (g.aperture_y[i,j+1]*vy[i,j+1]-g.aperture_y[i,j]*vy[i,j])*g.dx
    end
    (;fx,fy,wallx,wally,ux,uy,vx,vy)
end
