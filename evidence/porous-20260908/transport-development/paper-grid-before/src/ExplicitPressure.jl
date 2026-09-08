# Explicit pressure relaxation on merged cut volumes. Density/volume scaling
# replaces the Cartesian pointwise divergence. Outlet pressure remains fixed;
# subtracting a pressure mean here would change the outlet boundary condition.
function explicit_pressure!(p,b,a,g,cv,rho,dt;iterations=2,sound_cfl=.4)
    iterations>=1 && dt>0 && sound_cfl>0 || error("Invalid EPP parameters")
    density=volume_integrals(rho,g,cv)./cv.volumes
    d=density./cv.volumes
    cs=sound_cfl*min(g.dx,g.dy)/dt
    # D*A is similar to SPD sqrt(D)*A*sqrt(D). Its eigenvalues lie in
    # (0, 2max(D_ii*A_ii)] by the row Gershgorin bound. Require relaxation
    # eigenvalues <=1.8, avoiding unstable reuse of a Cartesian sound CFL.
    bound=2maximum(d.*a.diagonal)
    scale=min(1.,sqrt(1.8/(dt^2*cs^2*bound)))
    cs*=scale;alpha=dt^2*cs^2
    ap=similar(p)
    for _ in 1:iterations
        graph_apply!(ap,p,a)
        p.+=alpha.*d.*(b.-ap)
    end
    graph_apply!(ap,p,a)
    (;residual=norm(b-ap),sound_speed=cs,sound_speed_scale=scale,
      relaxation_spectral_upper_bound=alpha*bound)
end

"""Weakly-compressible pressure correction of supplied provisional faces.
This leaves a finite divergence, as intended for EPP. It is a pressure component,
not an accumulated-subcycling two-phase timestep implementation.
"""
function project_cut_weak!(uf,vf,p,dt,g,cv,rho;inlet_velocity=.003,iterations=2,sound_cfl=.4)
    uf[1,:].=inlet_velocity;vf[:,1].=0;vf[:,end].=0
    uf[g.aperture_x.==0].=0;vf[g.aperture_y.==0].=0
    a=pressure_graph(g,cv,rho);b=-integrated_divergence(uf,vf,g,cv)./dt
    report=explicit_pressure!(p,b,a,g,cv,rho,dt;iterations=iterations,sound_cfl=sound_cfl)
    for k in eachindex(a.left)
        i,j=a.left[k],a.right[k];dp=(j==0 ? 0.0 : p[j])-p[i]
        if a.face_axis[k]==1
            uf[a.face_i[k],a.face_j[k]]-=dt*a.mobility[k]*dp
        else
            vf[a.face_i[k],a.face_j[k]]-=dt*a.mobility[k]*dp
        end
    end
    report
end
