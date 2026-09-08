# Fluid-side cubic reconstruction avoids differentiating the contact-angle
# extension twice across a curved wall. This is a documented implementation
# choice; it is not specified in the target paper's appendix.
struct BoundaryFit
    i::Int
    j::Int
    samples::Vector{CartesianIndex{2}}
    coefficients::Matrix{Float64}
end
function boundary_fits(g::CutGeometry;radius=5)
    h=min(g.dx,g.dy);fits=BoundaryFit[]
    for j in 1:g.ny,i in 1:g.nx
        -2h<g.distance[i,j]<2h || continue
        samples=CartesianIndex{2}[];rows=Vector{Float64}[]
        for jj in max(1,j-radius):min(g.ny,j+radius),ii in max(1,i-radius):min(g.nx,i+radius)
            g.distance[ii,jj]>0 || continue
            x=(ii-i)*g.dx/h;y=(jj-j)*g.dy/h
            push!(samples,CartesianIndex(ii,jj));push!(rows,[1.,x,y,x*x/2,x*y,y*y/2,x^3/6,x*x*y/2,x*y*y/2,y^3/6])
        end
        length(rows)>=12 || error("Insufficient fluid samples for boundary curvature")
        a=reduce(vcat,transpose.(rows));s=svdvals(a)
        s[end]>1e-8s[1] || error("Ill-conditioned fluid-side boundary stencil")
        push!(fits,BoundaryFit(i,j,samples,pinv(a)))
    end
    fits
end
function boundary_curvature!(k,psi,g,fits)
    h=min(g.dx,g.dy)
    for fit in fits
        values=[psi[p] for p in fit.samples];a=fit.coefficients*values
        gx,gy=a[2],a[3];r2=gx*gx+gy*gy
        k[fit.i,fit.j]=r2>1e-40 ? -(a[4]*gy^2-2a[5]*gx*gy+a[6]*gx^2)/(h*r2^1.5) : 0.0
    end
    k
end
