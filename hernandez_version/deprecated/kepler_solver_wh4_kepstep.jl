# Wisdom & Hernandez version of Kepler solver, but with quartic
# convergence.

function kep_eqn!(x0::Array{Float64,1},v0::Array{Float64,1},k::Float64,h::Float64,s0::Float64,ks::Kepler_step)
# Solves equation (35) from Wisdom & Hernandez

# First, compute beta:
ks.r0 = norm(x0)
ks.beta = 2.0*k/ks.r0-dot(v0,v0)
ks.dr0dt = dot(x0,v0)/ks.r0
# Now, solve for s in elliptical Kepler case:
if ks.beta > 1e-15
  if s0 == 0.0
# Initial guess:
    s0 = h/r0
  end
  ks.sqb = sqrt(beta)
  ks.y = 0.0; ks.yp = 1.0
  ks.iter = 0
  ks.ds = Inf
  while ks.iter == 0 || (abs(ks.ds) > 1e-8*s0 && ks.iter < 10)
    ks.xx = 0.5*ks.sqb*s0
    ks.sx = sin(ks.xx) 
    ks.cx = cos(ks.xx)
# Third derivative:
    ks.yppp = (k-ks.r0*ks.beta)*(ks.cx*ks.cx-ks.sx*ks.sx) - ks.sqb*ks.r0*ks.dr0dt*2.0*ks.*sx*ks.cx
# Take derivative:
    ks.yp = (-ks.yppp+ k)/ks.beta
# Second derivative:
    ks.ypp = -(ks.r0-k/ks.beta)*2.*ks.sx*ks.cx*ks.sqb  + ks.r0*ks.dr0dt*(ks.cx*ks.cx-ks.sx*ks.sx)
    ks.y  = (-ks.ypp +ks.r0*ks.dr0dt +k*s0)/ks.beta - h  # eqn 35
# Now, compute fourth-order estimate:
    ks.ds = -ks.y/(ks.yp-ks.y/ks.yp*ks.ypp*.5)
    ks.ds = -ks.y/(ks.yp+ks.ds*.5*(ks.ypp+ks.ds*ks.yppp/3.))
    s0 += ks.ds
    ks.iter +=1
  end
  ks.xx = 0.5*ks.sqb*s0; ks.sx = sin(ks.xx) ; ks.cx = cos(ks.xx)
# Now, compute final values:
  ks.g1bs = 2.0*ks.sx*ks.cx/ks.sqb
  ks.g2bs = 2.0*ks.sx^2/ks.beta
  ks.f = 1.0 - k/r0*g2bs # eqn (25)
  ks.g = r0*g1bs + r0*dr0dt*g2bs # eqn (27)
  for j=1:3
    ks.x[j] = x0[j]*ks.f+v0[j]*g
  end
  ks.r = norm(ks.x)
  ks.dfdt = -k/ks.r/ks.r0*ks.g1bs
  ks.dgdt = ks.r0/ks.r*(1.0-ks.beta*ks.g2bs+ks.dr0dt*ks.g1bs)
  for j=1:3
    ks.v[j] = x0[j]*ks.dfdt+v0[j]*ks.dgdt
  end
else
  if ks.beta < -1e-15
# Hyperbolic case: (see 8/9/17 notes)
    if s0 == 0.0
# Initial guess:
      s0 = h/ks.r0
    end
    ks.sqb = sqrt(-ks.beta)
    ks.y = 0.0; ks.yp = 1.0
    ks.iter = 0
    ks.ds = Inf
    while ks.iter == 0 || (abs(ks.ds) > 1e-8*s0 && ks.iter < 10)
      ks.xx = 0.5*ks.sqb*s0; ks.sx = sinh(ks.xx) ; ks.cx = cosh(ks.xx)
# Third derivative:
      ks.yppp = (k-ks.r0*ks.beta)*(ks.cx*ks.cx+ks.sx*ks.sx) + ks.sqb*ks.r0*ks.dr0dt*2.0*ks.sx*ks.cx
# Take derivative:
      ks.yp = (-ks.yppp+ k)/ks.beta
# Second derivative:
      ks.ypp = 2.0*(ks.r0-k/ks.beta)*ks.cx*ks.sx*ks.sqb  + ks.r0*ks.dr0dt*(ks.cx*ks.cx+ks.sx*ks.sx)
      ks.y  = (-ks.ypp +ks.r0*ks.dr0dt +k*s0)/ks.beta - h  # eqn 35
# Now, compute fourth-order estimate:
      ks.ds = -ks.y/(ks.yp-ks.y/ks.yp*ks.ypp*.5)
      ks.ds = -ks.y/(ks.yp+ks.ds*.5*(ks.ypp+ks.ds*ks.yppp/3.))
      s0 += ks.ds
      ks.iter +=1
    end
    ks.xx = 0.5*ks.sqb*s0; ks.sx = sinh(ks.xx) ; ks.cx = cosh(ks.xx)
# Now, compute final values:
    ks.g1bs = 2.0*ks.sx*ks.cx/ks.sqb
    ks.g2bs = -2.0*ks.sx^2/ks.beta
    ks.f = 1.0 - k/ks.r0*ks.g2bs # eqn (25)
    ks.g = ks.r0*ks.g1bs + ks.r0*ks.dr0dt*ks.g2bs # eqn (27)
    for j=1:3
      ks.x[j] = x0[j]*ks.f+v0[j]*ks.g
    end
    ks.r = norm(ks.x)
    ks.dfdt = -k/ks.r/ks.r0*ks.g1bs
    ks.dgdt = ks.r0/ks.r*(1.0-ks.beta*ks.g2bs+ks.dr0dt*ks.g1bs)
    for j=1:3
      ks.v[j] = x0[j]*ks.dfdt+v0[j]*ks.dgdt
    end
  else
# Parabolic case: cubic equation
    ks.Q3 = (ks.r0*ks.dr0dt/k)^2 - 2.0*ks.r0/k # = (r0/k)^2*(dr0dt^2-2k/r0) = (r0/k)^2*(dr0dt^2-v0^2) < 0
    ks.R3 = (ks.r0*ks.dr0dt/k)^3-3.0*ks.r0^2*ks.dr0dt/k^2-3.0*h/k # = (r0/k)^3*dr0dt*(dr0dt^2 - 3/2*v0^2 -3*h*(k/r0)^3/dr0dt/k)
    if ks.R3^2 >= ks.Q3^3  # This should *always* be the case since Q <= 0
#      println("One real root:  Q: ",Q3," R: ",R3)
      ks.A3 = -sign(ks.R3)*(abs(ks.R3)+sqrt(abs(ks.R3^2-ks.Q3^3)))^(1//3)
      if ks.A3 == 0.0
        ks.B3 = 0.0
      else
        ks.B3 = ks.Q3/ks.A3
      end
      s0 = ks.A3 + ks.B3 - ks.r0*ks.dr0dt/k
    else
      println("Error: Q should be negative")
    end
    ks.g1bs = s0
    ks.g2bs = 0.5*s0^2
    ks.f = 1.0 - k/ks.r0*ks.g2bs # eqn (25)
    ks.g = ks.r0*ks.g1bs + ks.r0*ks.dr0dt*ks.g2bs # eqn (27)
    for j=1:3
      ks.x[j] = x0[j]*ks.f+v0[j]*ks.g
    end
    ks.r = norm(x)
    ks.dfdt = -k/ks.r/ks.r0*ks.g1bs
    ks.dgdt = ks.r0/ks.r*(1.0-ks.beta*ks.g2bs+ks.dr0dt*ks.g1bs)
    for j=1:3
      ks.v[j] = x0[j]*ks.dfdt+v0[j]*ks.dgdt
    end
# recompute beta:
#    beta_final = 2.0*k/r-dot(v,v)
#    println("final beta: ",beta_final," rdot: ",dr0dt)
  end
end
#println("s0: ",s0)
return ks.x,ks.v,s0
end
