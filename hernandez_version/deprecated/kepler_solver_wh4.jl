# Wisdom & Hernandez version of Kepler solver, but with quartic
# convergence.

function kep_eqn!(x0::Array{Float64,1},v0::Array{Float64,1},k::Float64,h::Float64,s0::Float64)
# Solves equation (35) from Wisdom & Hernandez

# First, compute beta:
r0 = norm(x0)
beta = 2.0*k/r0-dot(v0,v0)
dr0dt = dot(x0,v0)/r0
x = zeros(3)
v = zeros(3)
#println("r0: ",r0," dr0dt: ",dr0dt)
# Now, solve for s in elliptical Kepler case:
if beta > 1e-15
  if s0 == 0.0
# Initial guess:
    s0 = h/r0
  end
  sqb = sqrt(beta)
# h = r0 * 2.*s2*c2/sqrt(beta) + r0*dr0dt*2*s2^2/beta + k*(s-2*s2*c2/sqrt(beta))/beta
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  while iter == 0 || (abs(ds) > 1e-8*s0 && iter < 10)
    xx = 0.5*sqb*s0
    sx = sin(xx) 
    cx = cos(xx)
#    println((r0-k/beta)*2.0*sx*cx/sqb," ",r0*dr0dt*2.*sx^2/beta," ",k*s0/beta," ",h)
# Third derivative:
    yppp = (k-r0*beta)*(cx*cx-sx*sx) - sqb*r0*dr0dt*2.*sx*cx
# Take derivative:
    yp = (-yppp+ k)/beta
# Second derivative:
    ypp = -(r0-k/beta)*2.*sx*cx*sqb  + r0*dr0dt*(cx*cx-sx*sx)
#    y  = (r0-k/beta)*2.0*sx*cx/sqb+r0*dr0dt*2.*sx^2/beta+k*s0/beta -h  # eqn 35
    y  = (-ypp +r0*dr0dt +k*s0)/beta - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = -y/(yp-y/yp*ypp*.5)
    ds = -y/(yp+ds*.5*(ypp+ds*yppp/3.))
    s0 += ds
#    println("s0: ",s0," x: ",x," y: ",y," yp: ",yp," ds: ",ds)
    iter +=1
  end
#  println("iter: ",iter," ds/s: ",ds/s0)
  xx = 0.5*sqb*s0; sx = sin(xx) ; cx = cos(xx)
# Now, compute final values:
  g1bs = 2.*sx*cx/sqb
  g2bs = 2.*sx^2/beta
  f = 1.0 - k/r0*g2bs # eqn (25)
  g = r0*g1bs + r0*dr0dt*g2bs # eqn (27)
  for j=1:3
    x[j] = x0[j]*f+v0[j]*g
  end
  r = norm(x)
  dfdt = -k/r/r0*g1bs
  dgdt = r0/r*(1.0-beta*g2bs+dr0dt*g1bs)
  for j=1:3
    v[j] = x0[j]*dfdt+v0[j]*dgdt
  end
else
  if beta < -1e-15
#    println("Negative beta ",beta)
# Hyperbolic case: (see 8/9/17 notes)
    if s0 == 0.0
# Initial guess:
      s0 = h/r0
    end
    sqb = sqrt(-beta)
# h = r0 * 2.*s2*c2/sqrt(beta) - r0*dr0dt*2*s2^2/beta + k*(s-2*s2*c2/sqrt(beta))/beta
    y = 0.0; yp = 1.0
    iter = 0
    ds = Inf
    while iter == 0 || (abs(ds) > 1e-8*s0 && iter < 10)
      xx = 0.5*sqb*s0; sx = sinh(xx) ; cx = cosh(xx)
#    println((r0-k/beta)*2.0*sx*cx/sqb," ",r0*dr0dt*2.*sx^2/beta," ",k*s0/beta," ",h)
# Third derivative:
      yppp = (k-r0*beta)*(cx*cx+sx*sx) + sqb*r0*dr0dt*2.0*sx*cx
# Take derivative:
# yp = k/beta - (k/beta-r0)*cx - sqrt(-beta)*r0*dr0dt*2*sx*cx
      yp = (-yppp+ k)/beta
# Second derivative:
      ypp = 2.0*(r0-k/beta)*cx*sx*sqb  + r0*dr0dt*(cx*cx+sx*sx)
#    y  = (r0-k/beta)*2.0*sx*cx/sqb+r0*dr0dt*2.*sx^2/beta+k*s0/beta -h  # eqn 35
      y  = (-ypp +r0*dr0dt +k*s0)/beta - h  # eqn 35
# Now, compute fourth-order estimate:
      ds = -y/(yp-y/yp*ypp*.5)
      ds = -y/(yp+ds*.5*(ypp+ds*yppp/3.))
      s0 += ds
#    println("s0: ",s0," x: ",x," y: ",y," yp: ",yp," ds: ",ds)
      iter +=1
    end
#    println("iter: ",iter," ds/s: ",ds/s0)
    xx = 0.5*sqb*s0; sx = sinh(xx) ; cx = cosh(xx)
# Now, compute final values:
    g1bs = 2.0*sx*cx/sqb
    g2bs = -2.0*sx^2/beta
    f = 1.0 - k/r0*g2bs # eqn (25)
    g = r0*g1bs + r0*dr0dt*g2bs # eqn (27)
    for j=1:3
      x[j] = x0[j]*f+v0[j]*g
    end
    r = norm(x)
    dfdt = -k/r/r0*g1bs
    dgdt = r0/r*(1.0-beta*g2bs+dr0dt*g1bs)
    for j=1:3
      v[j] = x0[j]*dfdt+v0[j]*dgdt
    end
  else
#    println("beta = 0")
# Parabolic case: cubic equation
    Q3 = (r0*dr0dt/k)^2 - 2r0/k # = (r0/k)^2*(dr0dt^2-2k/r0) = (r0/k)^2*(dr0dt^2-v0^2) < 0
    R3 = (r0*dr0dt/k)^3-3r0^2*dr0dt/k^2-3*h/k # = (r0/k)^3*dr0dt*(dr0dt^2 - 3/2*v0^2 -3*h*(k/r0)^3/dr0dt/k)
    if R3^2 >= Q3^3  # This should *always* be the case since Q <= 0
#      println("One real root:  Q: ",Q3," R: ",R3)
      A3 = -sign(R3)*(abs(R3)+sqrt(abs(R3^2-Q3^3)))^(1//3)
      if A3 == 0.0
        B3 = 0.0
      else
        B3 = Q3/A3
      end
      s0 = A3 + B3 - r0*dr0dt/k
    else
      println("Error: Q should be negative")
    end
    g1bs = s0
    g2bs = 0.5*s0^2
    f = 1.0 - k/r0*g2bs # eqn (25)
    g = r0*g1bs + r0*dr0dt*g2bs # eqn (27)
    for j=1:3
      x[j] = x0[j]*f+v0[j]*g
    end
    r = norm(x)
    dfdt = -k/r/r0*g1bs
    dgdt = r0/r*(1.0-beta*g2bs+dr0dt*g1bs)
    for j=1:3
      v[j] = x0[j]*dfdt+v0[j]*dgdt
    end
# recompute beta:
#    beta_final = 2.0*k/r-dot(v,v)
#    println("final beta: ",beta_final," rdot: ",dr0dt)
  end
end
#println("s0: ",s0)
return x,v,s0
end
