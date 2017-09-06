# Wisdom & Hernandez version of Kepler solver.

function kep_eqn(x0::Array{Float64,1},v0::Array{Float64,1},k::Float64,h::Float64,s0::Float64)
# Solve equation (35):

# First, compute beta:
r0 = norm(x0)
beta = 2.0*k/r0-dot(v0,v0)
dr0dt = dot(x0,v0)/r0
x = zeros(3)
v = zeros(3)
#println("r0: ",r0," dr0dt: ",dr0dt)
# Now, solve for s in elliptical Kepler case:
if beta > 0.0
  if s0 == 0.0
# Initial guess:
    s0 = h/r0
  end
  sqb = sqrt(beta)
# h = r0 * 2.*s2*c2/sqrt(beta) + r0*dr0dt*2*s2^2/beta + k*(s-2*s2*c2/sqrt(beta))/beta
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  while iter == 0 || (abs(ds) > 1e-9*s0 && iter < 10)
    x = 0.5*sqb*s0; sx = sin(x) ; cx = cos(x)
#    println((r0-k/beta)*2.0*sx*cx/sqb," ",r0*dr0dt*2.*sx^2/beta," ",k*s0/beta," ",h)
# Third derivative:
    yppp = (k-r0*beta)*(1.0-2.0*sx*sx) - sqb*r0*dr0dt*2.0*sx*cx
# Take derivative:
    yp = (-yppp+ k)/beta
# Second derivative:
    ypp = -2.0*(r0-k/beta)*cx*sx*sqb  + r0*dr0dt*(1.0-2.*sx^2)
#    y  = (r0-k/beta)*2.0*sx*cx/sqb+r0*dr0dt*2.*sx^2/beta+k*s0/beta -h  # eqn 35
    y  = (-ypp +r0*dr0dt +k*s0)/beta - h  # eqn 35
# Now, compute fourth-order estimate:
#    ds = -y/(yp-y/yp*ypp*.5)
#    ds = -y/(yp+ds*.5*(ypp+ds*yppp/3.))
    ds = -y/yp
    s0 += ds
#    println("s0: ",s0," x: ",x," y: ",y," yp: ",yp," ds: ",ds)
    iter +=1
  end
  println("iter: ",iter," ds: ",ds," ds/s: ",ds/s0)
  x = 0.5*sqb*s0; sx = sin(x) ; cx = cos(x)
# Now, compute final values:
  g1bs = 2.0*sx*cx/sqb
  g2bs = 2.0*sx^2/beta
  f = 1.0 - k/r0*g2bs # eqn (25)
  g = r0*g1bs + r0*dr0dt*g2bs # eqn (27)
  x = x0*f+v0*g
  r = norm(x)
  dfdt = -k/r/r0*g1bs
  dgdt = r0/r*(1.0-beta*g2bs+dr0dt*g1bs)
  v = x0*dfdt+v0*dgdt
else
  if beta < 0
# Hyperbolic case:
  else
# Parabolic case:
  end
end
#println("s0: ",s0)
return x,v,s0
end
