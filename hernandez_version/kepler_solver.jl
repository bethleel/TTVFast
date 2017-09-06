# Wisdom & Hernandez version of Kepler solver, but with quartic
# convergence.

function calc_ds_opt(y,yp,ypp,yppp)
# Computes quartic Newton's update to equation y=0 using first through 3rd derivatives.
# Uses techniques outlined in Murray & Dermott for Kepler solver.
# Rearrange to reduce number of divisions:
num = y*yp
den1 = yp*yp-y*ypp*.5
den12 = den1*den1
den2 = yp*den12-num*.5*(ypp*den1-third*num*yppp)
return -y*den12/den2
end

#function kep_elliptic(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,s0::Float64,beta0::Float64)
function kep_elliptic!(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,beta0::Float64,s0::Float64,state::Array{Float64,1})
# Solves equation (35) from Wisdom & Hernandez for the elliptic case.

# First, compute beta:
#r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
#beta = 2.0*k/r0-(v0[1]*v0[1]+v0[2]*v0[2]+v0[3]*v0[3])
#dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
# Define output phase-space variables:
#x=zeros(Float64,3)
#v=zeros(Float64,3)
#x=copy(x0)
#v=copy(v0)
r0inv = inv(r0)
beta0inv = inv(beta0)
# Now, solve for s in elliptical Kepler case:
if beta0 > 1e-15
# Initial guess (if s0 = 0):
  if s0 == 0.0
    s = h*r0inv
  else
    s = copy(s0)
  end
  sqb = sqrt(beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
#  while iter == 0 || (abs(ds) > 1e-4*s && iter < 10)
  while iter == 0 || (abs(ds) > 1e-8*s && iter < 10)
#    xx = 0.5*sqb*s
    xx = sqb*s
#    println("xx: ",xx," r0: ",r0)
    sx = sqb*sin(xx)
    cx = cos(xx)
# Third derivative:
#    yppp = (k-r0*beta0)*(cx*cx-sx*sx) - sqb*r0*dr0dt*2.*sx*cx
#    yppp = (k-r0*beta0)*cx - sqb*r0*dr0dt*sx
    yppp = fac1*cx - fac2*sx
# Take derivative:
    yp = (-yppp+ k)*beta0inv
# Second derivative:
#    ypp = -(r0-k/beta0)*2.*sx*cx*sqb  + r0*dr0dt*(cx*cx-sx*sx)
#    ypp = -(r0-k/beta0)*sx*sqb  + r0*dr0dt*cx
    ypp = fac1*beta0inv*sx + fac2*cx
    y  = (-ypp + fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
#    ds = -y/(yp-y/yp*ypp*.5)
#    ds = -y/(yp+ds*.5*(ypp+ds*yppp/3.))
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
    iter +=1
  end
#  if iter > 1
#    println(iter," ",s," ",s/s0-1," ds: ",ds)
#  end
# Compute error (this is useful if value of tol is adjusted):
#  xx = sqb*s; sx = sin(xx); cx = cos(xx)
#  yppp = (k-r0*beta0)*cx - sqb*r0*dr0dt*sx
#  yp = (-yppp+ k)/beta0
#  ypp = -(r0-k/beta0)*sx*sqb  + r0*dr0dt*cx
#  y  = (-ypp +r0*dr0dt +k*s)/beta0 - h  # eqn 35
#  error = calc_ds(y,yp,ypp,yppp)
#  println("error: ",error)  
# Since we updated s, need to recompute:
  xx = 0.5*sqb*s; sx = sin(xx) ; cx = cos(xx)
# Now, compute final values:
  g1bs = 2.*sx*cx/sqb
  g2bs = 2.*sx^2*beta0inv
  f = 1.0 - k*r0inv*g2bs # eqn (25)
  g = r0*g1bs + fac2*g2bs # eqn (27)
  for j=1:3
#    x[j] = x0[j]*f+v0[j]*g
# Position is components 2-4 of state:
    state[1+j] = x0[j]*f+v0[j]*g
  end
#  r = sqrt(x[1]*x[1]+x[2]*x[2]+x[3]*x[3])
  r = sqrt(state[2]*state[2]+state[3]*state[3]+state[4]*state[4])
  rinv = inv(r)
  dfdt = -k*g1bs*rinv*r0inv
  dgdt = r0*(1.0-beta0*g2bs+dr0dt*g1bs)*rinv
  for j=1:3
# Velocity is components 5-7 of state:
    state[4+j] = x0[j]*dfdt+v0[j]*dgdt
  end
else
  println("Not elliptic ",beta0," x0 ",x0)
end
# recompute beta:
state[8]= r
#drdt = (x[1]*v[1]+x[2]*v[2]+x[3]*v[3])/r
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
#beta= 2.0*k/r-(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
#return x,v,r,drdt,beta,s
return iter
#return nothing
end
# Wisdom & Hernandez version of Kepler solver, but with quartic
# convergence.

function kep_hyperbolic!(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,beta0::Float64,s0::Float64,state::Array{Float64,1})
# Solves equation (35) from Wisdom & Hernandez for the hyperbolic case.

r0inv = inv(r0)
beta0inv = inv(beta0)
# Now, solve for s in hyperbolic Kepler case:
if beta0 < -1e-15
# Initial guess (if s0 = 0):
  if s0 == 0.0
    s = h*r0inv
  else
    s = copy(s0)
  end
  sqb = sqrt(-beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
#  while iter == 0 || (abs(ds) > 1e-4*s && iter < 10)
  while iter == 0 || (abs(ds) > 1e-8*s && iter < 10)
    xx = sqb*s; cx = cosh(xx); sx = sqb*(exp(xx)-cx)
# Third derivative:
    yppp = fac1*cx + fac2*sx
# Take derivative:
# yp = k/beta - (k/beta-r0)*cx - sqrt(-beta)*r0*dr0dt*2*sx*cx
    yp = (-yppp+ k)*beta0inv
# Second derivative:
    ypp = -fac1*beta0inv*sx  + fac2*cx
#    y  = (r0-k/beta)*2.0*sx*cx/sqb+r0*dr0dt*2.*sx^2/beta+k*s0/beta -h  # eqn 35
    y  = (-ypp +fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
#    println("s0: ",s0," x: ",x," y: ",y," yp: ",yp," ds: ",ds)
    iter +=1
  end
#    println("iter: ",iter," ds/s: ",ds/s0)
  xx = 0.5*sqb*s; cx = cosh(xx); sx = exp(xx)-cx
# Now, compute final values:
  g1bs = 2.0*sx*cx/sqb
  g2bs = -2.0*sx^2*beta0inv
  f = 1.0 - k*r0inv*g2bs # eqn (25)
  g = r0*g1bs + fac2*g2bs # eqn (27)
  for j=1:3
    state[1+j] = x0[j]*f+v0[j]*g
  end
  # r = norm(x)
  r = sqrt(state[2]*state[2]+state[3]*state[3]+state[4]*state[4])
  rinv = inv(r)
  dfdt = -k*g1bs*rinv*r0inv
  dgdt = r0*(1.0-beta0*g2bs+dr0dt*g1bs)*rinv
  for j=1:3
# Velocity is components 5-7 of state:
    state[4+j] = x0[j]*dfdt+v0[j]*dgdt
  end
else
  println("Not hyperbolic",beta0," x0 ",x0)
end
# recompute beta:
state[8]= r
#drdt = (x[1]*v[1]+x[2]*v[2]+x[3]*v[3])/r
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
#beta= 2.0*k/r-(v[1]*v[1]+v[2]*v[2]+v[3]*v[3])
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
#return x,v,r,drdt,beta,s
return iter
#return nothing
end
