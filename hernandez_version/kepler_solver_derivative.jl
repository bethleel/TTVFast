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

function kep_elliptic!(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,beta0::Float64,s0::Float64,state::Array{Float64,1})
# Solves equation (35) from Wisdom & Hernandez for the elliptic case.

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
  s0 = copy(s)
  sqb = sqrt(beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
  while iter == 0 || (abs(ds) > 1e-8 && iter < 10)
    xx = sqb*s
    sx = sqb*sin(xx)
    cx = cos(xx)
# Third derivative:
    yppp = fac1*cx - fac2*sx
# Take derivative:
    yp = (-yppp+ k)*beta0inv
# Second derivative:
    ypp = fac1*beta0inv*sx + fac2*cx
    y  = (-ypp + fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
    iter +=1
  end
  if iter > 2
    println(iter," ",s," ",s/s0-1," ds: ",ds)
  end
# Since we updated s, need to recompute:
  xx = 0.5*sqb*s; sx = sin(xx) ; cx = cos(xx)
# Now, compute final values:
  g1bs = 2.*sx*cx/sqb
  g2bs = 2.*sx^2*beta0inv
  f = 1.0 - k*r0inv*g2bs # eqn (25)
  g = r0*g1bs + fac2*g2bs # eqn (27)
  for j=1:3
# Position is components 2-4 of state:
    state[1+j] = x0[j]*f+v0[j]*g
  end
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
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
return iter
end

function kep_elliptic!(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,beta0::Float64,s0::Float64,state::Array{Float64,1},jacobian::Array{Float64,2})
# Computes the Jacobian as well
# Solves equation (35) from Wisdom & Hernandez for the elliptic case.

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
  s0 = copy(s)
  sqb = sqrt(beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
  while iter == 0 || (abs(ds) > 1e-8 && iter < 10)
    xx = sqb*s
    sx = sqb*sin(xx)
    cx = cos(xx)
# Third derivative:
    yppp = fac1*cx - fac2*sx
# Take derivative:
    yp = (-yppp+ k)*beta0inv
# Second derivative:
    ypp = fac1*beta0inv*sx + fac2*cx
    y  = (-ypp + fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
    iter +=1
  end
  if iter > 2
    println(iter," ",s," ",s/s0-1," ds: ",ds)
  end
# Since we updated s, need to recompute:
  xx = 0.5*sqb*s; sx = sin(xx) ; cx = cos(xx)
# Now, compute final values:
  g1bs = 2.*sx*cx/sqb
  g2bs = 2.*sx^2*beta0inv
  f = 1.0 - k*r0inv*g2bs # eqn (25)
  g = r0*g1bs + fac2*g2bs # eqn (27)
  for j=1:3
# Position is components 2-4 of state:
    state[1+j] = x0[j]*f+v0[j]*g
  end
  r = sqrt(state[2]*state[2]+state[3]*state[3]+state[4]*state[4])
  rinv = inv(r)
  dfdt = -k*g1bs*rinv*r0inv
  dgdt = r0*(1.0-beta0*g2bs+dr0dt*g1bs)*rinv
  for j=1:3
# Velocity is components 5-7 of state:
    state[4+j] = x0[j]*dfdt+v0[j]*dgdt
  end
# Now, compute the jacobian:
  fill!(jacobian,0.0)
  compute_jacobian!(h,k,x0,v0,beta0,s,f,g,dfdt,dgdt,cx,sx,g1bs,g2bs,r0,dr0dt,r,jacobian)
else
  println("Not elliptic ",beta0," x0 ",x0)
end
# recompute beta:
state[8]= r
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
# Compute the Jacobian.  jacobian[i,j] is derivative of final state variable q[i]
# with respect to initial state variable q0[j], where q = {x,v} & q0 = {x0,v0}.

return iter
end

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
  s0 = copy(s)
  sqb = sqrt(-beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
  while iter == 0 || (abs(ds) > 1e-8 && iter < 10)
    xx = sqb*s; cx = cosh(xx); sx = sqb*(exp(xx)-cx)
# Third derivative:
    yppp = fac1*cx + fac2*sx
# Take derivative:
    yp = (-yppp+ k)*beta0inv
# Second derivative:
    ypp = -fac1*beta0inv*sx  + fac2*cx
    y  = (-ypp +fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
    iter +=1
  end
  if iter > 2
    #println("iter: ",iter," ds/s: ",ds/s0)
    println(iter," ",s," ",s/s0-1," ds: ",ds)
  end
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
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
return iter
end

function kep_hyperbolic!(x0::Array{Float64,1},v0::Array{Float64,1},r0::Float64,dr0dt::Float64,k::Float64,h::Float64,beta0::Float64,s0::Float64,state::Array{Float64,1},jacobian::Array{Float64,2})
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
  s0 = copy(s)
  sqb = sqrt(-beta0)
  y = 0.0; yp = 1.0
  iter = 0
  ds = Inf
  fac1 = k-r0*beta0
  fac2 = r0*dr0dt
  while iter == 0 || (abs(ds) > 1e-8 && iter < 10)
    xx = sqb*s; cx = cosh(xx); sx = sqb*(exp(xx)-cx)
# Third derivative:
    yppp = fac1*cx + fac2*sx
# Take derivative:
    yp = (-yppp+ k)*beta0inv
# Second derivative:
    ypp = -fac1*beta0inv*sx  + fac2*cx
    y  = (-ypp +fac2 +k*s)*beta0inv - h  # eqn 35
# Now, compute fourth-order estimate:
    ds = calc_ds_opt(y,yp,ypp,yppp)
    s += ds
    iter +=1
  end
  if iter > 2
    #println("iter: ",iter," ds/s: ",ds/s0)
    println(iter," ",s," ",s/s0-1," ds: ",ds)
  end
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
# Now, compute the jacobian:
  fill!(jacobian,0.0)
  compute_jacobian!(h,k,x0,v0,beta0,s,f,g,dfdt,dgdt,cx,sx,g1bs,g2bs,r0,dr0dt,r,jacobian)
else
  println("Not hyperbolic",beta0," x0 ",x0)
end
# recompute beta:
state[8]= r
state[9] = (state[2]*state[5]+state[3]*state[6]+state[4]*state[7])*rinv
# beta is element 10 of state:
state[10] = 2.0*k*rinv-(state[5]*state[5]+state[6]*state[6]+state[7]*state[7])
# s is element 11 of state:
state[11] = s
# ds is element 12 of state:
state[12] = ds
return iter
end

function compute_jacobian!(h,k,x0,v0,beta0,s,f,g,dfdt,dgdt,cx,sx,g1,g2,r0,dr0dt,r,jacobian)
# Compute the Jacobian.  jacobian[i,j] is derivative of final state variable q[i]
# with respect to initial state variable q0[j], where q = {x,v,k} & q0 = {x0,v0,k}.
# Now, compute the Jacobian: (9/18/2017 notes)
g0 = cx^2-sx^2
g3 = (s-g1)/beta0
dsdbeta = 0.5*(2h-r0*g1-r0*s*g0-r0*dr0dt*s*g1-k*(g1-s*g0)/beta0)/(beta0*r)
dsdr0 = -(2k/r0^2*dsdbeta+g1/r+g2*dr0dt/r)
dsdv0 = -(2dr0dt*dsdbeta+g2*r0/r)
dsdk = 2/r0*dsdbeta-g3/r
dbetadr0 = -2k/r0^2
dbetadv0 = -2dr0dt
dbetadk  = 2/r0
# "p" for partial derivative:
pxpr0 = k/r0^2*g2*x0+(g1+dr0dt*g2)*v0
pxpv0 = r0*g2*v0
pxpk  = -g2/r0*v0
pxps  = -k/r0*g1*x0+(r0*g0+r0*dr0dt*g1)*v0
pxpbeta = k/(2beta0*r0)*(s*g1-2g2)*x0+r0/(2beta0)*(s*g0-g1+s*dr0dt*g1-2*dr0dt*g2)*v0
pvpr0 = k*g1/(r*r0^2)*x0+(g0+dr0dt*g1)*v0
pvpv0 = r0/r*g1*v0
pvpk = -g1/(r*r0)*x0
pvps = -k*g0/(r*r0)*x0+r0/r*(-beta0*g1+dr0dt*g0)*v0
pvpbeta = -k/(2beta0*r*r0)*(s*g0-g1)*x0+r0/(2beta0*r)*(-s*beta0*g1+dr0dt*s*g0-dr0dt*g1)*v0
dxdr0 = pxps*dsdr0 + pxpbeta*dbetadr0+pxpr0
dxdv0 = pxps*dsdv0 + pxpbeta*dbetadv0+pxpv0
dxdk  = pxps*dsdk  + pxpbeta*dbetadk +pxpk
dvdr0 = pvps*dsdr0 + pvpbeta*dbetadr0+pvpr0
dvdv0 = pvps*dsdv0 + pvpbeta*dbetadv0+pvpv0
dvdk  = pvps*dsdk  + pvpbeta*dbetadk +pvpk
# Now, compute Jacobian:
for i=1:3
  jacobian[  i,  i] = f
  jacobian[  i,3+i] = g
  jacobian[3+i,  i] = dfdt
  jacobian[3+i,3+i] = dgdt
  jacobian[  i,7] = dxdk[i]
  jacobian[3+i,7] = dvdk[i]
  for j=1:3
    jacobian[  i,  j] += dxdr0[i]*x0[j]/r0
    if dr0dt != 0.0 
      jacobian[  i,3+j] += dxdv0[i]*v0[j]/dr0dt
    end
    jacobian[3+i,  j] += dvdr0[i]*x0[j]/r0
    if dr0dt != 0.0
      jacobian[3+i,3+j] += dvdr0[i]*v0[j]/dr0dt
    end
  end
  jacobian[7,7]=1.0
end
return
end
