# Translation of David Hernandez's nbody.c for integrating hiercharical
# system with BH15 integrator.  Please cite Hernandez & Bertschinger (2015)
# if using this in a paper.

const YEAR  = 365.242
const GNEWT = 39.4845/YEAR^2
const NDIM  = 3
const KEPLER_TOL = 1e-8
const third = 1./3.
const alpha0 = 0.0
include("kepler_step.jl")
include("init_nbody.jl")
const pxpr0 = zeros(Float64,3);const  pxpa0=zeros(Float64,3);const  pxpk=zeros(Float64,3);const  pxps=zeros(Float64,3);const  pxpbeta=zeros(Float64,3)
const dxdr0 = zeros(Float64,3);const  dxda0=zeros(Float64,3);const  dxdk=zeros(Float64,3);const  dxdv0 =zeros(Float64,3)
const prvpr0 = zeros(Float64,3);const  prvpa0=zeros(Float64,3);const  prvpk=zeros(Float64,3);const  prvps=zeros(Float64,3);const  prvpbeta=zeros(Float64,3)
const drvdr0 = zeros(Float64,3);const  drvda0=zeros(Float64,3);const  drvdk=zeros(Float64,3);const  drvdv0=zeros(Float64,3)
const vtmp = zeros(Float64,3);const  dvdr0 = zeros(Float64,3);const  dvda0=zeros(Float64,3);const  dvdv0=zeros(Float64,3);const  dvdk=zeros(Float64,3)

const state0 = zeros(Float64,12); state = zeros(Float64,12); delx = zeros(Float64,NDIM); delv = zeros(Float64,NDIM); xcm = zeros(Float64,NDIM); vcm = zeros(Float64,NDIM)

function ttv!(n::Int64,t0::Float64,h::Float64,tmax::Float64,elements::Array{Float64,2},tt::Array{Float64,2},count::Array{Int64,1},dtdq0::Array{Float64,4})
#fcons = open("fcons.txt","w");
m=zeros(Float64,n)
x=zeros(Float64,NDIM,n)
v=zeros(Float64,NDIM,n)
# Fill the transit-timing & jacobian arrays with zeros:
fill!(tt,0.0)
fill!(dtdq0,0.0)
# Counter for transits of each planet:
fill!(count,0)
for i=1:n
  m[i] = elements[i,1]
end
# Initialize the N-body problem using nested hierarchy of Keplerians:
x,v = init_nbody(elements,t0,n)
xprior = copy(x)
vprior = copy(v)
xtransit = copy(x)
vtransit = copy(v)
# Set the time to the initial time:
t = t0
# Set step counter to zero:
istep = 0
# Jacobian for each step (7- 6 elements+mass, n_planets, 7 - 6 elements+mass, n planets):
#jac_step = zeros(Float64,7*n,7*n)
jac_prior = zeros(Float64,7*n,7*n)
jac_transit = zeros(Float64,7*n,7*n)
dtdq = zeros(Float64,7,n)
# Initialize the Jacobian to the identity matrix:
jac_step = eye(Float64,7*n)

# Save the g function, which computes the relative sky velocity dotted with relative position
# between the planets and star:
gsave = zeros(Float64,n)
# Loop over time steps:
dt::Float64 = 0.0
gi = 0.0
while t < t0+tmax
  # Carry out a phi^2 mapping step:
#  phi2!(x,v,h,m,n)
  dh17!(x,v,h,m,n,jac_step)
  # Check to see if a transit may have occured.  Sky is x-y plane; line of sight is z.
  # Star is body 1; planets are 2-nbody (note that this could be modified to see if
  # any body transits another body):
  for i=2:n
    # Compute the relative sky velocity dotted with position:
    gi = g!(i,1,x,v)
    ri = sqrt(x[1,i]^2+x[2,i]^2+x[3,i]^2)
    # See if sign switches, and if planet is in front of star (by a good amount):
    if gi > 0 && gsave[i] < 0 && x[3,i] > 0.25*ri
      # A transit has occurred between the time steps.
      # Approximate the planet-star motion as a Keplerian, weighting over timestep:
      count[i] += 1
#      tt[i,count[i]]=t+findtransit!(i,h,gi,gsave[i],m,xprior,vprior,x,v)
      dt = -gsave[i]*h/(gi-gsave[i])
#      dt = findtransit2!(1,i,h,dt,m,xprior,vprior)
      xtransit .= xprior
      vtransit .= vprior
      jac_transit .= jac_prior
      dt = findtransit2!(1,i,h,dt,m,xtransit,vtransit,jac_transit,dtdq) # 20%
      tt[i,count[i]]=t+dt
      # Save for posterity:
      for k=1:7, p=1:n
        dtdq0[i,count[i],k,p] = dtdq[k,p]
      end
    end
    gsave[i] = gi
  end
  # Save the current state as prior state:
  xprior .= x
  vprior .= v
  jac_prior .= jac_step
  # Increment time by the time step:
  t += h
  # Increment counter by one:
  istep +=1
end
return 
end

function ttv!(n::Int64,t0::Float64,h::Float64,tmax::Float64,elements::Array{Float64,2},tt::Array{Float64,2},count::Array{Int64,1},dlnq::Float64,iq::Int64,jq::Int64)
#fcons = open("fcons.txt","w");
m=zeros(Float64,n)
x=zeros(Float64,NDIM,n)
v=zeros(Float64,NDIM,n)
# Fill the transit-timing array with zeros:
fill!(tt,0.0)
# Counter for transits of each planet:
fill!(count,0)
for i=1:n
  m[i] = elements[i,1]
end
# Allow for perturbations to initial conditions: jq labels body; iq labels phase-space element (or mass)
# iq labels phase-space element (1-3: x; 4-6: v; 7: m)
dq = 0.0
if iq == 7 && dlnq != 0.0
  dq = m[jq]*dlnq
  m[jq] += dq
end
# Initialize the N-body problem using nested hierarchy of Keplerians:
x,v = init_nbody(elements,t0,n)
# Perturb the initial condition by an amount dlnq (if it is non-zero):
if dlnq != 0.0 && iq > 0 && iq < 7
  if iq < 4
    if x[iq,jq] != 0
      dq = x[iq,jq]*dlnq
    else
      dq = dlnq
    end
    x[iq,jq] += dq
  else
  # Same for v
    if v[iq-3,jq] != 0
      dq = v[iq-3,jq]*dlnq
    else
      dq = dlnq
    end
    v[iq-3,jq] += dq
  end
end
xprior = copy(x)
vprior = copy(v)
xtransit = copy(x)
vtransit = copy(v)
# Set the time to the initial time:
t = t0
# Set step counter to zero:
istep = 0
# Jacobian for each step (7 elements+mass, n_planets, 7 elements+mass, n planets):
# Save the g function, which computes the relative sky velocity dotted with relative position
# between the planets and star:
gsave = zeros(Float64,n)
gi  = 0.0
dt::Float64 = 0.0
# Loop over time steps:
while t < t0+tmax
  # Carry out a phi^2 mapping step:
#  phi2!(x,v,h,m,n)
  dh17!(x,v,h,m,n)
  # Check to see if a transit may have occured.  Sky is x-y plane; line of sight is z.
  # Star is body 1; planets are 2-nbody:
  for i=2:n
    # Compute the relative sky velocity dotted with position:
    gi = g!(i,1,x,v)
    ri = sqrt(x[1,i]^2+x[2,i]^2+x[3,i]^2)
    # See if sign switches, and if planet is in front of star (by a good amount):
    if gi > 0 && gsave[i] < 0 && x[3,i] > 0.25*ri
      # A transit has occurred between the time steps.
      # Approximate the planet-star motion as a Keplerian, weighting over timestep:
      count[i] += 1
#      tt[i,count[i]]=t+findtransit!(i,h,gi,gsave[i],m,xprior,vprior,x,v)
      dt = -gsave[i]*h/(gi-gsave[i])
      xtransit .= xprior
      vtransit .= vprior
      dt = findtransit2!(1,i,h,dt,m,xtransit,vtransit)
      tt[i,count[i]]=t+dt
#      tt[i,count[i]]=t+findtransit2!(1,i,h,gi,gsave[i],m,xprior,vprior)
    end
    gsave[i] = gi
  end
  # Save the current state as prior state:
  xprior .=x
  vprior .=v
  # Increment time by the time step:
  t += h
  # Increment counter by one:
  istep +=1
end
return dq
end

# Advances the center of mass of a binary
function centerm!(m::Array{Float64,1},mijinv::Float64,x::Array{Float64,2},v::Array{Float64,2},vcm::Array{Float64,1},delx::Array{Float64,1},delv::Array{Float64,1},i::Int64,j::Int64,h::Float64)
for k=1:NDIM
  x[k,i] +=  m[j]*mijinv*delx[k] + h*vcm[k]
  x[k,j] += -m[i]*mijinv*delx[k] + h*vcm[k]
  v[k,i] +=  m[j]*mijinv*delv[k]
  v[k,j] += -m[i]*mijinv*delv[k]
end
return
end


# Drifts bodies i & j
function driftij!(x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64)
for k=1:NDIM
  x[k,i] += h*v[k,i]
  x[k,j] += h*v[k,j]
end
return
end

# Drifts bodies i & j
function driftij!(x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64,jac_step::Array{Float64,2},nbody::Int64)
indi = (i-1)*7
indj = (j-1)*7
for k=1:NDIM
  x[k,i] += h*v[k,i]
  x[k,j] += h*v[k,j]
  # Now for Jacobian:
  for m=1:7*nbody
    jac_step[indi+k,m] += h*jac_step[indi+3+k,m]
    jac_step[indj+k,m] += h*jac_step[indj+3+k,m]
  end    
end
return
end

# Carries out a Kepler step for bodies i & j
function keplerij!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64)
# The state vector has: 1 time; 2-4 position; 5-7 velocity; 8 r0; 9 dr0dt; 10 beta; 11 s; 12 ds
# Initial state:
state0 = zeros(Float64,12)
# Final state (after a step):
state = zeros(Float64,12)
delx = zeros(Float64,NDIM)
delv = zeros(Float64,NDIM)
#println("Masses: ",i," ",j)
for k=1:NDIM
  state0[1+k     ] = x[k,i] - x[k,j]
  state0[1+k+NDIM] = v[k,i] - v[k,j]
end
gm = GNEWT*(m[i]+m[j])
if gm == 0
  for k=1:NDIM
    x[k,i] += h*v[k,i]
    x[k,j] += h*v[k,j]
  end
else
  # predicted value of s
  state0[11]=0.0
  kepler_step!(gm, h, state0, state)
  for k=1:NDIM
    delx[k] = state[1+k] - state0[1+k]
    delv[k] = state[1+NDIM+k] - state0[1+NDIM+k]
  end
# Advance center of mass:
# Compute COM coords:
  mijinv =1.0/(m[i] + m[j])
  vcm = zeros(Float64,NDIM)
  for k=1:NDIM
    vcm[k] = (m[i]*v[k,i] + m[j]*v[k,j])*mijinv
  end
  centerm!(m,mijinv,x,v,vcm,delx,delv,i,j,h)
end
return
end

# Carries out a Kepler step for bodies i & j
function keplerij!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64,jac_ij::Array{Float64,2})
# The state vector has: 1 time; 2-4 position; 5-7 velocity; 8 r0; 9 dr0dt; 10 beta; 11 s; 12 ds
# Initial state:
state0 = zeros(Float64,12)
# Final state (after a step):
state = zeros(Float64,12)
delx = zeros(Float64,NDIM)
delv = zeros(Float64,NDIM)
# jac_ij should be the Jacobian for going from (x_{0,i},v_{0,i},m_i) &  (x_{0,j},v_{0,j},m_j)
# to  (x_i,v_i,m_i) &  (x_j,v_j,m_j), a 14x14 matrix for the 3-dimensional case. 
# Fill with zeros for now:
fill!(jac_ij,0.0)
for k=1:NDIM
  state0[1+k     ] = x[k,i] - x[k,j]
  state0[1+k+NDIM] = v[k,i] - v[k,j]
end
gm = GNEWT*(m[i]+m[j])
# The following jacobian is just computed for the Keplerian coordinates (i.e. doesn't include
# center-of-mass motion, or scale to motion of bodies about their common center of mass):
jac_kepler = zeros(Float64,7,7)
kepler_step!(gm, h, state0, state, jac_kepler)
for k=1:NDIM
  delx[k] = state[1+k] - state0[1+k]
  delv[k] = state[1+NDIM+k] - state0[1+NDIM+k]
end
# Compute COM coords:
mijinv =1.0/(m[i] + m[j])
xcm = zeros(Float64,NDIM)
vcm = zeros(Float64,NDIM)
mi = m[i]*mijinv # Normalize the masses
mj = m[j]*mijinv
#println("Masses: ",i," ",j)
for k=1:NDIM
  xcm[k] = mi*x[k,i] + mj*x[k,j]
  vcm[k] = mi*v[k,i] + mj*v[k,j]
end
# Compute the Jacobian:
jac_ij[ 7, 7] = 1.0  # the masses don't change with time!
jac_ij[14,14] = 1.0
for k=1:NDIM
   jac_ij[   k,   k] +=   mi
   jac_ij[   k, 3+k] += h*mi
   jac_ij[   k, 7+k] +=   mj
   jac_ij[   k,10+k] += h*mj
   jac_ij[ 3+k, 3+k] +=   mi
   jac_ij[ 3+k,10+k] +=   mj
   jac_ij[ 7+k,   k] +=   mi
   jac_ij[ 7+k, 3+k] += h*mi
   jac_ij[ 7+k, 7+k] +=   mj
   jac_ij[ 7+k,10+k] += h*mj
   jac_ij[10+k, 3+k] +=   mi
   jac_ij[10+k,10+k] +=   mj
end
for l=1:NDIM, k=1:NDIM
# Compute derivatives of \delta x_i with respect to initial conditions:
     jac_ij[   k,   l] += mj*jac_kepler[  k,  l]
     jac_ij[   k, 3+l] += mj*jac_kepler[  k,3+l]
     jac_ij[   k, 7+l] -= mj*jac_kepler[  k,  l]
     jac_ij[   k,10+l] -= mj*jac_kepler[  k,3+l]
# Compute derivatives of \delta v_i with respect to initial conditions:
     jac_ij[ 3+k,   l] += mj*jac_kepler[3+k,  l]
     jac_ij[ 3+k, 3+l] += mj*jac_kepler[3+k,3+l]
     jac_ij[ 3+k, 7+l] -= mj*jac_kepler[3+k,  l]
     jac_ij[ 3+k,10+l] -= mj*jac_kepler[3+k,3+l]
# Compute derivatives of \delta x_j with respect to initial conditions:
     jac_ij[ 7+k,   l] -= mi*jac_kepler[  k,  l]
     jac_ij[ 7+k, 3+l] -= mi*jac_kepler[  k,3+l]
     jac_ij[ 7+k, 7+l] += mi*jac_kepler[  k,  l]
     jac_ij[ 7+k,10+l] += mi*jac_kepler[  k,3+l]
# Compute derivatives of \delta v_j with respect to initial conditions:
     jac_ij[10+k,   l] -= mi*jac_kepler[3+k,  l]
     jac_ij[10+k, 3+l] -= mi*jac_kepler[3+k,3+l]
     jac_ij[10+k, 7+l] += mi*jac_kepler[3+k,  l]
     jac_ij[10+k,10+l] += mi*jac_kepler[3+k,3+l]
end
for k=1:NDIM
# Compute derivatives of \delta x_i with respect to the masses:
   jac_ij[   k, 7] += (x[k,i]+h*v[k,i]-xcm[k]-h*vcm[k]-mj*state[1+k])*mijinv + GNEWT*mj*jac_kepler[  k,7]
   jac_ij[   k,14] += (x[k,j]+h*v[k,j]-xcm[k]-h*vcm[k]+mi*state[1+k])*mijinv + GNEWT*mj*jac_kepler[  k,7]
# Compute derivatives of \delta v_i with respect to the masses:
   jac_ij[ 3+k, 7] += (v[k,i]-vcm[k]-mj*state[4+k])*mijinv + GNEWT*mj*jac_kepler[3+k,7]
   jac_ij[ 3+k,14] += (v[k,j]-vcm[k]+mi*state[4+k])*mijinv + GNEWT*mj*jac_kepler[3+k,7]
# Compute derivatives of \delta x_j with respect to the masses:
   jac_ij[ 7+k, 7] += (x[k,i]+h*v[k,i]-xcm[k]-h*vcm[k]-mj*state[1+k])*mijinv - GNEWT*mi*jac_kepler[  k,7]
   jac_ij[ 7+k,14] += (x[k,j]+h*v[k,j]-xcm[k]-h*vcm[k]+mi*state[1+k])*mijinv - GNEWT*mi*jac_kepler[  k,7]
# Compute derivatives of \delta v_j with respect to the masses:
   jac_ij[10+k, 7] += (v[k,i]-vcm[k]-mj*state[4+k])*mijinv - GNEWT*mi*jac_kepler[3+k,7]
   jac_ij[10+k,14] += (v[k,j]-vcm[k]+mi*state[4+k])*mijinv - GNEWT*mi*jac_kepler[3+k,7]
end
# Advance center of mass & individual Keplerian motions:
centerm!(m,mijinv,x,v,vcm,delx,delv,i,j,h)
return
end

# Drifts all particles:
function drift!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,n::Int64)
for i=1:n, j=1:NDIM
  x[j,i] += h*v[j,i]
end
return
end

# Drifts all particles:
function drift!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,n::Int64,jac_step::Array{Float64,2})
indi = 0
for i=1:n
  indi = (i-1)*7
  for j=1:NDIM
    x[j,i] += h*v[j,i]
    # Now for Jacobian:
    for k=1:7*n
      jac_step[indi+j,k] += h*jac_step[indi+3+j,k]
    end    
  end
end
return
end

function phisalpha!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},alpha::Float64,n::Int64)
# Computes the 4th-order correction:
#function [v] = phisalpha(x,v,h,m,alpha)
#n = size(m,2);
a = zeros(Float64,3,n)
rij = zeros(Float64,3)
aij = zeros(Float64,3)
coeff = alpha*h^3/96*2*GNEWT
fac = 0.0; fac1 = 0.0; fac2 = 0.0; r1 = 0.0; r2 = 0.0; r3 = 0.0
for i=1:n
  for j = i+1:n
    for k=1:3
      rij[k] = x[k,i] - x[k,j]
    end
    r2 = rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3]
    r3 = r2*sqrt(r2)
    for k=1:3
      fac = GNEWT*rij[k]/r3
      a[k,i] -= m[j]*fac
      a[k,j] += m[i]*fac
    end
  end
end
# Next, compute \tilde g_i acceleration vector (this is rewritten
# slightly to avoid reference to \tilde a_i):
for i=1:n
  for j=i+1:n
    for k=1:3
      aij[k] = a[k,i] - a[k,j]
#      aij[k] = 0.0
      rij[k] = x[k,i] - x[k,j]
    end
    r2 = rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3]
    r1 = sqrt(r2)
    ardot = aij[1]*rij[1]+aij[2]*rij[2]+aij[3]*rij[3]
    fac1 = coeff/r1^5
    fac2 = (2*GNEWT*(m[i]+m[j])/r1 + 3*ardot) 
    for k=1:3
#      fac = coeff/r1^5*(rij[k]*(2*GNEWT*(m[i]+m[j])/r1 + 3*ardot) - r2*aij[k])
      fac = fac1*(rij[k]*fac2- r2*aij[k])
      v[k,i] += m[j]*fac
      v[k,j] -= m[i]*fac
    end
  end
end
return
end

function phisalpha!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},alpha::Float64,n::Int64,jac_step::Array{Float64,2})
# Computes the 4th-order correction:
#function [v] = phisalpha(x,v,h,m,alpha)
#n = size(m,2);
a = zeros(Float64,3,n)
dadq = zeros(Float64,3,n,4,n)  # There is no velocity dependence
dotdadq = zeros(Float64,4,n)  # There is no velocity dependence
rij = zeros(Float64,3)
aij = zeros(Float64,3)
coeff = alpha*h^3/96*2*GNEWT
fac = 0.0; fac1 = 0.0; fac2 = 0.0; fac3 = 0.0; r1 = 0.0; r2 = 0.0; r3 = 0.0
for i=1:n-1
  for j = i+1:n
    for k=1:3
      rij[k] = x[k,i] - x[k,j]
    end
    r2 = rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3]
    r3 = r2*sqrt(r2)
    for k=1:3
      fac = GNEWT*rij[k]/r3
      a[k,i] -= m[j]*fac
      a[k,j] += m[i]*fac
      # Mass derivative of acceleration vector (10/6/17 notes):
      # Since there is no velocity dependence, this is fourth parameter.
      # Acceleration of ith particle depends on mass of jth particle:
      dadq[k,i,4,j] -= fac
      dadq[k,j,4,i] += fac
      # x derivative of acceleration vector:
      fac *= 3.0/r2
      # Dot product x_ij.\delta x_ij means we need to sum over components:
      for p=1:3
        dadq[k,i,p,i] += fac*m[j]*rij[p]
        dadq[k,i,p,j] -= fac*m[j]*rij[p]
        dadq[k,j,p,j] += fac*m[i]*rij[p]
        dadq[k,j,p,i] -= fac*m[i]*rij[p]
      end
      # Final term has no dot product, so just diagonal:
      fac = GNEWT/r3
      dadq[k,i,k,i] -= fac*m[j]
      dadq[k,i,k,j] += fac*m[j]
      dadq[k,j,k,j] -= fac*m[i]
      dadq[k,j,k,i] += fac*m[i]
    end
  end
end
# Delete this when finished debugging:
#fill!(dadq,0.0)
# Next, compute \tilde g_i acceleration vector (this is rewritten
# slightly to avoid reference to \tilde a_i):
fill!(jac_step,0.0)
# Note that jac_step[(i-1)*7+k,(j-1)*7+p] is the derivative of the kth coordinate
# of planet i with respect to the pth coordinate of planet j.
indi = 0; indj=0; indd = 0
for i=1:n-1
  indi = (i-1)*7
  for j=i+1:n
    indj = (j-1)*7
    for k=1:3
      aij[k] = a[k,i] - a[k,j]
#      aij[k] = 0.0
      rij[k] = x[k,i] - x[k,j]
    end
    # Compute dot product of r_ij with \delta a_ij:
    fill!(dotdadq,0.0)
    for d=1:n, p=1:4, k=1:3
      dotdadq[p,d] += rij[k]*(dadq[k,i,p,d]-dadq[k,j,p,d])
    end
    r2 = rij[1]*rij[1]+rij[2]*rij[2]+rij[3]*rij[3]
    r1 = sqrt(r2)
    ardot = aij[1]*rij[1]+aij[2]*rij[2]+aij[3]*rij[3]
    fac1 = coeff/r1^5
    fac2 = (2*GNEWT*(m[i]+m[j])/r1 + 3*ardot) 
    for k=1:3
#      fac = coeff/r1^5*(rij[k]*(2*GNEWT*(m[i]+m[j])/r1 + 3*ardot) - r2*aij[k])
      fac = fac1*(rij[k]*fac2- r2*aij[k])
      v[k,i] += m[j]*fac
      v[k,j] -= m[i]*fac
      # Mass derivative (first part is easy):
      jac_step[indi+3+k,indj+7] += fac
      jac_step[indj+3+k,indi+7] -= fac
      # Position derivatives:
      fac *= 5.0/r2
      for p=1:3
        jac_step[indi+3+k,indi+p] -= fac*m[j]*rij[p]
        jac_step[indi+3+k,indj+p] += fac*m[j]*rij[p]
        jac_step[indj+3+k,indj+p] -= fac*m[i]*rij[p]
        jac_step[indj+3+k,indi+p] += fac*m[i]*rij[p]
      end
      # Second mass derivative:
      fac = 2*GNEWT*fac1*rij[k]/r1
      jac_step[indi+3+k,indi+7] += fac*m[j]
      jac_step[indi+3+k,indj+7] += fac*m[j]
      jac_step[indj+3+k,indj+7] -= fac*m[i]
      jac_step[indj+3+k,indi+7] -= fac*m[i]
      #  (There's also a mass term in dadq [x]. See below.)
      # Diagonal position terms:
      fac = fac1*fac2
      jac_step[indi+3+k,indi+k] += fac*m[j]
      jac_step[indi+3+k,indj+k] -= fac*m[j]
      jac_step[indj+3+k,indj+k] += fac*m[i]
      jac_step[indj+3+k,indi+k] -= fac*m[i]
      # Dot product \delta rij terms:
      fac = -2*fac1*(rij[k]*GNEWT*(m[i]+m[j])/(r2*r1)+aij[k])
      for p=1:3
        fac3 = fac*rij[p] + fac1*3.0*rij[k]*aij[p]
        jac_step[indi+3+k,indi+p] += m[j]*fac3
        jac_step[indi+3+k,indj+p] -= m[j]*fac3
        jac_step[indj+3+k,indj+p] += m[i]*fac3
        jac_step[indj+3+k,indi+p] -= m[i]*fac3
      end
      # Diagonal acceleration terms:
      fac = -fac1*r2
      # Duoh.  For dadq, have to loop over all other parameters!
      for d=1:n
        indd = (d-1)*7
        for p=1:3
          jac_step[indi+3+k,indd+p] += fac*m[j]*(dadq[k,i,p,d]-dadq[k,j,p,d])
          jac_step[indj+3+k,indd+p] -= fac*m[i]*(dadq[k,i,p,d]-dadq[k,j,p,d])
        end
        # Don't forget mass-dependent term:
        jac_step[indi+3+k,indd+7] += fac*m[j]*(dadq[k,i,4,d]-dadq[k,j,4,d])
        jac_step[indj+3+k,indd+7] -= fac*m[i]*(dadq[k,i,4,d]-dadq[k,j,4,d])
      end
      # Now, for the final term:  (\delta a_ij . r_ij ) r_ij
      fac = 3.*fac1*rij[k]
      for d=1:n
        indd = (d-1)*7
        for p=1:3
          jac_step[indi+3+k,indd+p] += fac*m[j]*dotdadq[p,d]
          jac_step[indj+3+k,indd+p] -= fac*m[i]*dotdadq[p,d]
        end
        jac_step[indi+3+k,indd+7] += fac*m[j]*dotdadq[4,d]
        jac_step[indj+3+k,indd+7] -= fac*m[i]*dotdadq[4,d]
      end
    end
  end
end
for i=1:n
  indi = (i-1)*7
  for k=1:3
  # Position remains unchanged, so Jacobian of position should be identity matrix:
    jac_step[indi+  k,indi+  k] += 1.0
  # Jacobian of velocity has linear dependence on initial velocity
    jac_step[indi+3+k,indi+3+k] += 1.0
  end
  # Mass remains unchanged:
  jac_step[indi+7,indi+7] += 1.0
end
return
end

# Carries out the DH17 mapping
function dh17!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},n::Int64)
alpha = alpha0
h2 = 0.5*h
# alpha = 0. is similar in precision to alpha=0.25
if alpha != 0.0
  phisalpha!(x,v,h,m,alpha,n)
end
drift!(x,v,h2,n)
for i=1:n-1
  for j=i+1:n
    driftij!(x,v,i,j,-h2)
    keplerij!(m,x,v,i,j,h2)
  end
end
if alpha != 1.0
  phisalpha!(x,v,h,m,2.*(1.-alpha),n)
end
for i=n-1:-1:1
  for j=n:-1:i+1
    keplerij!(m,x,v,i,j,h2)
    driftij!(x,v,i,j,-h2)
  end
end
drift!(x,v,h2,n)
if alpha != 0.0
  phisalpha!(x,v,h,m,alpha,n)
end
return
end

# Used in computing the transit time:
function g!(i::Int64,j::Int64,x::Array{Float64,2},v::Array{Float64,2})
# See equation 8-10 Fabrycky (2008) in Seager Exoplanets book
g = (x[1,j]-x[1,i])*(v[1,j]-v[1,i])+(x[2,j]-x[2,i])*(v[2,j]-v[2,i])
return g
end

# Carries out the DH17 mapping & computes the jacobian:
function dh17!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},n::Int64,jac_step::Array{Float64,2})
h2 = 0.5*h
alpha = alpha0
jac_phi = zeros(Float64,7*n,7*n)
jac_tmp = zeros(Float64,7*n,7*n)
jac_step_ij = zeros(Float64,14,7*n)
# alpha = 0. is similar in precision to alpha=0.25
if alpha != 0.0
  phisalpha!(x,v,h,m,alpha,n,jac_phi)
  jac_step .= jac_phi*jac_step # < 1%
end
drift!(x,v,h2,n,jac_step)
jac_ij = zeros(Float64,14,14)
indi = 0; indj = 0
for i=1:n-1
  indi = (i-1)*7
  for j=i+1:n
    indj = (j-1)*7
    driftij!(x,v,i,j,-h2,jac_step,n)
    keplerij!(m,x,v,i,j,h2,jac_ij) # 21%
    for k1=1:7, k2=1:7*n
      jac_step_ij[  k1,k2]=jac_step[indi+k1,k2]
      jac_step_ij[7+k1,k2]=jac_step[indj+k1,k2]
    end
    jac_step_ij .= *(jac_ij,jac_step_ij)
    for k1=1:7, k2=1:7*n
      jac_step[indi+k1,k2]=jac_step_ij[  k1,k2]
      jac_step[indj+k1,k2]=jac_step_ij[7+k1,k2]
    end
  end
end
if alpha != 1.0
  phisalpha!(x,v,h,m,2.*(1.-alpha),n,jac_phi) # 10%
  jac_step .= jac_phi*jac_step # < 1%
end
indi=0; indj=0
for i=n-1:-1:1
  indi=(i-1)*7
  for j=n:-1:i+1
    indj=(j-1)*7
    keplerij!(m,x,v,i,j,h2,jac_ij) # 23%
    for k2=1:7*n, k1=1:7
      jac_step_ij[  k1,k2]=jac_step[indi+k1,k2]
      jac_step_ij[7+k1,k2]=jac_step[indj+k1,k2]
    end
    jac_step_ij .= *(jac_ij,jac_step_ij)
    for k2=1:7*n, k1=1:7
      jac_step[indi+k1,k2]=jac_step_ij[  k1,k2]
      jac_step[indj+k1,k2]=jac_step_ij[7+k1,k2]
    end
    driftij!(x,v,i,j,-h2,jac_step,n) 
  end
end
drift!(x,v,h2,n,jac_step)
if alpha != 0.0
  phisalpha!(x,v,h,m,alpha,n,jac_phi)
  jac_step .= jac_phi*jac_step # < 1%
end
return
end

function findtransit!(i,h,g1,g2,m,x1,v1,x2,v2)
# Computes the transit time, approximating the motion
# as a Keplerian forward & backward in time, weighted by location in the timestep.
# Initial guess using linear interpolation:
tt = -g1*h/(g2-g1)
dt = 1.0

# Setup state vectors for kepler_step:
s10 = zeros(Float64,12)
s20 = zeros(Float64,12)
# Final state (after a step):
s1 = zeros(Float64,12)
s2 = zeros(Float64,12)
s = zeros(Float64,12)
for k=1:NDIM
  s10[1+k     ] = x1[k,i] - x1[k,1]
  s10[1+k+NDIM] = v1[k,i] - v1[k,1]
  s20[1+k     ] = x2[k,i] - x2[k,1]
  s20[1+k+NDIM] = v2[k,i] - v2[k,1]
end
gm = GNEWT*(m[i]+m[1])
iter = 0
accel1= 0.
accel2= 0.
accel = zeros(Float64,3)
while abs(dt) > 1e-8 && iter < 20
  # Advance planet state at start of step to estimated transit time:
  kepler_step!(gm,    tt, s10, s1)
  # Reverse planet state at end of step to estimated transit time:
  kepler_step!(gm, -h+tt, s20, s2)
  # Weight:
  w = tt/h
  # Compute weighting of states:
  for j=2:7
    s[j] = (1.0-w)*s1[j]+w*s2[j]
  end
  # Compute time offset:
  g = s[2]*s[5]+s[3]*s[6]
  # Compute gravitational acceleration
  r1_3 = norm(s1[2:4])^3
  r2_3 = norm(s2[2:4])^3
  for k=1:3
    accel1 = -gm*s1[k+1]/r1_3
    accel2 = -gm*s2[k+1]/r2_3
    accel[k] = accel1*(1.0-w)+accel2*w
  end
  # Compute derivative of g with respect to time:
  gdot = s[5]^2+s[6]^2+s[2]*accel[1]+s[3]*accel[2]
  # Include time derivatives of interpolation (10/4/17 notes):
  gdot += ((s2[2]-s1[2])*s[5] + (s2[5]-s1[5])*s[2] +(s2[3]-s1[3])*s[6]+(s2[6]-s1[6])*s[3])/h
  # Refine estimate of transit time with Newton's method:
  dt = -g/gdot
  # Add refinement to estimated time:
  tt += dt
  iter +=1
end
# Note: this is the time elapsed *after* the beginning of the timestep:
return tt
end

function findtransit!(i,h,g1,g2,m,x1,v1,x2,v2,dtdqn)
# Computes the transit time, approximating the motion
# as a Keplerian forward & backward in time, weighted by location in the timestep.
# Initial guess using linear interpolation:
tt = -g1*h/(g2-g1)
dt = 1.0

# Vector for computing derivative with respect to the initial and final elements
# of the planet & star -  2 planets with 7 elements/masses at two times (initial & final):
dtdqn = zeros(7,2,2)

# Setup state vectors for kepler_step:
s10 = zeros(Float64,12)
s20 = zeros(Float64,12)
# Final state (after a step):
s1 = zeros(Float64,12)
s2 = zeros(Float64,12)
s = zeros(Float64,12)
for k=1:NDIM
  s10[1+k     ] = x1[k,i] - x1[k,1]
  s10[1+k+NDIM] = v1[k,i] - v1[k,1]
  s20[1+k     ] = x2[k,i] - x2[k,1]
  s20[1+k+NDIM] = v2[k,i] - v2[k,1]
end
gm = GNEWT*(m[i]+m[1])
iter = 0
accel1= 0.
accel2= 0.
accel = zeros(Float64,3)
while abs(dt) > 1e-8 && iter < 20
  # Advance planet state at start of step to estimated transit time:
  kepler_step!(gm,    tt, s10, s1)
  # Reverse planet state at end of step to estimated transit time:
  kepler_step!(gm, -h+tt, s20, s2)
  # Weight:
  w = tt/h
  # Compute weighting of states:
  for j=2:7
    s[j] = (1.0-w)*s1[j]+w*s2[j]
  end
  # Compute time offset:
  g = s[2]*s[5]+s[3]*s[6]
  # Compute gravitational acceleration
  r1_3 = norm(s1[2:4])^3
  r2_3 = norm(s2[2:4])^3
  for k=1:3
    accel1 = -gm*s1[k+1]/r1_3
    accel2 = -gm*s2[k+1]/r2_3
    accel[k] = accel1*(1.0-w)+accel2*w
  end
  # Compute derivative of g with respect to time:
  gdot = s[5]^2+s[6]^2+s[2]*accel[1]+s[3]*accel[2]
  # Include time derivatives of interpolation (10/4/17 notes):
  gdot += ((s2[2]-s1[2])*s[5] + (s2[5]-s1[5])*s[2] +(s2[3]-s1[3])*s[6]+(s2[6]-s1[6])*s[3])/h
  # Refine estimate of transit time with Newton's method:
  dt = -g/gdot
  # Add refinement to estimated time:
  tt += dt
  iter +=1
end
# Now compute derivative of transit time with respect to initial (& final)
# positions & masses of the planet/star:


# Note: this is the time elapsed *after* the beginning of the timestep:
return tt
end

function findtransit2!(i::Int64,j::Int64,h::Float64,tt::Float64,m::Array{Float64,1},x1::Array{Float64,2},v1::Array{Float64,2})
# Computes the transit time, approximating the motion as a fraction of a DH17 step forward in time.
# Initial guess using linear interpolation:
dt = 1.0
iter = 0
r3 = 0.0
gdot = 0.0
x = copy(x1)
v = copy(v1)
while abs(dt) > 1e-8 && iter < 20
  x .= x1
  v .= v1
  # Advance planet state at start of step to estimated transit time:
  dh17!(x,v,tt,m,n)
  # Compute time offset:
  gsky = g!(i,j,x,v)
  # Compute gravitational acceleration in sky plane dotted with sky position:
  gdot = 0.0
  for k=1:n
    if k != i
      r3 = sqrt((x[1,k]-x[1,i])^2+(x[2,k]-x[2,i])^2 +(x[3,k]-x[3,i])^2)^3
      gdot += GNEWT*m[k]*((x[1,k]-x[1,i])*(x[1,j]-x[1,i])+(x[2,k]-x[2,i])*(x[2,j]-x[2,i]))/r3
    end
    if k != j
      r3 = sqrt((x[1,k]-x[1,j])^2+(x[2,k]-x[2,j])^2 +(x[3,k]-x[3,j])^2)^3
      gdot -= GNEWT*m[k]*((x[1,k]-x[1,j])*(x[1,j]-x[1,i])+(x[2,k]-x[2,j])*(x[2,j]-x[2,i]))/r3
    end
  end
  # Compute derivative of g with respect to time:
  gdot += (v[1,j]-v[1,i])^2+(v[2,j]-v[2,i])^2
  # Refine estimate of transit time with Newton's method:
  dt = -gsky/gdot
  # Add refinement to estimated time:
  tt += dt
  iter +=1
end
# Compute derivatives:
#  dh17!(x,v,tt,m,n,jac_step)
# Note: this is the time elapsed *after* the beginning of the timestep:
return tt::Float64
end

function findtransit2!(i::Int64,j::Int64,h::Float64,tt::Float64,m::Array{Float64,1},x1::Array{Float64,2},v1::Array{Float64,2},jac_step::Array{Float64,2},dtdq::Array{Float64,2})
# Computes the transit time, approximating the motion as a fraction of a DH17 step forward in time.
# Also computes the Jacobian of the transit time with respect to the initial parameters, dtdq[7,n].
# Initial guess using linear interpolation:
dt = 1.0
iter = 0
r3 = 0.0
gdot = 0.0
x = copy(x1)
v = copy(v1)
while abs(dt) > 1e-8 && iter < 20
  x .= x1
  v .= v1
  # Advance planet state at start of step to estimated transit time:
  dh17!(x,v,tt,m,n)
  # Compute time offset:
  gsky = g!(i,j,x,v)
  # Compute gravitational acceleration in sky plane dotted with sky position:
  gdot = 0.0
  for k=1:n
    if k != i
      r3 = sqrt((x[1,k]-x[1,i])^2+(x[2,k]-x[2,i])^2 +(x[3,k]-x[3,i])^2)^3
      gdot += GNEWT*m[k]*((x[1,k]-x[1,i])*(x[1,j]-x[1,i])+(x[2,k]-x[2,i])*(x[2,j]-x[2,i]))/r3
    end
    if k != j
      r3 = sqrt((x[1,k]-x[1,j])^2+(x[2,k]-x[2,j])^2 +(x[3,k]-x[3,j])^2)^3
      gdot -= GNEWT*m[k]*((x[1,k]-x[1,j])*(x[1,j]-x[1,i])+(x[2,k]-x[2,j])*(x[2,j]-x[2,i]))/r3
    end
  end
  # Compute derivative of g with respect to time:
  gdot += (v[1,j]-v[1,i])^2+(v[2,j]-v[2,i])^2
  # Refine estimate of transit time with Newton's method:
  dt = -gsky/gdot
  # Add refinement to estimated time:
  tt += dt
  iter +=1
end
# Compute time derivatives:
x = copy(x1)
v = copy(v1)
# Compute dgdt with the updated time step.
dh17!(x,v,tt,m,n,jac_step)
# Compute time offset:
gsky = g!(i,j,x,v)
# Compute gravitational acceleration in sky plane dotted with sky position:
gdot = 0.0
for k=1:n
  if k != i
    r3 = sqrt((x[1,k]-x[1,i])^2+(x[2,k]-x[2,i])^2 +(x[3,k]-x[3,i])^2)^3
    gdot += GNEWT*m[k]*((x[1,k]-x[1,i])*(x[1,j]-x[1,i])+(x[2,k]-x[2,i])*(x[2,j]-x[2,i]))/r3
  end
  if k != j
    r3 = sqrt((x[1,k]-x[1,j])^2+(x[2,k]-x[2,j])^2 +(x[3,k]-x[3,j])^2)^3
    gdot -= GNEWT*m[k]*((x[1,k]-x[1,j])*(x[1,j]-x[1,i])+(x[2,k]-x[2,j])*(x[2,j]-x[2,i]))/r3
  end
end
# Compute derivative of g with respect to time:
gdot += (v[1,j]-v[1,i])^2+(v[2,j]-v[2,i])^2
# Set dtdq to zero:
fill!(dtdq,0.0)
indj = (j-1)*7+1
indi = (i-1)*7+1
for p=1:n
  indp = (p-1)*7
  for k=1:7
    # Compute derivatives:
    dtdq[k,p] = -((jac_step[indj,indp+k]-jac_step[indi,indp+k])*(v[1,j]-v[1,i])+(jac_step[indj+1,indp+k]-jac_step[indi+1,indp+k])*(v[2,j]-v[2,i])+
                  (jac_step[indj+3,indp+k]-jac_step[indi+3,indp+k])*(x[1,j]-x[1,i])+(jac_step[indj+4,indp+k]-jac_step[indi+4,indp+k])*(x[2,j]-x[2,i]))/gdot
  end
end
# Note: this is the time elapsed *after* the beginning of the timestep:
return tt::Float64
end