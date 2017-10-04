# Translation of David Hernandez's nbody.c for integrating hiercharical
# system with BH15 integrator.  Please cite Hernandez & Bertschinger (2015)
# if using this in a paper.

const YEAR  = 365.242
const GNEWT = 39.4845/YEAR^2
const NDIM  = 3
const third = 1./3.
include("kepler_step.jl")
include("init_nbody.jl")

function ttv!(n::Int64,t0::Float64,h::Float64,tmax::Float64,elements::Array{Float64,2},tt::Array{Float64,2},count::Array{Int64,1})
fcons = open("fcons.txt","w");
m=zeros(n)
x=zeros(NDIM,n)
v=zeros(NDIM,n)
ssave = zeros(n,n,3)
# Fill the transit-timing array with zeros:
fill!(tt,0.0)
# Counter for transits of each planet:
fill!(count,0)
for i=1:n
  m[i] = elements[i,1]
end
# Initialize the N-body problem using nested hierarchy of Keplerians:
x,v = init_nbody(elements,t0,n)
xprior = copy(x)
vprior = copy(v)
# Set the time to the initial time:
t = t0
# Set step counter to zero:
i=0
# Save the g function, which computes the relative sky velocity dotted with relative position
# between the planets and star:
gsave = zeros(n)
# Loop over time steps:
while t < t0+tmax
  # Increment time by the time step:
  t += h
  # Increment counter by one:
  i +=1
  # Carry out a phi^2 mapping step:
  phi2!(x,v,h,m,n)
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
      tt[i,count[i]]=t+findtransit!(i,h,gi,gsave[i],m,xprior,vprior,x,v)
    end
    gsave[i] = gi
  end
  # Save the current state as prior state:
  for i=1:NDIM
    for j=1:n
      xprior[i,j]=x[i,j]
      vprior[i,j]=v[i,j]
    end
  end
end
return 
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

# Carries out a Kepler step for bodies i & j
function keplerij!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64,jac_ij::Array{Float64,2})
# The state vector has: 1 time; 2-4 position; 5-7 velocity; 8 r0; 9 dr0dt; 10 beta; 11 s; 12 ds
# Initial state:
s0 = zeros(Float64,12)
# Final state (after a step):
s = zeros(Float64,12)
delx = zeros(NDIM)
delv = zeros(NDIM)
# jac_ij should be the Jacobian for going from (x_{0,i},v_{0,i},m_i) &  (x_{0,j},v_{0,j},m_j)
# to  (x_i,v_i,m_i) &  (x_j,v_j,m_j), a 14x14 matrix for the 3-dimensional case. 
# Fill with zeros for now:
fill!(jac_ij,0.0)
for k=1:NDIM
  s0[1+k     ] = x[k,i] - x[k,j]
  s0[1+k+NDIM] = v[k,i] - v[k,j]
end
gm = GNEWT*(m[i]+m[j])
# The following jacobian is just computed for the Keplerian coordinates (i.e. doesn't include
# center-of-mass motion, or scale to motion of bodies about their common center of mass):
jac_kepler = zeros(7,7)
kepler_step!(gm, h, s0, s, jac_kepler)
for k=1:NDIM
  delx[k] = s[1+k] - s0[1+k]
  delv[k] = s[1+NDIM+k] - s0[1+NDIM+k]
end
# Compute COM coords:
mijinv =1.0/(m[i] + m[j])
xcm = zeros(NDIM)
vcm = zeros(NDIM)
mi = m[i]*mijinv # Normalize the masses
mj = m[j]*mijinv
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
   for l=1:NDIM
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
# Compute derivatives of \delta x_i with respect to the masses:
   jac_ij[   k, 7] += (x[k,i]+h*v[k,i]-xcm[k]-mj*s[1+k])*mijinv + GNEWT*mj*jac_kepler[  k,7]
   jac_ij[   k,14] += (x[k,j]+h*v[k,j]-xcm[k]+mi*s[1+k])*mijinv + GNEWT*mj*jac_kepler[  k,7]
# Compute derivatives of \delta v_i with respect to the masses:
   jac_ij[ 3+k, 7] += (v[k,i]-vcm[k]-mj*s[4+k])*mijinv + GNEWT*mj*jac_kepler[3+k,7]
   jac_ij[ 3+k,14] += (v[k,j]-vcm[k]+mi*s[4+k])*mijinv + GNEWT*mj*jac_kepler[3+k,7]
# Compute derivatives of \delta x_j with respect to the masses:
   jac_ij[ 7+k, 7] += (x[k,i]+h*v[k,i]-xcm[k]-mj*s[1+k])*mijinv - GNEWT*mi*jac_kepler[  k,7]
   jac_ij[ 7+k,14] += (x[k,j]+h*v[k,j]-xcm[k]+mi*s[1+k])*mijinv - GNEWT*mi*jac_kepler[  k,7]
# Compute derivatives of \delta v_j with respect to the masses:
   jac_ij[10+k, 7] += (v[k,i]-vcm[k]-mj*s[4+k])*mijinv - GNEWT*mi*jac_kepler[3+k,7]
   jac_ij[10+k,14] += (v[k,j]-vcm[k]+mi*s[4+k])*mijinv - GNEWT*mi*jac_kepler[3+k,7]
end
# Advance center of mass & individual Keplerian motions:
centerm!(m,mijinv,x,v,vcm,delx,delv,i,j,h)
return
end

# Carries out a Kepler step for bodies i & j
function keplerij!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},i::Int64,j::Int64,h::Float64)
# The state vector has: 1 time; 2-4 position; 5-7 velocity; 8 r0; 9 dr0dt; 10 beta; 11 s; 12 ds
# Initial state:
s0 = zeros(Float64,12)
# Final state (after a step):
s = zeros(Float64,12)
delx = zeros(NDIM)
delv = zeros(NDIM)
for k=1:NDIM
  s0[1+k     ] = x[k,i] - x[k,j]
  s0[1+k+NDIM] = v[k,i] - v[k,j]
end
gm = GNEWT*(m[i]+m[j])
kepler_step!(gm, h, s0, s)
for k=1:NDIM
  delx[k] = s[1+k] - s0[1+k]
  delv[k] = s[1+NDIM+k] - s0[1+NDIM+k]
end
# Advance center of mass:
# Compute COM coords:
mijinv =1.0/(m[i] + m[j])
vcm = zeros(NDIM)
for k=1:NDIM
  vcm[k] = (m[i]*v[k,i] + m[j]*v[k,j])*mijinv
end
centerm!(m,mijinv,x,v,vcm,delx,delv,i,j,h)
return
end

# Drifts all particles:
function drift!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,n::Int64)
for i=1:n
  for j=1:NDIM
    x[j,i] += h*v[j,i]
  end
end
return
end

# Carries out the phi^2 mapping
function phi2!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},n::Int64)
drift!(x,v,h/2,n)
for i=1:n-1
  for j=i+1:n
    driftij!(x,v,i,j,-h/2)
    keplerij!(m,x,v,i,j,h/2)
  end
end
for i=n-1:-1:1
  for j=n:-1:i+1
    keplerij!(m,x,v,i,j,h/2)
    driftij!(x,v,i,j,-h/2)
  end
end
drift!(x,v,h/2,n)
return
end

# Used in computing the transit time:
function g!(i,j,x,v)
# See equation 8-10 Fabrycky (2008) in Seager Exoplanets book
g = (x[1,j]-x[1,i])*(v[1,j]-v[1,i])+(x[2,j]-x[2,i])*(v[2,j]-v[2,i])
return g
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
while abs(dt) > 1e-8 && iter < 20
  # Advance planet state at start of step to estimated transit time:
  kepler_step!(gm,    tt, s10, s1)
  # Reverse planet state at end of step to estimated transit time:
  kepler_step!(gm, -h+tt, s20, s2)
  # Compute linear weighting of states:
  for j=2:7
    s[j] = (s1[j]*(h-tt)+s2[j]*tt)/h
  end
  # Compute time offset:
  g = s[2]*s[5]+s[3]*s[6]
  # Compute gravitational acceleration
  accel1 = -gm*s1[2:4]/norm(s1[2:4])^3
  accel2 = -gm*s2[2:4]/norm(s2[2:4])^3
  accel = (accel1*(h-tt)+accel2*tt)/h
  # Compute derivative of g with respect to time:
  gdot = s[5]^2+s[6]^2+s[2]*accel[1]+s[3]*accel[2]
  # Refine estimate of transit time with Newton's method:
  dt = -g/gdot
  # Add refinement to estimated time:
  tt += dt
  iter +=1
end
# Note: this is the time elapsed *after* the beginning of the timestep:
return tt
end
