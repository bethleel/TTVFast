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
#      tt[i,count[i]]=t+findtransit2!(i,h,gi,gsave[i],m,xprior,vprior,x,v)
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
function centerm!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},delx::Array{Float64,1},delv::Array{Float64,1},i::Int64,j::Int64,h::Float64)
vcm=zeros(NDIM)
mij =m[i] + m[j]
if mij == 0 
  for k=1:NDIM
    vcm[k] = (v[k,i]+v[k,j])/2
    x[k,i] += h*vcm[k]
    x[k,j] += h*vcm[k]
  end
else
  for k=1:NDIM
    vcm[k] = (m[i]*v[k,i] + m[j]*v[k,j])/mij
    x[k,i] +=  m[j]/mij*delx[k] + h*vcm[k]
    x[k,j] += -m[i]/mij*delx[k] + h*vcm[k]
    v[k,i] +=  m[j]/mij*delv[k]
    v[k,j] += -m[i]/mij*delv[k]
  end
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
if gm == 0
  for k=1:NDIM
    x[k,i] += h*v[k,i]
    x[k,j] += h*v[k,j]
  end
else
  kepler_step!(gm, h, s0, s)
  for k=1:NDIM
    delx[k] = s[1+k] - s0[1+k]
    delv[k] = s[1+NDIM+k] - s0[1+NDIM+k]
  end
# Advance center of mass:
  centerm!(m,x,v,delx,delv,i,j,h)
end
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

function findtransit2!(i,h,g1,g2,m,x1,v1,x2,v2)
# Computes the transit time, approximating the motion
# as a Keplerian forward & backward in time, weighted by location in the timestep.
# Initial guess using linear interpolation:
tt = -g1*h/(g2-g1)
dt1 = 1.0
dt2 = 1.0

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
# estimate transit time going forward:
tt1 = tt
# estimate transit time going backwards:
tt2 = tt
while (abs(dt1) > 1e-8 || abs(dt2) > 1e-8) && iter < 20
  # Advance planet state at start of step to estimated transit time:
  kepler_step!(gm,    tt1, s10, s1)
  # Reverse planet state at end of step to estimated transit time:
  kepler_step!(gm, -h+tt2, s20, s2)
  # Compute linear weighting of states:
#  for j=2:7
#    s[j] = (s1[j]*(h-tt)+s2[j]*tt)/h
#  end
  # Compute time offset:
  g1 = s1[2]*s1[5]+s1[3]*s1[6]
  g2 = s2[2]*s2[5]+s2[3]*s2[6]
  # Compute gravitational acceleration
  accel1 = -gm*s1[2:4]/norm(s1[2:4])^3
  accel2 = -gm*s2[2:4]/norm(s2[2:4])^3
#  accel = (accel1*(h-tt)+accel2*tt)/h
  # Compute derivative of g with respect to time:
  gdot1 = s1[5]^2+s1[6]^2+s1[2]*accel1[1]+s1[3]*accel1[2]
  gdot2 = s2[5]^2+s2[6]^2+s2[2]*accel2[1]+s2[3]*accel2[2]
  # Refine estimate of transit time with Newton's method:
  dt1 = -g1/gdot1
  # Add refinement to estimated time:
  tt1 += dt1
  tt2 += dt1
  iter +=1
end
# Weight the forwards & backwards transit times:
w1,dw1 = weight_time(tt1/h)
w2,dw2 = weight_time(1.0-tt2/h)
#tt = (tt1*tt2+(h-tt2)*tt1)/(h-tt2+tt1)
tt = (w1*tt1+t2*tt2)/(w1+w2)
# If tt1 ~ 0, then (h-tt2)*tt1/(h-tt2) ~ tt1
# If tt2 ~ h, then tt1*tt2/tt1 ~ tt2
# Note: this is the time elapsed *after* the beginning of the timestep:
return tt
end

function weight_time(x)
@assert(x >= 0.0)
@assert(x <= 1.0)
if x < 0.5
  return 1.0-2.*x^2,-4.*x
else
  return 2.*(1.0-x)^2,-4.*(1.0-x)
end
end
