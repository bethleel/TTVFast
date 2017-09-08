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
L0=zeros(NDIM)
p0=zeros(NDIM)
xcm0=zeros(NDIM)
ssave = zeros(n,n,3)
# Fill the transit-timing array with zeros:
fill!(tt,0.0)
# Counter for transits of each planet:
fill!(count,0)
# Read in the initial orbital elements from a file (this
# should be moved to an input parameter):
#elements = readdlm("elements.txt",',')
for i=1:n
  m[i] = elements[i,1]
end
# Initialize the N-body problem using nested hierarchy of Keplerians:
x,v = init_nbody(elements,t0,n)
xprior = copy(x)
vprior = copy(v)
#    infig1!(m,x,v,n)
# Compute the conserved quantities:
#E0=consq!(m,x,v,n,L0,p0,xcm0)
#println("E0: ",E0," L0: ",L0," p0: ",p0," xcm0: ",xcm0)
#read(STDIN,Char)
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
#      println("Transit occured: ",i," g1: ",gsave[i]," g2: ",gi," z: ",x[:,i])
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
# Compute the final values of the conserved quantities:
L=zeros(NDIM)
p=zeros(NDIM)
xcm=zeros(NDIM)
#E=consq!(m,x,v,n,L,p,xcm)
#println("dE/E=",(E0-E)/E0)
return 
end

# Advances the center of mass of a binary
function centerm!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},delx::Array{Float64,1},delv::Array{Float64,1},i::Int64,j::Int64,h::Float64)
vcm=zeros(NDIM)
mij =m[i] + m[j]
if mij == 0 
  for k=1:NDIM
    vcm[k] = (v[k,i]+v[k,j])/h
    x[k,i] += h*vcm[k]
    v[k,j] += h*vcm[k]
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


# Initializes the system for a specific dynamical problem: eccentric binary orbiting
# a larger body.  Units in ~AU & yr.
function infig1!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},n::Int64)
a0 = 0.0125
a1 = 1.0
e0 = 0.6
e1 = 0
m[1] = 1e-3
m[2] = 1e-3
m[3] = 1.0
mu0 = (m[1]*m[2])/(m[1]+m[2])
mt0 = m[1]+m[2]
mu1 = (m[3]*mt0)/(m[3]+mt0)
mt1 = m[3]+mt0
x0 = a0*(1+e0)
x1 = a1*(1+e1)
v0 = sqrt(GNEWT*mt0/a0*(1-e0)/(1+e0))
v1 = sqrt(GNEWT*mt1/a1*(1-e1)/(1+e1))
fa0 = m[1]/mt0
fb0 = m[2]/mt0
fa1 = mt0/mt1
fb1 = m[3]/mt1
x[1,1] = -fb1*x1-fb0*x0
x[1,2] = -fb1*x1+fa0*x0
x[1,3] = fa1*x1
v[2,1] = -fb1*v1-fb0*v0
v[2,2] = -fb1*v1+fa0*v0
v[2,3] = fa1*v1
return
end

# Computes conserved quantities:
function consq!(m::Array{Float64,1},x::Array{Float64,2},v::Array{Float64,2},n::Int64,
              L::Array{Float64,1},p::Array{Float64,1},xcm::Array{Float64,1})
fill!(p,0.0) # =zeros(NDIM)
fill!(xcm,0.0) # =zeros(NDIM)
fill!(L,0.0) #  =zeros(NDIM)
E = 0.0
for i=1:n
  if m[i] != 0
    for j=1:NDIM
      p[j] += m[i]*v[j,i];
      xcm[j] += m[i]*x[j,i];
    end
    L[1] += m[i]*(x[2,i]*v[3,i]-x[3,i]*v[2,i]);
    L[2] += m[i]*(x[3,i]*v[1,i]-x[1,i]*v[3,i]);
    L[3] += m[i]*(x[1,i]*v[2,i]-x[2,i]*v[1,i]);
  end
  for j=1:NDIM
    E += 1.0/2.0*m[i]*v[j,i]*v[j,i];
  end
  for j=i+1:n
    if m[j] != 0
      rij = 0.;
      for k=1:NDIM
        rij += (x[k,i]-x[k,j])*(x[k,i]-x[k,j]);
      end
      rij = sqrt(rij);
      E -= GNEWT*m[i]*m[j]/rij;
    end
  end
end
mt=sum(m)
if mt != 0
  for i=1:NDIM
    xcm[i] /= mt;
    p[i] /= mt;
  end
end
println("L: ",L)
println("p: ",p)
println("xcm: ",xcm)
return E
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
#println("Found transit time: ",i," ",tt," iter: ",iter," ",s[2]*s[5]+s[3]*s[6]," ",dt)
# Note: this is the time *after* the beginning of the timestep:
return tt
end
