include("kepler.jl")

function kepler_init(time::Float64,mass::Float64,elements::Array{Float64,1})
# Takes orbital elements of a single Keplerian; returns positions & velocities.
# This is 3D), so 6 orbital elements specified, the code returns 3D.  For
# Inclination = pi/2, motion is in X-Z plane; sky plane is X-Y.
# Elements are given by: period, t0, e*cos(omega), e*sin(omega), Inclination, Omega
period = elements[1]
# Compute the semi-major axis in AU (or other units specified by GNEWT):
semi = (GNEWT*mass*period^2/4/pi^2)^third
# Convert to eccentricity & longitude of periastron:
ecc=sqrt(elements[3]^2+elements[4]^2)
omega = atan2(elements[4],elements[3])
# The true anomaly at the time of transit:
f1 = 1.5*pi-omega
# Compute the time of periastron passage:
sqrt1mecc2 = sqrt(1.0-ecc^2)
tp=(elements[2]+period*sqrt1mecc2/2.0/pi*(ecc*sin(f1)/(1.0+ecc*cos(f1))
    -2.0/sqrt1mecc2*atan2(sqrt1mecc2*tan(0.5*f1),1.0+ecc)))
# Compute the mean anomaly
n = 2pi/period
m=n*(time-tp)
# Kepler solver:
f=kepler(m,ecc)
ecosfp1 = 1.0+ecc*cos(f)
fdot = n*ecosfp1^2/sqrt1mecc2^3
# Compute the radial distance:
r=semi*(1.0-ecc^2)/ecosfp1
rdot = semi*n/sqrt1mecc2*ecc*sin(f)
# For now assume plane-parallel:
#inc = pi/2
#capomega = pi
inc = elements[5]
capomega = elements[6]
# Now, compute the positions
x = zeros(Float64,3)
v = zeros(Float64,3)
coscapomega = cos(capomega) ; sincapomega = sin(capomega)
#coscapomega = -1.0 ; sincapomega = 0.0
cosomegapf = cos(omega+f) ; sinomegapf = sin(omega+f) 
cosinc = cos(inc) ; sininc = sin(inc)
#cosinc = 0.0 ; sininc = 1.0
x[1]=r*(coscapomega*cosomegapf-sincapomega*sinomegapf*cosinc)
x[2]=r*(sincapomega*cosomegapf+coscapomega*sinomegapf*cosinc)
x[3]=r*sinomegapf*sininc
rdotonr = rdot/r
# Compute the velocities:
rfdot = r*fdot
v[1]=x[1]*rdotonr+rfdot*(-coscapomega*sinomegapf-sincapomega*cosomegapf*cosinc)
v[2]=x[2]*rdotonr+rfdot*(-sincapomega*sinomegapf+coscapomega*cosomegapf*cosinc)
v[3]=x[3]*rdotonr+rfdot*cosomegapf*sininc
return x,v
end

function kepler_init_ekep(time::Float64,mass::Float64,elements::Array{Float64,1},jac_init::Array{Float64,2})
# Takes orbital elements of a single Keplerian; returns positions & velocities.
# This is 3D), so 6 orbital elements specified, the code returns 3D.  For
# Inclination = pi/2, motion is in X-Z plane; sky plane is X-Y.
# Elements are given by: period, t0, e*cos(omega), e*sin(omega), Inclination, Omega
# Returns the 
period = elements[1]
# Compute the semi-major axis in AU (or other units specified by GNEWT):
semi = (GNEWT*mass*period^2/4/pi^2)^third
dsemidp = 2third*semi/period
dsemidm = third*semi/mass
# Convert to eccentricity & longitude of periastron:
ecc=sqrt(elements[3]^2+elements[4]^2)
omega = atan2(elements[4],elements[3])
# The true anomaly at the time of transit:
f1 = 1.5*pi-omega
# Compute the time of periastron passage:
sqrt1mecc2 = sqrt(1.0-ecc^2)
tp=(elements[2]+period*sqrt1mecc2/2.0/pi*(ecc*sin(f1)/(1.0+ecc*cos(f1))
    -2.0/sqrt1mecc2*atan2(sqrt1mecc2*tan(0.5*f1),1.0+ecc)))
# Compute the mean anomaly
n = 2pi/period
m=n*(time-tp)
# Kepler solver: instead of true anomly, return eccentric anomaly:
ekep=ekepler(m,ecc)
cosekep = cos(ekep); sinekep = sin(ekep)
# Compute the radial distance:
r=semi*(1.0-ecc*cosekep)
inc = elements[5]
capomega = elements[6]
# Now, compute the positions
coscapomega = cos(capomega) ; sincapomega = sin(capomega)
cosomega = cos(omega); sinomega = sin(omega)
cosinc = cos(inc) ; sininc = sin(inc)
# Define rotation matrices (M&D 2.119-2.120):
P1 = [cosomega -sinomega 0.0; sinomega cosomega 0.0; 0.0 0.0 1.0]
P2 = [1.0 0.0 0.0; 0.0 cosinc -sininc; 0.0 sininc cosinc]
P3 = [coscapomega -sincapomega 0.0; sincapomega coscapomega 0.0; 0.0 0.0 1.0]
P321 = P3*P2*P1
# Express x & v in terms of eccentric anomaly (2.121, 2.41, 2.68):
x = P321*semi*[cosekep-ecc;  sqrt1mecc2*sinekep;  0.0]
# Compute the velocities:
v = P321*n*semi^2/r*[-sinekep;  sqrt1mecc2*cosekep;  0.0]
# Now, take derivatives (11/15/2017 notes):

return x,v
end
