include("kepler.jl")

function kepler_init(time,mass,elements)
# Takes orbital elements of a single Keplerian; returns positions & velocities.
# Right now this is plane-parallel (2D), so only 4 orbital elements specified
# (although the code returns three dimensions:  motion is in X-Z plane).
# Elements are given by: period, t0, e*cos(omega), e*sin(omega)
period = elements[1]
# Compute the semi-major axis in AU (or other units specified by GNEWT):
semi = (GNEWT*mass*period^2/4/pi^2)^third
# Convert to eccentricity & longitude of periastron:
ecc=sqrt(elements[3]^2+elements[4]^2)
omega = atan2(elements[4],elements[3])
# The true anomaly at the time of transit:
f1 = 1.5*pi-omega*pi/180.0
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
inc = pi/2
capomega = pi
# Now, compute the positions
x = zeros(Float64,3)
v = zeros(Float64,3)
#coscapomega = cos(capomega) & sincapomega = sin(capomega)
coscapomega = -1.0 ; sincapomega = 0.0
cosomegapf = cos(omega+f) ; sinomegapf = sin(omega+f) 
#cosinc = cos(inc) ; sininc = sin(inc)
cosinc = 0.0 ; sininc = 1.0
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
