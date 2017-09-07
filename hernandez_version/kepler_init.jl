include("kepler.jl")

function kepler_init(time,mass,elements)
# Takes orbital elements; returns positions & velocities.
# Right now this is plane-parallel (2D), so only 4 elements specified.
# Compute the semi-major axis in AU:
period = elements[1]
semi = (GNEWT*mass*period^2/4/pi^2)^third
# Elements are given by: period, t0, e*cos(omega), e*sin(omega)
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
fdot = n*(1.0+ecc*cos(f))^2/sqrt1mecc2^3
# Compute the radial distance:
r=semi*(1.0-ecc^2)/(1.0+ecc*cos(f))
rdot = semi*n/sqrt1mecc2*ecc*sin(f)
# For now assume plane-parallel:
inc = pi/2
capomega = pi
# Now, compute the positions
x = zeros(Float64,3)
v = zeros(Float64,3)
x[1]=r*(cos(capomega)*cos(omega+f)-sin(capomega)*sin(omega+f)*cos(inc))
x[2]=r*(sin(capomega)*cos(omega+f)+cos(capomega)*sin(omega+f)*cos(inc))
x[3]=r*sin(omega+f)*sin(inc)
# Compute the velocities:
v[1]=x[1]*rdot/r+r*fdot*(-cos(capomega)*sin(omega+f)-sin(capomega)*cos(omega+f)*cos(inc))
v[2]=x[2]*rdot/r+r*fdot*(-sin(capomega)*sin(omega+f)+cos(capomega)*cos(omega+f)*cos(inc))
v[3]=x[3]*rdot/r+r*fdot*cos(omega+f)*sin(inc)
return x,v
end
