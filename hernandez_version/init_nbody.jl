include("kepler_init.jl")

# Initializes N-body integration for a planet-parallel hierarchical system 
# (see Hamers & Portugies-Zwart 2016 (HPZ16); Beust 2003).
# We want to define a "simplex" hierarchy which can be decomposed into N-1 Keplerian binaries.
# This can be diagramed with a "mobile" diagram, pioneered by Evans (1968).  Here's an example:
#   Level             |
#    4         _______|_______
#    3   _____|____           |     
#    2  |          |      ____|_____     
#    1  |          |     |       ___|____
#       |          |     |      |        |
#       5          4     3      2        1

# Number of levels:  n_level
# Number of bodies:  n_body
#  - Problem is divided up into N-1 Kepler problems:  there is a single Kepler problem at each level.
#  - For example, in the above "mobile" diagram, 1-2 could be a binary star,
#    while 3 is an interior planets orbiting the stars, and 4/5 are a planet/moon orbiting exterior.
#  - We compute the N-body positions with each Keplerian connection (one at each level), starting
#    at the bottom and moving upwards.
#  
function init_nbody(n_body,elements,t0)
# the "_plane" is to remind us that this is currently plane-parallel, so inclination & Omega are zero
n_level = n_body-1
# Input -
# elements: masses & orbital elements for each Keplerian (in this case, each planet plus star)
# Output -
# x: n_body x 3 array of positions  of each planet.
# v: n_body x 3 array of velocities "   "      "
#
# Read in the orbital elements:
# elements = readdlm("elements.txt",',')
# Initialize N-body for each Keplerian:
# Get the indices:
indices = get_indices_planetary(n_body)
# Set up "A" matrix (Hamers & Portegies-Zwart 2016) which transforms from
# cartesian coordinates to Keplerian orbits (we are using opposite sign
# convention of HPZ16, for example, r_1 = R_2-R_1).
amat = zeros(Float64,n_body,n_body)
# Mass vector:
mass = vcat(elements[:,1])
# Set up array for orbital positions of each Keplerian:
rkepler = zeros(3,nbody)
rdotkepler = zeros(3,nbody)
# Fill in the A matrix & compute the Keplerian elements:
for i=1:n_body-1
  # Sums of masses for two components of Keplerian:
  m1 = 0.0
  m2 = 0.0
  for j=1:n_body
    if indices[i,j] == 1
      m1 += mass[j]
    end
    if indices[i,j] == -1
      m2 += mass[j]
    end
  end
  # Compute Kepler problem: r is a vector of positions of "body" 2 with respect to "body" 1; rdot is velocity vector
  r,rdot = kepler_init(t0,m1+m2,elements[i+1,2:5])
  for j=1:3
    rkepler[j,i] = r[j]
    rdotkepler[j,i] = rdot[j]
  end
  # Now, fill in the A matrix
  for j=1:n_body
    if indices[i,j] == 1
      amat[i,j] = -mass[j]/m1
    end
    if indices[i,j] == -1
      amat[i,j] =  mass[j]/m2
    end
  end
end
mtot = sum(mass)
for j=1:n_body
  amat[n_body,j] = mass[j]/mtot
end
ainv = inv(amat)
# Now, compute the Cartesian coordinates (eqn A6 from HPZ16):
x = zeros(Float64,3,n_body)
v = zeros(Float64,3,n_body)
for i=1:n_body
  for j=1:3
    for k=1:n_body
      x[j,i] += ainv[i,k]*rkepler[j,k]
      v[j,i] += ainv[i,k]*dotrkepler[j,k]
    end
  end
end
#v = *(ainv,rdotkepler)
return x,v
end

function get_indices_planetary(n_body)
# This sets up a planetary-hierarchy index matrix
indices = zeros(Int64,n_body,n_body)
for i=1:n_body-1
 for j=1:i
  indices[i,j ]=-1
 end
 indices[i,i+1]= 1
 indices[n_body,i]=1
end
# This is an example for TRAPPIST-1
# indices = [[-1, 1, 0, 0, 0, 0, 0, 0],  # first two bodies orbit in a binary
#            [-1,-1, 1, 0, 0, 0, 0, 0],  # next planet orbits about these
#            [-1,-1,-1, 1, 0, 0, 0, 0],  # etc...
#            [-1,-1,-1,-1, 1, 0, 0, 0],
#            [-1,-1,-1,-1,-1, 1, 0, 0],
#            [-1,-1,-1,-1,-1,-1, 1, 0],
#            [-1,-1,-1,-1,-1,-1,-1, 1],
#            [ 1, 1, 1, 1, 1, 1, 1, 1]  # center of mass of the system
return indices
end
