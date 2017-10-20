using PyPlot
#include("kepler_solver_derivative.jl")
include("ttv.jl")

# Define a constant of 1/3:
#const third = 1.0/3.0

function test_elliptic_derivative(dlnq)
# Call as: save,jac_num,jacobian=test_elliptic_derivative(1e-4)
# This routine runs a test of the kep_elliptic_jacobian function in kepler_solver.jl
# Define the central force constant in terms of AU and days:
k = (2.0*pi/365.25)^2
# Initial position at 1 AU:
x0 = [0.01,1.0,0.01]
#r0 = norm(x0)
# Circular velocity:
r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
vcirc = sqrt(k/r0)
# Define initial velocity at apastron:
v0 = [.9*vcirc,0.01*vcirc,0.01*vcirc]  # The eccentricity is about ~2(1-v0/vcirc).
dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
h = 100.0 # 18-day timesteps

s0::Float64 = 0.0
x = zeros(Float64,3)
v = zeros(Float64,3)
s::Float64 = 0.0
t = 0.0
ssave = zeros(2)
nsteps = 1
xsave=zeros(Float64,12,nsteps)  # First is time (1,:); next three are position (2-4,:); next three velocity (5-7,:); then r (8,:), drdt (9,:); then beta (10,:); finally s & ds (11:12,:)
# Create a variable to store the state at the end of a step:
state=zeros(Float64,12)
state_diff=zeros(Float64,12)

xsave[1,1]=0.0
for jj=1:3
  xsave[jj+1,1]=x0[jj]
  xsave[jj+4,1]=v0[jj]
end
# Save beta:
xsave[8,1]=r0
xsave[9,1]=dr0dt
beta0 = 2.0*k/r0-dot(v0,v0)
xsave[10,1]=beta0
jacobian=zeros(Float64,7,7)
iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state,jacobian)
state_save = copy(state)
println("Initial conditions: ",x0,v0)
iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state)
println("Final state: ",state[2:7])

#read(STDIN,Char)
s = state[11]
# Now, compute the Jacobian numerically:
jac_num = zeros(Float64,7,7)
#dlnq = 1e-4
# jac_num[i,j]: derivative of (x_i,v_i,k) with respect to (x_{0,j},v_{0,j},k):
for j=1:3
  x0save = copy(x0)
  dq = dlnq * x0[j]
  if x0[j] != 0.0
    x0[j] +=  dq
  else
    dq = dlnq
    x0[j] = dq
  end
  # Recompute quantities:
  r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
  dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
  beta0 = 2.0*k/r0-dot(v0,v0)
  iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state_diff)
  x0 = copy(x0save)
  for i=1:3
    jac_num[  i,  j] = (state_diff[1+i]-state[1+i])/dq
    jac_num[3+i,  j] = (state_diff[4+i]-state[4+i])/dq
  end
  v0save = copy(v0)
  dq = dlnq * v0[j]
  if v0[j] != 0.0
    v0[j] +=  dq
  else
    dq = dlnq
    v0[j] = dq
  end
  r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
  dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
  beta0 = 2.0*k/r0-dot(v0,v0)
  iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state_diff)
  v0 = copy(v0save)
  for i=1:3
    jac_num[  i,3+j] = (state_diff[1+i]-state[1+i])/dq
    jac_num[3+i,3+j] = (state_diff[4+i]-state[4+i])/dq
  end
  # Now vary mass:
  ksave = copy(k)
  dq = k*dlnq
  k += dq
  r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
  dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
  beta0 = 2.0*k/r0-dot(v0,v0)
  iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state_diff)
  for i=1:3
    jac_num[  i,7] = (state_diff[1+i]-state[1+i])/dq
    jac_num[3+i,7] = (state_diff[4+i]-state[4+i])/dq
  end
  k = copy(ksave)
  jac_num[7,  7] = 1.0
end
return xsave,jac_num,jacobian
end

# First try:
xsave,jac_num1,jacobian=test_elliptic_derivative(1e-7)
# Second try:
xsave,jac_num2,jacobian=test_elliptic_derivative(1e-6)


println("Jac dlnq=1e-7 ")
#jac_num1
println("Jac dlnq=1e-6 ")
#jac_num2
#jacobian
#jac_num1-jac_num2
println(jac_num1./jacobian-1.)

# Next, try computing two-body Keplerian Jacobian:

n = 2
t0 = 7257.93115525
h  = 0.05
tmax = 600.0
dlnq = 1e-5

elements = readdlm("elements.txt",',')
elements[2,1] = 0.75
#elements[1,1] = 1e-5

m =zeros(n)
x0=zeros(NDIM,n)
v0=zeros(NDIM,n)

for k=1:n
  m[k] = elements[k,1]
end

x0,v0 = init_nbody(elements,t0,n)
# Tilt the orbits a bit:
x0[2,1] = 5e-1*sqrt(x0[1,1]^2+x0[3,1]^2)
x0[2,2] = -5e-1*sqrt(x0[1,2]^2+x0[3,2]^2)
v0[2,1] = 5e-1*sqrt(v0[1,1]^2+v0[3,1]^2)
v0[2,2] = -5e-1*sqrt(v0[1,2]^2+v0[3,2]^2)
# Reduce the masses to make it hyperbolic:
m *= 1e-1


jac_ij = zeros(14,14)
i=1 ; j=2
x = copy(x0) ; v=copy(v0)
# Predict values of s:
keplerij!(m,x,v,i,j,h,jac_ij)
x0 = copy(x) ; v0 = copy(v)
keplerij!(m,x,v,i,j,h,jac_ij)

# xtest = copy(x0) ; vtest=copy(v0)
# keplerij!(m,xtest,vtest,i,j,h,jac_ij)
# println("Test of jacobian vs. none: ",maximum(abs(x-xtest)),maximum(abs(v-vtest)))

# Now, compute the derivatives numerically:
jac_ij_num = zeros(14,14)
xsave = copy(x)
vsave = copy(v)
msave = copy(m)

for jj=1:3
  # Initial positions, velocities & masses:
  x = copy(x0)
  v = copy(v0)
  m = copy(msave)
  dq = dlnq * x[jj,i]
  if x[jj,i] != 0.0
    x[jj,i] +=  dq
  else
    dq = dlnq
    x[jj,i] = dq
  end
  keplerij!(m,x,v,i,j,h)
  # Now x & v are final positions & velocities after time step
  for k=1:3
    jac_ij_num[   k,  jj] = (x[k,i]-xsave[k,i])/dq
    jac_ij_num[ 3+k,  jj] = (v[k,i]-vsave[k,i])/dq
    jac_ij_num[ 7+k,  jj] = (x[k,j]-xsave[k,j])/dq
    jac_ij_num[10+k,  jj] = (v[k,j]-vsave[k,j])/dq
  end
  x=copy(x0)
  v=copy(v0)
  m=copy(msave)
  dq = dlnq * v[jj,i]
  if v[jj,i] != 0.0
    v[jj,i] +=  dq
  else
    dq = dlnq
    v[jj,i] = dq
  end
  keplerij!(m,x,v,i,j,h)
  for k=1:3
    jac_ij_num[   k,3+jj] = (x[k,i]-xsave[k,i])/dq
    jac_ij_num[ 3+k,3+jj] = (v[k,i]-vsave[k,i])/dq
    jac_ij_num[ 7+k,3+jj] = (x[k,j]-xsave[k,j])/dq
    jac_ij_num[10+k,3+jj] = (v[k,j]-vsave[k,j])/dq
  end
end
# Now vary mass of inner planet:
x=copy(x0)
v=copy(v0)
m=copy(msave)
dq = m[i]*dlnq
m[i] += dq
keplerij!(m,x,v,i,j,h)
for k=1:3
  jac_ij_num[   k,7] = (x[k,i]-xsave[k,i])/dq
  if k == 2
    println("varying mass of i: ",x[k,i]," ",xsave[k,i]," ",dq," ",jac_ij[k,7]," ",jac_ij_num[k,7])
  end
  jac_ij_num[ 3+k,7] = (v[k,i]-vsave[k,i])/dq
  jac_ij_num[ 7+k,7] = (x[k,j]-xsave[k,j])/dq
  if k == 2
    println("varying mass of i: ",x[k,j]," ",xsave[k,j]," ",dq," ",jac_ij[7+k,7]," ",jac_ij_num[7+k,7])
  end
  jac_ij_num[10+k,7] = (v[k,j]-vsave[k,j])/dq
end
# The mass doesn't change:
jac_ij_num[7,7] =  1.0
for jj=1:3
  # Now vary parameters of outer planet:
  x = copy(x0)
  v = copy(v0)
  m = copy(msave)
  dq = dlnq * x[jj,j]
  if x[jj,j] != 0.0
    x[jj,j] +=  dq
  else
    dq = dlnq
    x[jj,j] = dq
  end
  keplerij!(m,x,v,i,j,h)
  for k=1:3
    jac_ij_num[   k,7+jj] = (x[k,i]-xsave[k,i])/dq
    jac_ij_num[ 3+k,7+jj] = (v[k,i]-vsave[k,i])/dq
    jac_ij_num[ 7+k,7+jj] = (x[k,j]-xsave[k,j])/dq
    jac_ij_num[10+k,7+jj] = (v[k,j]-vsave[k,j])/dq
  end
  x=copy(x0)
  v=copy(v0)
  m=copy(msave)
  dq = dlnq * v[jj,j]
  if v[jj,j] != 0.0
    v[jj,j] +=  dq
  else
    dq = dlnq
    v[jj,j] = dq
  end
  keplerij!(m,x,v,i,j,h)
  for k=1:3
    jac_ij_num[   k,10+jj] = (x[k,i]-xsave[k,i])/dq
    jac_ij_num[ 3+k,10+jj] = (v[k,i]-vsave[k,i])/dq
    jac_ij_num[ 7+k,10+jj] = (x[k,j]-xsave[k,j])/dq
    jac_ij_num[10+k,10+jj] = (v[k,j]-vsave[k,j])/dq
  end
end
# Now vary mass of outer planet:
x = copy(x0)
v = copy(v0)
m = copy(msave)
dq = m[j]*dlnq
m[j] += dq
keplerij!(m,x,v,i,j,h)
for k=1:3
  jac_ij_num[   k,14] = (x[k,i]-xsave[k,i])/dq
  if k == 2
    println("varying mass of j: ",x[k,i]," ",xsave[k,i]," ",dq," ",jac_ij[k,14]," ",jac_ij_num[k,14])
  end
  jac_ij_num[ 3+k,14] = (v[k,i]-vsave[k,i])/dq
  jac_ij_num[ 7+k,14] = (x[k,j]-xsave[k,j])/dq
  if k == 2
    println("varying mass of j: ",x[k,j]," ",xsave[k,j]," ",dq," ",jac_ij[7+k,14]," ",jac_ij_num[7+k,14])
  end
  jac_ij_num[10+k,14] = (v[k,j]-vsave[k,j])/dq
end
# The mass doesn't change:
jac_ij_num[14,14] =  1.0

println(jac_ij)
println(jac_ij_num)
