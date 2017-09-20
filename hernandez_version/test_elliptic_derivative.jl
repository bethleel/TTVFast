using PyPlot
include("kepler_solver_derivative.jl")

# Define a constant of 1/3:
const third = 1.0/3.0

function test_elliptic_derivative()
# This routine runs a test of the kep_elliptic_jacobian function in kepler_solver.jl
# Define the central force constant in terms of AU and days:
k = (2.0*pi/365.25)^2
# Initial position at 1 AU:
x0 = [0.01,1.0,0.01]
#r0 = norm(x0)
r0 = sqrt(x0[1]*x0[1]+x0[2]*x0[2]+x0[3]*x0[3])
# Circular velocity:
vcirc = sqrt(k/r0)
# Define initial velocity at apastron:
v0 = [.9*vcirc,0.01*vcirc,0.01*vcirc]  # The eccentricity is about ~2(1-v0/vcirc).
dr0dt = (x0[1]*v0[1]+x0[2]*v0[2]+x0[3]*v0[3])/r0
h = 10.0 # 18-day timesteps

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
for j=1:3
  xsave[j+1,1]=x0[j]
  xsave[j+4,1]=v0[j]
end
# Save beta:
xsave[8,1]=r0
xsave[9,1]=dr0dt
beta0 = 2.0*k/r0-dot(v0,v0)
xsave[10,1]=beta0
jacobian=zeros(Float64,7,7)
iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state,jacobian)
println("Initial conditions: ",x0,v0)
iter = kep_elliptic!(x0,v0,r0,dr0dt,k,h,beta0,s0,state)
println("Final state: ",state[2:7])
read(STDIN,Char)
s = state[11]
# Now, compute the Jacobian numerically:
jac_num = zeros(Float64,7,7)
dlnq = 1e-4
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
  for i=1:3
    jac_num[  i,7] = (state_diff[1+i]-state[1+i])/dq
    jac_num[3+i,7] = (state_diff[4+i]-state[4+i])/dq
  end
  k = copy(ksave)
  jac_num[7,  7] = 1.0
end
return xsave,jac_num,jacobian
end
