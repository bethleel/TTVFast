# Tests the routine dh17 jacobian:

include("ttv.jl")

#function dh17!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},n::Int64,jac_step::Array{Float64,4})


# Next, try computing three-body Keplerian Jacobian:

n = 3
#n = 2
t0 = 7257.93115525
#h  = 0.05
h  = 0.15
tmax = 600.0
dlnq = 1e-6

nstep = 10

elements = readdlm("elements.txt",',')
elements[2,1] = 1.0
elements[3,1] = 1.0

m =zeros(n)
x0=zeros(3,n)
v0=zeros(3,n)

# Predict values of s:

jac_step = zeros(7,n,7,n)
# Initialize with identity matrix:
for i=1:7
  for j=1:n
    jac_step[i,j,i,j] = 1.0
  end
end

for k=1:n
  m[k] = elements[k,1]
end
m0 = copy(m)

x0,v0 = init_nbody(elements,t0,n)

# Tilt the orbits a bit:
x0[2,1] = 5e-1*sqrt(x0[1,1]^2+x0[3,1]^2)
x0[2,2] = -5e-1*sqrt(x0[1,2]^2+x0[3,2]^2)
x0[2,3] = -5e-1*sqrt(x0[1,2]^2+x0[3,2]^2)
v0[2,1] = 5e-1*sqrt(v0[1,1]^2+v0[3,1]^2)
v0[2,2] = -5e-1*sqrt(v0[1,2]^2+v0[3,2]^2)
v0[2,3] = -5e-1*sqrt(v0[1,2]^2+v0[3,2]^2)

# Take a single step (so that we aren't at initial coordinates):
dh17!(x0,v0,h,m,n)

# Now, copy these to compute Jacobian (so that I don't step
# x0 & v0 forward in time):
x = copy(x0)
v = copy(v0)
m = copy(m0)
# Compute jacobian exactly:
for istep=1:nstep
  dh17!(x,v,h,m,n,jac_step)
end
# Save these so that I can compute derivatives numerically:
xsave = copy(x)
vsave = copy(v)
msave = copy(m)
## Check that we have agreement:
#xtest = copy(x0)
#vtest = copy(v0)
#m = copy(m0)
#dh17!(xtest,vtest,h,m,n)
#println("x/v difference: ",x-xtest,v-vtest)

# Now compute numerical derivatives:
jac_step_num = zeros(7,n,7,n)
# Vary the initial parameters of planet j:
for j=1:n
  # Vary the initial phase-space elements:
  for jj=1:3
  # Initial positions, velocities & masses:
    x = copy(x0)
    v = copy(v0)
    m = copy(m0)
    dq = dlnq * x[jj,j]
    if x[jj,j] != 0.0
      x[jj,j] +=  dq
    else
      dq = dlnq
      x[jj,j] = dq
    end
    for istep=1:nstep
      dh17!(x,v,h,m,n)
    end
  # Now x & v are final positions & velocities after time step
    for i=1:n
      for k=1:3
        jac_step_num[  k, i, jj, j] = (x[k,i]-xsave[k,i])/dq
        jac_step_num[3+k, i, jj, j] = (v[k,i]-vsave[k,i])/dq
      end
    end
    x=copy(x0)
    v=copy(v0)
    m=copy(m0)
    dq = dlnq * v[jj,j]
    if v[jj,j] != 0.0
      v[jj,j] +=  dq
    else
      dq = dlnq
      v[jj,j] = dq
    end
    for istep=1:nstep
      dh17!(x,v,h,m,n)
    end
    for i=1:n
      for k=1:3
        jac_step_num[  k,i,3+jj,j] = (x[k,i]-xsave[k,i])/dq
        jac_step_num[3+k,i,3+jj,j] = (v[k,i]-vsave[k,i])/dq
      end
    end
  end
# Now vary mass of planet:
  x=copy(x0)
  v=copy(v0)
  m=copy(m0)
  dq = m[j]*dlnq
  m[j] += dq
  for istep=1:nstep
    dh17!(x,v,h,m,n)
  end
  for i=1:n
    for k=1:3
      jac_step_num[  k,i,7,j] = (x[k,i]-xsave[k,i])/dq
      jac_step_num[3+k,i,7,j] = (v[k,i]-vsave[k,i])/dq
    end
  end
  # Mass unchanged -> identity
  jac_step_num[7,j,7,j] = 1.0
end

# Now, compare the results:
#println(jac_step)
#println(jac_step_num)

for j=1:n
  for i=1:7
    for k=1:n
      println(i," ",j," ",k," ",jac_step[i,j,:,k]," ",jac_step_num[i,j,:,k]," ",jac_step[i,j,:,k]./jac_step_num[i,j,:,k]-1.)
#      println(i," ",j," 7 ",k," ",jac_step[i,j,7,k]," ",jac_step_num[i,j,7,k]," ",jac_step[i,j,7,k]./jac_step_num[i,j,7,k]-1.)
    end
  end
end

jacmax = 0.0

for i=1:7, j=1:3, k=1:7, l=1:3
  if jac_step[i,j,k,l] != 0
    diff = abs(jac_step_num[i,j,k,l]/jac_step[i,j,k,l]-1.0)
    if diff > jacmax
      jacmax = diff
    end
  end
end

println("Maximum fractional error: ",jacmax)

