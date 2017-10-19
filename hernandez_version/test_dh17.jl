# Tests the routine dh17 jacobian:

include("ttv.jl")

#function dh17!(x::Array{Float64,2},v::Array{Float64,2},h::Float64,m::Array{Float64,1},n::Int64,jac_step::Array{Float64,4})


# Next, try computing three-body Keplerian Jacobian:

n = 3
t0 = 7257.93115525
#h  = 0.05
h  = 10.0
tmax = 600.0
dlnq = 1e-4

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

# Take a step:
dh17!(x0,v0,h,m,n)

# Now, copy these to compute Jacobian (so that I don't step
# x0 & v0 forward in time):
x = copy(x0)
v = copy(v0)
m = copy(m0)
# Compute jacobian exactly:
dh17!(x,v,h,m,n,jac_step)
# Save these so that I can compute derivatives numerically:
xsave = copy(x)
vsave = copy(v)
msave = copy(m)

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
    dh17!(x,v,h,m,n)
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
    dh17!(x,v,h,m,n)
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
  dh17!(x,v,h,m,n)
  for i=1:n
    for k=1:3
      jac_step_num[  k,i,7,j] = (x[k,i]-xsave[k,i])/dq
      jac_step_num[3+k,i,7,j] = (v[k,i]-vsave[k,i])/dq
    end
    # Mass unchanged -> identity
    jac_step_num[7,i,7,i] = 1.0
  end
end

# Now, compare the results:
#println(jac_step)
#println(jac_step_num)

for j=1:3
  for i=1:7
    for k=1:3
      println(jac_step[i,j,:,k]," ",jac_step_num[i,j,:,k]," ",jac_step[i,j,:,k]./jac_step_num[i,j,:,k]-1.)
    end
  end
end
