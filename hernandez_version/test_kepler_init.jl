# Runs a simple test of kepler_init.jl

const GNEWT = 39.4845/365.242^2
const third = 1.0/3.0

include("kepler_init.jl")

t0 = 2.4
mass = 1.0
period = 1.5
elements = [1.5,rand()*period,.1*randn(),0.1*randn(),pi/2,pi]
ecc = sqrt(elements[3]^2+elements[4]^2)
ntime = 10000
time = linspace(t0,t0+period,ntime)
xvec = zeros(3,ntime)
vvec = zeros(3,ntime)
vfvec = zeros(3,ntime)
dt = 1e-8
for i=1:ntime
  x,v = kepler_init(time[i],mass,elements)
  if i==1 
    x_ekep,v_ekep = kepler_init_ekep(time[i],mass,elements)
  # Check that these agree:
    println("x-x_ekep: ",maximum(abs.(x-x_ekep))," v-v_ekep: ",maximum(abs.(v-v_ekep)))
  end
  xvec[:,i]=x
  vvec[:,i]=v
  # Compute finite difference velocity:
  xf,vf = kepler_init(time[i]+dt,mass,elements)
  vfvec[:,i]=(xf-x)/dt
end
period = elements[1]
omega=atan2(elements[4],elements[3])
semi = (GNEWT*mass*period^2/4/pi^2)^third
# Compute the focus:
focus = [semi*ecc,0.,0.]
xfocus = semi*ecc*(-cos(omega)-sin(omega))
zfocus = semi*ecc*sin(omega)
# Check that we have an ellipse (we're assuming that the motion is in the x-z plane):
Atot = 0.
dAdt = zeros(ntime)
dx = xvec[1,1]-xvec[1,ntime-1]
dz = xvec[3,1]-xvec[3,ntime-1]
dA = 0.5*(xvec[3,1]*dx-xvec[1,1]*dz)
dAdt[1] = dA/(time[2]-time[1])
Atot += dA
for i=2:ntime
  dx = xvec[1,i]-xvec[1,i-1]
  dz = xvec[3,i]-xvec[3,i-1]
  dA = 0.5*(xvec[3,i]*dx-xvec[1,i]*dz)
  dAdt[i] = dA/(time[i]-time[i-1])
  Atot += dA
end  
println("Total area: ",Atot," ",pi*sqrt(1.-ecc^2)*semi^2," ratio: ",Atot/pi/semi^2/sqrt(1.-ecc^2))
using PyPlot
clf()
plot(time,dAdt)
axis([minimum(time),maximum(time),0.,1.5*maximum(dAdt)])

read(STDIN,Char)
clf()
plot(time,vvec[1,:])
plot(time,vfvec[1,:],".")
plot(time,vvec[1,:]-vfvec[1,:],".")
plot(time,vvec[3,:])
plot(time,vfvec[3,:],".")
plot(time,vvec[3,:]-vfvec[3,:],".")
# Check that velocities match finite difference values
