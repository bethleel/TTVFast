# Test Kepler solver:
using CGS
using PyPlot

type Kepler_step
 x0::Vector{Float64}
 v0::Vector{Float64}
 r0::Float64
 dr0dt::Float64
 x::Vector{Float64}
 v::Vector{Float64}
 r::Float64
 h::Float64
 sqb::Float64
 y::Float64
 yp::Float64
 ypp::Float64
 yppp::Float64
 xx::Float64
 cx::Float64
 sx::Float64
 beta::Float64
 beta_final::Float64
 f::Float64
 dfdt::Float64
 g1bs::Float64
 g2bs::Float64
 g::Float64
 dgdt::Float64
 A3::Float64
 B3::Float64
 Q3::Float64
 R3::Float64
 ds::Float64
 iter::Int64
 Kepler_step() = new(zeros(3),zeros(3),0.0,0.0,zeros(3),zeros(3),0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0)
end

function run_test()

nsteps = 1000000
xsave = zeros(3,nsteps)
vsave = zeros(3,nsteps)
ssave = zeros(nsteps)
tsave = zeros(nsteps)

#y=0.0
#yp=0.0
#ypp=0.0
#yppp=0.0
#xx=0.0
#cx=0.0
#sx=0.0
#r0=0.0
#dr0dt=0.0
#beta=0.0
#beta_final=0.0
#sqb=0.0
#f=0.0
#g=0.0
#g1bs=0.0
#g2bs=0.0
#dfdt=0.0
#dgdt=0.0
#A3=0.0
#B3=0.0
#Q3=0.0
#R3=0.0
#iter =0
#ds = 0.0
#x=zeros(3)
#v=zeros(3)
include("kepler_solver_wh4.jl")
#include("kepler_solver_wh4_kepstep.jl")

#function test_kepler!(h::Float64,nsteps::Int64,x0::Vector{Float64},v0::Vector{Float64},k::Float64,ks::Kepler_step)
function test_kepler!(h::Float64,nsteps::Int64,x0::Vector{Float64},v0::Vector{Float64},k::Float64)
for j=1:3
  xsave[j,1]=x0[j]
  vsave[j,1]=v0[j]
end
ssave[1]=0.0

#h = 10.0*24.0*3600.0 # 10-day timesteps

# Now, loop over steps:
s0 = 0.0
tsave[1]=0.0
x = zeros(3)
v = zeros(3)
s = 0.0


@inbounds for i=2:nsteps
  tsave[i]=tsave[i-1]+h
#  x,v,s = kep_eqn!(x0,v0,k,h,s0,ks)
  x,v,s = kep_eqn!(x0,v0,k,h,s0)
#  println("s0: ",s0," s: ",s," diff: ",s/s0-1)
  ssave[i]=s
  for j=1:3
    xsave[j,i] =x[j]
    vsave[j,i] =v[j]
  end
  x0 = x
  v0 = v 
  if i >= 4
    s0 = 3.0*s - 3.0*ssave[i-1] + ssave[i-2]
  end
#  println(i,x,v)
end

r = norm(xsave[:,nsteps])
beta_final = 2.0*k/r-dot(vsave[:,nsteps],vsave[:,nsteps])

#using PyPlot
#plot(tsave,xsave[1,:])
#plot(tsave,xsave[2,:])
#return tsave,xsave,vsave
#return tsave,xsave,vsave
return beta_final
end



m1 = 1.0
m2 = MEARTH/MSUN
m = 1./(1./m1+1./m2)
#k = GRAV*m1*m2/m
# v = 2piAU/YR; v^2*AU = k
k = (2pi/365.25)^2*(m1+m2)
#x0 = [0,AU,0]
x0 = [0,1.0,0]
r0 = norm(x0)
#v0 = [30e5,0,0]
#v0 = [1.4*2pi/365.25,0,0]
#v0 = [sqrt(2.0*k/r0),0,0]
v0 = [1.1*sqrt(2.0*k/r0),0,0]
h = 1.0 # 5-day timesteps

#@time tsave,xsave,vsave = test_kepler!(h,nsteps,x0,v0,k)
s0 = 0.0
#kepstep=Kepler_step()
#beta_final = @time test_kepler!(h,nsteps,x0,v0,k,kepstep)
#beta_final = @time test_kepler!(h,nsteps,x0,v0,k,kepstep)
beta_final = @time test_kepler!(h,nsteps,x0,v0,k)
beta_final = @time test_kepler!(h,nsteps,x0,v0,k)

plot(tsave,xsave[1,:])
plot(tsave,xsave[2,:])

end
