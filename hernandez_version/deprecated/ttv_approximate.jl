using PyPlot
include("ttv.jl")
include("/Users/ericagol/Computer/Julia/regress.jl")

n = 8
t0 = 7257.93115525
h  = 0.05
tmax = 600.0

# Read in initial conditions:
elements0 = readdlm("elements.txt",',')

# Make an array, tt, to hold transit times:
# First, though, make sure it is large enough:
ntt = zeros(Int64,n)
for i=2:n
  ntt[i] = ceil(Int64,tmax/elements0[i,2])+2
end
println("ntt: ",ntt)
tt1 = zeros(n,maximum(ntt))
# Save a counter for the actual number of transit times of each planet:
count1 = zeros(Int64,n)
# Call the ttv function:
@time ttv!(n,t0,h,tmax,elements0,tt1,count1)
# Now call with half the timestep:
count2 = zeros(Int64,n)

# Now, pertub the elements, keep track of the results, and see if
# we can "learn" how to model the transit times based on the initial conditions:

# Number of cases to vary:
nvary = 64
elements_vary = zeros(n,5,nvary)
tt1_vary = zeros(n,maximum(ntt))
tt_vary = zeros(n,maximum(ntt),nvary)
mag_vary = zeros(n,5)
psig = [0.,1e-6,4e-6,1e-3,2e-5,3e-5,2e-4,8e-3]
t0sig = [0.,2e-4,4e-4,7e-4,8e-4,5e-4,4e-4,8e-4]
ecc = [0.,.02,.02,.015,.01,.015,.005,.02]
for i=2:n
  mag_vary[i,1] = 0.5*elements0[i,1]
  mag_vary[i,2] = psig[i]
  mag_vary[i,3] = t0sig[i]
  mag_vary[i,4] = ecc[i]
  mag_vary[i,5] = ecc[i]
end
count10 = copy(count1)
elements = copy(elements0)
for ivary=1:nvary
  fill!(count1,0)
  dcount = maximum(abs(count1-count10))
  while dcount > 0
    fill!(count1,0)
    min_mass = -1.0
    while min_mass < 0.0
      elements .= elements0 .+ randn(n,5).*mag_vary
      min_mass = Inf
      for k=2:n
        if elements[k,1] < min_mass
          min_mass = elements[k,1]
        end
      end
    end
    ttv!(n,t0,h,tmax,elements,tt1_vary,count1)
    dcount = maximum(abs(count1-count10))
  end
  tt_vary[:,:,ivary] = tt1_vary
  elements_vary[:,:,ivary] = elements
  println(ivary," ",elements)
end

using PyPlot

# Make a plot of some TTVs:

fig,axes = subplots(4,2)

ttv_vary = copy(tt_vary)

for i=2:8
  ax = axes[i-1]
  fn = zeros(Float64,2,count1[i])
  sig = ones(count1[i])
  tti1 = tt1[i,1:count1[i]]
  for j=1:count1[i]
    fn[1,j] = 1.0
    fn[2,j] = round(Int64,(tti1[j]-elements0[i,3])/elements0[i,2])
  end
  coeff,cov = regress(fn,tti1,sig)
  tt_ref1 = coeff[1]+coeff[2]*fn[2,:]
  ttv1 = (tti1-tt_ref1)*24.*60.
  ax[:plot](tti1,ttv1)
#  ax[:plot](tti2,ttv2)
#  ax[:plot](tti2,((ttv1-ttv2)-mean(ttv1-ttv2))*60.)
# Plot the varied TTVs:
  for ivary=1:nvary
    tti2 = tt_vary[i,1:count1[i],ivary]
    coeff,cov = regress(fn,tti2,sig)
    tt_ref2 = coeff[1]+coeff[2]*fn[2,:]
    ttv2 = (tti2-tt_ref2)*24.*60.
    ax[:plot](tti2,ttv2)
    ttv_vary[i,1:count1[i],ivary] = ttv2
  end
#  println(i," ",coeff," ",elements0[i,2:3]," ",coeff[1]-elements0[i,3]," ",coeff[2]-elements0[i,2])
#  println(i," ",maximum(ttv1-ttv2-mean(ttv1-ttv2))*60.," ", minimum(ttv1-ttv2-mean(ttv1-ttv2))*60.)
end

read(STDIN,Char)

nparam = (n-1)*5
# Set up grid of "functions" to regress against:
fn = zeros(nparam,nvary)
sig = ones(nvary)
# Now, see if we can make sense of these.
#for i=1:n
for i=2:8
# Loop over the transit times for the ith planet:
  for j=1:count1[i]
#  for j=1:1
  # Carry out a regression against the initial parameters:
    for ivary=1:nvary
      fn[1:nparam,ivary] = vcat(elements_vary[2:n,:,ivary]...)
#      fn[1,ivary] = tt_vary[i,1,ivary]
#      fn[2,ivary] = tt_vary[i,count1[i],ivary]
#      fn[3,ivary] = tt_vary[i,floor(Int64,count1[i]/2),ivary]
    end
#    coeff,cov = regress(fn,tt_vary[i,j,:],sig)
    coeff,cov = regress(fn,ttv_vary[i,j,:],sig)
  # Now, compute residuals:
    dt = zeros(nvary)
    tti = zeros(nvary)
    ttv_pred = zeros(nvary)
    for ivary=1:nvary
#      tt_pred[ivary] = coeff[1]*tt_vary[i,        1,ivary]
#      tt_pred[ivary] += coeff[2]*tt_vary[i,count1[i],ivary]
#      tt_pred[ivary] += coeff[3]*tt_vary[i,floor(Int64,count1[i]/2),ivary]
      param =  vcat(elements_vary[2:n,:,ivary]...)
      for k=1:nparam
        ttv_pred[ivary] += coeff[k]*param[k]
      end
      dt[ivary] = ttv_vary[i,j,ivary]-ttv_pred[ivary]
      tti[ivary] = ttv_vary[i,j,ivary]
    end
    println(i," ",j," ",mean(dt)," ",std(dt)," ",median(abs(dt))," ",std(tti))
    plot(tti,ttv_pred,".")
  end
end

# The following version regressed against the first & last transit times:
#nparam = (n-1)*5+3
## Set up grid of "functions" to regress against:
#fn = zeros(nparam,nvary)
#sig = ones(nvary)
## Now, see if we can make sense of these.  
##for i=1:n
#for i=7:7
## Loop over the transit times for the ith planet:
#  for j=1:count1[i]
##  for j=1:1
#  # Carry out a regression against the initial parameters:
#    for ivary=1:nvary
#      fn[4:nparam,ivary] = vcat(elements_vary[2:n,:,ivary]...)
#      fn[1,ivary] = tt_vary[i,1,ivary]
#      fn[2,ivary] = tt_vary[i,count1[i],ivary]
#      fn[3,ivary] = tt_vary[i,floor(Int64,count1[i]/2),ivary]
#    end
#    coeff,cov = regress(fn,tt_vary[i,j,:],sig)
#  # Now, compute residuals:
#    dt = zeros(nvary)
#    tti = zeros(nvary)
#    tt_pred = zeros(nvary)
#    for ivary=1:nvary
#      tt_pred[ivary] = coeff[1]*tt_vary[i,        1,ivary]
#      tt_pred[ivary] += coeff[2]*tt_vary[i,count1[i],ivary]
#      tt_pred[ivary] += coeff[3]*tt_vary[i,floor(Int64,count1[i]/2),ivary]
#      param =  vcat(elements_vary[2:n,:,ivary]...)
#      for k=4:nparam
#        tt_pred[ivary] += coeff[k]*param[k-3]
#      end
#      dt[ivary] = tt_vary[i,j,ivary]-tt_pred[ivary]
#      tti[ivary] = tt_vary[i,j,ivary]
#    end
#    println(i," ",j," ",mean(dt)," ",std(dt)," ",median(abs(dt))," ",std(tti))
#    plot(tti,tt_pred,".")
#  end
#end

