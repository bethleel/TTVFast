const YEAR  = 365.242
const GNEWT = 39.4845/YEAR^2  # Units of MSUN*AU^3/YEAR^2
const NDIM  = 3
const third = 1./3.

include("init_nbody.jl")

elements = readdlm("elements.txt",',')

n_body = 8
t0 = 7257.93115525
jac_init = zeros(Float64,7*n_body,7*n_body)
x,v = init_nbody(elements,t0,n_body,jac_init)

println(x)
println(v)
println(jac_init)

# Estimate orbital periods of planets about the star:
mass = elements[1,1]
for i=2:n_body
  mass += elements[i,1]
  ri = sqrt((x[1,i]-x[1,1])^2 + (x[2,i]-x[2,1])^2 + (x[3,i]-x[3,1])^2)
  vi2 = (v[1,i]-v[1,1])^2 + (v[2,i]-v[2,1])^2 + (v[3,i]-v[3,1])^2
  gm = GNEWT*mass
  semi = gm*ri/(2.*gm-vi2*ri)
#  println(i, " ", semi)
  period = 2pi*sqrt(semi^3/gm)
  println(i, " ", period, " ",elements[i,2])
end
