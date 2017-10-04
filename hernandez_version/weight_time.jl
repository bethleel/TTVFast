function weight_time(x)
@assert(x >= 0.0)
@assert(x <= 1.0)
if x < 0.5
  return 1.0-2.*x^2,-4.*x
else
  return 2.*(1.0-x)^2,-4.*(1.0-x)
end
end

nx = 1000
x = linspace(0.0,1.0,nx)
w = zeros(nx)
wr = zeros(nx)
dwdx = zeros(nx)
for i=1:nx
  w[i],dwdx[i] = weight_time(x[i])
end

using PyPlot
plot(x,w)
plot(x,cos(x*pi/2).^2)
plot(x,1.-w)
plot(x,sin(x*pi/2).^2)
plot(x,dwdx)
plot(x,-2*(cos(x*pi/2).*sin(x*pi/2)*pi/2))
