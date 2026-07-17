function grd = kellogg_grad_u_exact(x)
% function grd = kellogg_grad_u_exact(x)



pi_2 = pi/2;
pi_32 = 3*pi/2;

gamma = 0.1;
rho = 0.25*pi;
sigma  = -14.92256510455152;
cosR1 = cos((pi_2 - sigma)*gamma);
cosR2 = cos(rho*gamma);
cosR3 = cos(sigma*gamma);
cosR4 = cos((pi_2 - rho)*gamma);

r = norm(x);
phi = atan2(x(2),x(1))+(2*pi*(x(2)<0));

grd = 0*x;
if (r < 1.e-10) 
  return
end

ara2 = gamma*r^(gamma - 2.0);


  if (x(1)>=0 && x(2)>=0)
    c = cos((phi - pi_2 + rho)*gamma);
    s = sin((phi - pi_2 + rho)*gamma);
    grd(1) = ara2*cosR1*(x(1)*c + x(2)*s);
    grd(2) = ara2*cosR1*(x(2)*c - x(1)*s);
  elseif (x(1)<=0 && x(2)>=0)
    c = cos((phi - pi + sigma)*gamma);
    s = sin((phi - pi + sigma)*gamma);
    grd(1) = ara2*cosR2*(x(1)*c + x(2)*s);
    grd(2) = ara2*cosR2*(x(2)*c - x(1)*s);
  elseif (x(1)<=0 && x(2)<=0)
    c = cos((phi - pi - rho)*gamma);
    s = sin((phi - pi - rho)*gamma);
    grd(1) = ara2*cosR3*(x(1)*c + x(2)*s);
    grd(2) = ara2*cosR3*(x(2)*c - x(1)*s);
  else
    c = cos((phi - pi_32 - sigma)*gamma);
    s = sin((phi - pi_32 - sigma)*gamma);
    grd(1) = ara2*cosR4*(x(1)*c + x(2)*s);
    grd(2) = ara2*cosR4*(x(2)*c - x(1)*s);
  end

end