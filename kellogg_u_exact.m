function u = kellogg_u_exact(x)
% function u = kellogg_u_exact(x)

pi_2 = pi/2;
pi_32 = 3*pi/2;

gamma = 0.1;
rho = 0.25*pi;
sigma  = -14.92256510455152;
% gamma3 = gamma1 = 161.4476387975881;
% gamma4 = gamma2 = 1.0;
cosR1 = cos((pi_2 - sigma)*gamma);
cosR2 = cos(rho*gamma);
cosR3 = cos(sigma*gamma);
cosR4 = cos((pi_2 - rho)*gamma);

r = norm(x);
phi = atan2(x(2),x(1))+(2*pi*(x(2)<0));

if (r < 1.e-10) 
  u = 0;
  return
end

  if (x(1)>=0 && x(2)>=0)
    mu = cosR1*cos((phi - pi_2 + rho)*gamma);
  elseif (x(1)<=0 && x(2)>=0)
    mu = cosR2*cos((phi - pi + sigma)*gamma);
  elseif (x(1)<=0 && x(2)<=0)
    mu = cosR3*cos((phi - pi - rho)*gamma);
  else
    mu = cosR4*cos((phi - pi_32 - sigma)*gamma);
  end

  u = r^gamma * mu;
end