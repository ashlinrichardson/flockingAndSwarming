%May 25, 2010 - Added new force term(according to Reinhard).
% Note that initial conditions are passed in the format of a bird matrix
% B like the one described below. However, for some reason MATLAB
% converts this matrix into one very long column vector when it
% runs ode45. The first few lines of the code convert
% this column vector back into the matrix B so that it is easier to
% work with.

function dA = n_interac_3D(t,a,p)

% Set constants
d_cutoff    = p.d_cutoff;%10; %0.01
e_cutoff    = p.e_cutoff; 
v_cutoff    = p.v_cutoff;%2; % 2
ve_cutoff   = p.ve_cutoff;%2; % 2
R_coeff     = p.R_coeff;%0.2;
S_coeff     = p.S_coeff;%2;
angle_cut1  = p.angle_cut1;%1.57;  % 80deg : delta_1
angle_cut2  = p.angle_cut2;%1.047; % 60deg: delta_2
alpha       = p.alpha;%0.5;
T_coeff     = p.T_coeff;%8;
%C_coeff     = p.C_coeff;%100;
fric_alpha  = p.fric_alpha;%1;
fric_beta   = p.fric_beta;%0.01;
B_delta = p.B_delta;
B_epsilon = p.B_epsilon;
%B_lambda    = p.B_lambda;

% Get number of birds
global n;

% Now convert input vector a into bird matrix B. This is only to make
% things easier to read. Each bird has its coordinates
% (x,xdot,y,ydot,z,zdot) stored in a column of the matrix. Every
% column represents a different bird.

% Create bird matrix and derivative of bird matrix
B = zeros(6,n);
dB = zeros(6,n);


% Fill bird matrix with IC input vector
for j=1:(6*n), 
  B(mod(j-1,6)+1,floor((j-1)/6)+1) = a(j);
end

% Now loop to calculation motion for each ith bird
for i=1:n,
  
  % Calculate bird speed
  V = (B(2,i)*B(2,i)+B(4,i)*B(4,i)+B(6,i)*B(6,i)).^(1/2);
  if (V==0)
    V = 1e-100; % Quick fix to avoid division by zero
  end
  
  fricV = (fric_alpha - fric_beta*V*V);
  
  if(V <= v_cutoff) 
      s1 = 1;
  elseif(V >= v_cutoff + ve_cutoff) 
      s1 = 0;
  else
      Vt = (V-v_cutoff)/ve_cutoff;
      s1 = .5*(1.+tanh( (1./(abs(Vt)-0)) +  (1./(abs(Vt)-1)) ));
  end
 
  ns1 = 1-s1;
  % Set interactions at zero before getting contributions from
  % other birds
  ix = 0;
  iy = 0;
  iz = 0;
  
  
   
  % Loop over all other birds for interactions
  for j=1:n,
    
    % Calculate separation between bird i and bird j
    Rij = ( (B(1,i) - B(1,j)).^2 + (B(3,i) - B(3,j)).^2 + (B(5,i) - B(5,j)).^2 ).^(1/2);
    Rij = max( Rij, eps);

    if (Rij==0)
        Rij = eps; % Quick fix to avoid division by zero
    end

    if(Rij <= d_cutoff)
        s0 = 1;
    elseif(Rij >= d_cutoff + e_cutoff) 
        s0 = 0;
    else
        %scale input from 0 to 1
        Rijt = (Rij-d_cutoff)/e_cutoff;   
        s0 = .5*(1.+tanh( (1./(abs(Rijt)-0)) +  (1./(abs(Rijt)-1)) ));
    end

    costheta = (B(2,i)*(B(1,j)-B(1,i)) + B(4,i)*(B(3,j)-B(3,i)) + B(6,i)*(B(5,j)-B(5,i)))/((Rij)*(V));	

    if(costheta <= cos(angle_cut1))     
        s2 = 0;
    elseif(costheta >= cos(angle_cut2)) 
        s2 = 1;
    else
        %scale input from 0 to 1
        At = (costheta-cos(angle_cut1))/(cos(angle_cut2) - cos(angle_cut1));
        s2 = 1-(.5*(1.+tanh( (1./(abs(At)-0)) +  (1./(abs(At)-1)) )));
    end
%==========================================================================
%repulsive force R_i
    ix = ix + s0*(R_coeff)*(1/(Rij))*(B(1,i)-B(1,j))/(1+Rij.^2).^(alpha);
    iy = iy + s0*(R_coeff)*(1/(Rij))*(B(3,i)-B(3,j))/(1+Rij.^2).^(alpha);
    iz = iz + s0*(R_coeff)*(1/(Rij))*(B(5,i)-B(5,j))/(1+Rij.^2).^(alpha);
%==========================================================================
%attractive force A_i
    ix = ix + (s1+ns1*s2)*(S_coeff)*(B(2,j)-B(2,i))/(1+Rij.^2).^(alpha);
    iy = iy + (s1+ns1*s2)*(S_coeff)*(B(4,j)-B(4,i))/(1+Rij.^2).^(alpha); 
    iz = iz + (s1+ns1*s2)*(S_coeff)*(B(6,j)-B(6,i))/(1+Rij.^2).^(alpha); 
%==========================================================================

    
  end % summation of interactions loop (j loop)
  
  
  % Calculate xdotdot,ydotdot,zdotdot for bird
  dB(2,i)=ix/n; dB(4,i)=iy/n;  dB(6,i)=iz/n;
    
  % Get the density of birds around bird i
  bird_density = 0;
 
  
  for q=1:n,
    if (q ~= i),
      Riq = ( (B(1,i) - B(1,q)).^2 + (B(3,i)-B(3,q)).^2 + (B(5,i)-B(5,q)).^2 ).^(1/2);     
      bird_density = bird_density + (1./(1+abs((Riq)*abs(Riq))));
    end
  end
  bird_density = bird_density/(1.*n);
  sf =  bird_density;
  bd(i)=bird_density;
  s0=0;
  if(sf <= B_delta)
        s0 = 1;
  elseif(sf >= B_delta + B_epsilon) 
      	s0 = 0;
  else
            %scale input from 0 to 1
            sft = (sf-B_delta)/B_epsilon;   
            s0 = .5*(1.+tanh( (1./(abs(sft)-0)) +  (1./(abs(sft)-1)) ));
  end
  
 %DELETE THIS.......
  if(sf <= B_delta)
        s0 = 1;
  else
        s0 =0;
  end
  
  % Get turning rho-factor
  rho_factor = s0;%1/(1 + C_coeff*bird_density);
  
%==========================================================================
  % Turning force B_i
  T_force_x = (1.)*(rho_factor)*T_coeff*B(4,i);
  T_force_y = (-1.)*(rho_factor)*T_coeff*B(2,i);
  if( mod(i,2)==0  )
      T_force_x = T_force_x*(-1);
      T_force_y = T_force_y*(-1);

  end
  
%==========================================================================
  
  % Add turning force to interact contributions
  dB(2,i) = dB(2,i) + T_force_x;
  dB(4,i) = dB(4,i) + T_force_y;
  
  % Add friction term
  dB(2,i) = dB(2,i) + fricV*B(2,i);
  dB(4,i) = dB(4,i) + fricV*B(4,i);
  dB(6,i) = dB(6,i) + fricV*B(6,i);
  
 
end % loop over all birds  (i loop).
%[min(bd), max(bd), sum(bd)/n]


% Fill rest of derivative of bird matrix (xdot,ydot and zdot)
for k=1:n,
  dB(1,k)=B(2,k);
  dB(3,k)=B(4,k);
  dB(5,k)=B(6,k);
end


% Now must convert back to vector format to return
dA = zeros(6*n,1);

for j=1:(6*n),
  dA(j) = dB(mod(j-1,6)+1,floor((j-1)/6)+1);
end

