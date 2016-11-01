% Note that initial conditions are passed in the format of a bird matrix
% B like the one described below. However, for some reason MATLAB
% converts this matrix into one very long column vector when it
% runs ode45. The first few lines of the code convert
% this column vector back into the matrix B so that it is easier to
% work with.

function dA = n_interac_3D(t,a,p)

% Set constants
d_cutoff    = p.d_cutoff;%10; %0.01
v_cutoff    = p.v_cutoff;%2; % 2
R_coeff     = p.R_coeff;%0.2;
S_coeff     = p.S_coeff;%2;
angle_cut1  = p.angle_cut1;%1.57;  % 80deg
angle_cut2  = p.angle_cut2;%1.047; % 60deg
alpha       = p.alpha;%0.5;
T_coeff     = p.T_coeff;%8;
C_coeff     = p.C_coeff;%100;
fric_alpha  = p.fric_alpha;%1;
fric_beta   = p.fric_beta;%0.01;

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
  
  % Set interactions at zero before getting contributions from
  % other birds
  sum_interact_x = 0;
  sum_interact_y = 0;
  sum_interact_z = 0;
   
  % Loop over all other birds for interactions
  for j=1:n,
    
    % Calculate separation between bird i and bird j
    Rij = ( (B(1,i) - B(1,j)).^2 + (B(3,i) - B(3,j)).^2 + (B(5,i) - B(5,j)).^2 ).^(1/2);
    if (Rij==0)
      Rij = 1e-100; % Quick fix to avoid division by zero
    end
    
    % If bird i and j are too close, use repulsive acceleration
    if (Rij < d_cutoff)
      interact_x = (R_coeff)*(1/(Rij))*(B(1,i)-B(1,j))/(1+Rij.^2).^(alpha);
      interact_y = (R_coeff)*(1/(Rij))*(B(3,i)-B(3,j))/(1+Rij.^2).^(alpha);
      interact_z = (R_coeff)*(1/(Rij))*(B(5,i)-B(5,j))/(1+Rij.^2).^(alpha);
    
    % If bird i and j are far enough, use attractive acceleration
    else
      
      % If bird is moving slowly, it will orient itself
      if (V < v_cutoff)     
        interact_x = (S_coeff)*(B(2,j)-B(2,i))/(1+Rij.^2).^(alpha);
        interact_y = (S_coeff)*(B(4,j)-B(4,i))/(1+Rij.^2).^(alpha); 
        interact_z = (S_coeff)*(B(6,j)-B(6,i))/(1+Rij.^2).^(alpha); 
	
	% If bird is moving quickly, will have peripheral vision
      else
	% Calculate cosine of angle between bird i velocity and
        % vector Rij	
	costheta = (B(2,i)*(B(1,j)-B(1,i)) + B(4,i)*(B(3,j)-B(3,i)) ...
		    + B(6,i)*(B(5,j)-B(5,i)))/((Rij)*(V));	
	
	% If bird j is in field of vision of bird i
	if (costheta>=cos(angle_cut2))
	  interact_x = (S_coeff)*(B(2,j)-B(2,i))/(1+Rij.^2).^(alpha);
	  interact_y = (S_coeff)*(B(4,j)-B(4,i))/(1+Rij.^2).^(alpha); 
	  interact_z = (S_coeff)*(B(6,j)-B(6,i))/(1+Rij.^2).^(alpha); 
	
	% If bird j is in peripheral vision of bird i
	elseif ((cos(angle_cut1)<costheta) && (costheta<cos(angle_cut2)))
	  decreasing_factor = (costheta-cos(angle_cut1))/(cos(angle_cut2)-cos(angle_cut1));
	  interact_x = (S_coeff)*(decreasing_factor)*(B(2,j)-B(2,i))/(1+Rij.^2).^(alpha);
	  interact_y = (S_coeff)*(decreasing_factor)*(B(4,j)-B(4,i))/(1+Rij.^2).^(alpha); 
	  interact_z = (S_coeff)*(decreasing_factor)*(B(6,j)-B(6,i))/(1+Rij.^2).^(alpha); 

	% If bird j is outside of cone of vision of bird i
	else	  
	  interact_x=0;
	  interact_y=0;
	  interact_z=0;	 	  	  
	end
      end
    end   
         	  
    % Add interaction contribution to sum
    sum_interact_x = sum_interact_x + interact_x;
    sum_interact_y = sum_interact_y + interact_y;
    sum_interact_z = sum_interact_z + interact_z;
    
  end % summation of interactions loop
  
  
  % Calculate xdotdot,ydotdot,zdotdot for bird
  dB(2,i)=sum_interact_x/n;
  dB(4,i)=sum_interact_y/n;
  dB(6,i)=sum_interact_z/n;
    
  % Get the density of birds around bird i
  bird_density = 0;
  
  for q=1:n,
    if (q ~= i),
      Riq = ( (B(1,i) - B(1,q)).^2 + (B(3,i)-B(3,q)).^2 + (B(5,i)-B(5,q)).^2 ).^(1/2);     
      if (Riq==0)
	Riq=0.000001;
      end    
      bird_density = bird_density + 1/(Riq);      
    end
  end
  
  % Get turning rho-factor
  rho_factor = 1/(1 + C_coeff*bird_density);
  
  % Define turning force  
  T_force_x = (1)*(rho_factor)*T_coeff*B(4,i);
  T_force_y = (-1)*(rho_factor)*T_coeff*B(2,i);
  
  % Add turning force to interact contributions
  dB(2,i) = dB(2,i) + T_force_x;
  dB(4,i) = dB(4,i) + T_force_y;
  
  % Add friction term
  dB(2,i) = dB(2,i) + (fric_alpha - fric_beta*V*V)*B(2,i);
  dB(4,i) = dB(4,i) + (fric_alpha - fric_beta*V*V)*B(4,i);
  dB(6,i) = dB(6,i) + (fric_alpha - fric_beta*V*V)*B(6,i);
  
 
end % loop over all birds


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

%B
%dB

%%%%%%%%%%%%%%%%%%%%% rho %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Old code that may be useful. This counts the number of other
% birds that the bird i can see in its peripheral vision.
  
% Check if bird sees anyone
% birdcount = 0;
% 
% if (V<0.1)
%   birdcount = n;
%   disp('slow vel is slow');
%   
% else  
%   for q=1:n,
%     if (q ~= i),
%	Riq = ( (B(1,i) - B(1,q)).^2 + (B(3,i) - B(3,q)).^2 ).^(1/2);
%      
%	if (Riq==0)
%	  costheta = 0;
%	else
%	  costheta = (B(2,i)*(B(1,q)-B(1,i)) + B(4,i)*(B(3,q)-B(3, ...
%						  i)))/((Riq)*(V));
%	end
%
%	if (costheta>=cos(angle_cut2))
%	  birdcount = birdcount + 1;
%	end
%      
%      end      
%    end
%  end
%  %%%%%%%%%%%%%%%%%%%%%% rho %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
