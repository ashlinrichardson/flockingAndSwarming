% Note that initial conditions are passed in the format of a bird matrix
% B like the one described below. However, for some reason MATLAB
% converts this matrix into one very long column vector when it
% runs ode45. The first few lines of the code convert
% this column vector back into the matrix B so that it is easier to
% work with.
function dA = n_interac_3D_1(t,a,p)

% Set constants
d_cutoff    = p.d_cutoff;%10; %0.01           %d_0
v_cutoff    = p.v_cutoff;%2; % 2              %d_1
R_coeff     = p.R_coeff;%0.2;                 %rho_1
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
%B = zeros(6,n);
dB = zeros(6,n);

% Fill bird matrix with IC input vector
%for j=1:(6*n), 
%  B(mod(j-1,6)+1,floor((j-1)/6)+1) = a(j);
%end
B = reshape( a, 6, n);
global X_i; global X_j; global V_i; global V_j; global interact; global sum_interact;

X_i = [0 0 0];
X_j = [0 0 0];
 V_i =[0 0 0];
V_j = [0 0 0];
 interact = [0. 0. 0.]; 
sum_interact = [0. 0. 0.];

% Now loop to calculation motion for each ith bird
for i=1:n,
  X_i = [B(1,i) B(3,i) B(5,i)];
  % Calculate bird speed
  V_i = [ B(2,i) B(4,i) B(6,i)];  
  V = norm(V_i)+eps; %max( norm(V_i), eps);  %avoid division by zero.
  
  % Set interactions at zero before getting contributions from
  % other birds

  interact = [0. 0. 0.]; 
  sum_interact = [0. 0. 0.];
   
  % Loop over all other birds for interactions
  for j=1:n,
    %velocity of bird j
    V_j = [B(2,i) B(4,i) B(6,i)]; 
    % Calculate separation between bird i and bird j
    X_j = [B(1,j) B(3,j) B(5,j)];
    Rij = norm( X_j - X_i) + eps;%max( norm( X_j - X_i), eps);  %avoid division by 0
    
    %============================================================  
    % If bird i and j are too close, use repulsive acceleration
    if (Rij < d_cutoff)
      interact = (R_coeff)*(1/(Rij)).* (X_i - X_j) ./ ((1+(Rij^2))^alpha);
    
    % If bird i and j are far enough, use attractive acceleration
    else
      
      % If bird is moving slowly, it will orient itself
      if (V < v_cutoff)  
        interact = (S_coeff).*(V_j - V_i)./( (1+Rij^2.)^alpha);
	% If bird is moving quickly, will have peripheral vision
      else
	% Calculate cosine of angle between bird i velocity and
        % vector Rij	
            %costheta = (B(2,i)*(B(1,j)-B(1,i)) + B(4,i)*(B(3,j)-B(3,i))+ B(6,i)*(B(5,j)-B(5,i)))/((Rij)*(V));	
            costheta = sum(V_j .* ( X_j - X_i)) / ( Rij*V);
	
	% If bird j is in field of vision of bird i
        if (costheta>=cos(angle_cut2))
          interact = (S_coeff).*(V_j - V_i)./((1+(Rij.^2)).^(alpha));

        % If bird j is in peripheral vision of bird i
        elseif ((cos(angle_cut1)<costheta) && (costheta<cos(angle_cut2)))
          decreasing_factor = (costheta-cos(angle_cut1))/(cos(angle_cut2)-cos(angle_cut1));
          interact = (S_coeff)*(decreasing_factor).*(V_j - V_i)./((1+(Rij.^2)).^(alpha));

        % If bird j is outside of cone of vision of bird i
        else
           interact = [0. 0. 0.];
        end
      end
    end %============================================================  
         	  
    % Add interaction contribution to sum
    sum_interact = sum_interact + interact;
    
  end % (j-loop) summation of interactions loop
  
  
  % Calculate xdotdot,ydotdot,zdotdot for bird
   dB(2,i)=sum_interact(1)/n;
   dB(4,i)=sum_interact(2)/n;
   dB(6,i)=sum_interact(3)/n;
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
%dA = zeros(6*n,1);

%for j=1:(6*n),
%  dA(j) = dB(mod(j-1,6)+1,floor((j-1)/6)+1);
%end
dA = reshape(dB, prod(size(dB)),1);
%B
%dB

