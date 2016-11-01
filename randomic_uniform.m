%|=====================================|
%|Ash Richardson, BSc, MSc Candidate   |
%|Graduate Student Member, IEEE        |
%|                                     |
%|Mathematics and Statistics           |
%|University of Victoria               |
%|PO BOX 3060 STN CSC                  |
%|Victoria, B.C. / Canada V8W 3R4      |
%|                                     |
%|mailto:ashy@uvic.ca                  |
%|+12508533292 / SSM A518              |
%|=====================================|

%generate random initial conditions for sebastian picard's implementation
%of the "refined flocking and swarming model".
%
%parameters:    n - the number of birds.
%               d_cutoff - distance cutoff parameter
%               fric_alpha - friction parameter
%               fric_beta - friction parameter
function B = randomic_uniform(n)
    par = getparameters();
    d_cutoff = par.d_cutoff;
    fric_alpha = par.fric_alpha;
    fric_beta = par.fric_beta;
% Each column (x,xdot,y,ydot,z,zdot) represents one bird.
    B = rand(6,n);%1+rand(6,n);
    c1 = 200;%10*d_cutoff/sqrt(3.);
    c2 = (1/3.)*sqrt(fric_alpha/fric_beta)/sqrt(3.);
    
    B(1,:)=B(1,:)*c1; B(3,:)=B(3,:)*c1; B(5,:)=B(5,:)*(c1/5);
    B(2,:)=B(2,:)-0.5; B(4,:)=B(4,:)-0.5; B(6,:)=B(6,:)-0.5;
    B(2,:)=B(2,:)*c2; B(4,:)=B(4,:)*c2; B(6,:)=B(6,:)*(c2/12.);

