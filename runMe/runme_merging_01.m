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
function retval = runme_merging_01(n)
%======================================================================
%PARAMETERS
%======================================================================
    par.d_cutoff   = 1; %default 50 %0.5;%0.7; %0.01;%10; %0.01
    par.e_cutoff   = 1;%par.d_cutoff/10.; %0.5;%0.7; %0.01;%10; %0.01
    
    par.R_coeff    = 10.;%1.;% 50;%1.;%1; %10 %0.2   % default 1
    par.S_coeff    = 10;%5.;% 50;   %default 100 
    par.angle_cut1 = 1.57;  % 80deg
    par.angle_cut2 = 1.047; % 60deg    %CHECK NOTATION.... 
    par.alpha      = 0.5;
    par.T_coeff    =0.003;%05;% 0.0001;%0.001;%50;%0.05;  
    par.fric_alpha = 1;
    par.fric_beta  = 1;
    
    par.v_cutoff   = 0.5;%0.99*sqrt(par.fric_alpha/par.fric_beta);%0.5; % 2
    par.ve_cutoff  = 0.5;%par.v_cutoff/12.;
    
    par.B_delta    = 0.0015;%000005;%par.d_cutoff/0.5;%1/20;    
    par.B_epsilon  = 0.0003;%par.d_cutoff/0.5;%1/20;
%======================================================================
%INITIAL DATA
%======================================================================
%    d_cutoff = par.d_cutoff;
% Each column (x,xdot,y,ydot,z,zdot) represents one bird.

    A = 0 + rand(6,n);%1+rand(6,n);
    c1 = 200;%10*d_cutoff/sqrt(3.);
    c2 =sqrt(par.fric_alpha/par.fric_beta);
    c2 = c2/2;
    A(1,:)=A(1,:)*c1; A(3,:)=(4*c1)+A(3,:)*c1; A(5,:)=A(5,:)*(c1);
    A(2,:)=A(2,:)-0.5; A(4,:)=A(4,:)-0.5; A(6,:)=A(6,:)-0.5;
    A(2,:)=A(2,:)*c2; A(4,:)=A(4,:)*c2; A(6,:)=A(6,:)*(c2/5.);
    retval =0;
    
    
    B = 0 + rand(6,n);%1+rand(6,n);
    c1 = 200;%10*d_cutoff/sqrt(3.);
    c2 =sqrt(par.fric_alpha/par.fric_beta);% 
    
    B(1,:)=B(1,:)*c1/0.7; B(3,:)=B(3,:)*c1/0.7; B(5,:)=B(5,:)*(c1/5)/0.7;
    
    B(2,:)=0; B(4,:)=10; B(6,:)=0;
    retval =0;
    adims = size(A)
    C= cat( 2, A, B);
    csize = size(C) 
    asize = size(A);
    
%======================================================================
%RUN THE SIMULATION
%======================================================================
    initial_conditions = C;
    n=n*2;
    T = 1000;
    save runme_merging_last_data.mat;
    n_interac_3D_mov( initial_conditions, par, T, T/5);  
