%parameter file for May 25 version (new Force term - Reinhard).
function par = getparameters()

noturning = 1;
onlyturning = 0;

if( noturning==1)
    par.d_cutoff   = 1; %default 50 %0.5;%0.7; %0.01;%10; %0.01
    par.e_cutoff   = par.d_cutoff/10.; %0.5;%0.7; %0.01;%10; %0.01
    par.v_cutoff   = 1;%0.5; % 2
    par.ve_cutoff  = 1;
    par.R_coeff    = 1;%1.;%1; %10 %0.2   % default 1
    par.S_coeff    = 1;   %default 100 
    par.angle_cut1 = 1.57;  % 80deg
    par.angle_cut2 = 1.047; % 60deg    %CHECK NOTATION.... 
    par.alpha      = 0.5;
    par.T_coeff    = 0;%1;%0.05;  %default 0.05 set to 1 for overcurly.
    par.C_coeff    = 1;
    par.fric_alpha = 1;
    par.fric_beta  = 2;%0.01;
    par.B_lambda  = 1;    %ten interesting... smaller radius% one....large radius... %1;%100;%0.0001% .1;
end

if( onlyturning==1)
    par.d_cutoff   = 1; %default 50 %0.5;%0.7; %0.01;%10; %0.01
    par.e_cutoff   = par.d_cutoff/10.; %0.5;%0.7; %0.01;%10; %0.01
    par.v_cutoff   = 1;%0.5; % 2
    par.ve_cutoff  = 1;
    par.R_coeff    = 1;%1.;%1; %10 %0.2   % default 1
    par.S_coeff    = 1;   %default 100 
    par.angle_cut1 = 1.57;  % 80deg
    par.angle_cut2 = 1.047; % 60deg    %CHECK NOTATION.... 
    par.alpha      = 0.5;
    par.T_coeff    = 0;%1;%0.05;  %default 0.05 set to 1 for overcurly.
    par.C_coeff    = 1;
    par.fric_alpha = 1;
    par.fric_beta  = 0.01;%0.01;
    par.B_lambda  = 1;    %ten interesting... smaller radius% one....large radius... %1;%100;%0.0001% .1;
end


% Set constants
%d_cutoff = 10; %0.01
%v_cutoff = 2; % 2
%R_coeff = 0.2;
%S_coeff = 2;
%angle_cut1 = 1.57;  % 80deg
%angle_cut2 = 1.047; % 60deg
%alpha = 0.5;
%T_coeff = 8;
%C_coeff = 100;
%fric_alpha = 1;
%%fric_beta = 0.01;
