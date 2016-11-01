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
function retval = runme_single_02_restore()
    load runme_single_last_data.mat  
    initial_conditions = A;
    T = 3300;
    n_interac_3D_mov( initial_conditions, par, T, T/5);  
