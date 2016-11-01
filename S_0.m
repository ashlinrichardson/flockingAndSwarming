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

%infinitely differentiable transition from one constant state to another.
%parameters:    hival - the high constant state.
%               loval - the low constant state.
%               xbot  - the x value at which loval occurs.
%               xtop  - the x value at which hival occurs.

function y = S_0( r, d_0, epsilon_0)%lowval, hival, xbot, xtop, x)
    if( min(r)<0)
        'S_0: Warning: negative radius encountered.'
    end
    y = transition( 0, 1, d_0, (d_0+epsilon_0), r);
    

