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
function y = transition( lowval, hival, xbot, xtop, x)
    xorig=x;
    x = x - xbot;
    x = x ./ ( xtop - xbot);
    y = x;
    I = ( xorig >= xbot & xorig<xtop);
    y(I) = .5 + .5.*tanh( (1./(abs(x(I))-0)) +   (1./(abs(x(I))-1)) );
    y(I) = y(I).* (hival-lowval); 
    y = y + lowval;
    y( xorig < xbot )=hival;
    y( xorig >= xtop )=lowval;
 
   
