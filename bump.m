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

function ret = bump(x,c)
    x = abs(x);
    ret = x;
    I = (abs(x)<1);
    ret(I) = exp( -1./(c.*(1-(x(I).*x(I)))));
    
    I = (abs(x)>=1);
    ret( I ) = 0.;

