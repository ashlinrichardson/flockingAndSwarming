% Plot the frame. Input format is: 
function new_scale = show_frame(subplot_num,i,n,W,scale, scale2)

h1 = subplot(1,2,1);

% Setup x-y coordinates for plotting
XCoords = zeros(1,n);
YCoords = zeros(1,n);
ZCoords = zeros(1,n);
XVeloci = zeros(1,n);
YVeloci = zeros(1,n);
%ZVeloci = zeros(1,n);


% Show/make movie for reference frame
%if (subplot_num ~=0)
%  subplot(1,2,subplot_num)
%end

  
% Get x-y coordinates at time t
for j=1:n,  
    XCoords(j) = W(i,(1+6*(j-1))); 
    YCoords(j) = W(i,(3+6*(j-1)));
    ZCoords(j) = W(i,(5+6*(j-1)));
    XVeloci(j) = W(i,(2+6*(j-1))); 
    YVeloci(j) = W(i,(4+6*(j-1)));
%    ZVeloci(j) = W(i,(6+6*(j-1)));
end 
  
% Modify scaling if necessary (Gonna have to fix this for 3D)
new_scale = scale_axis_check(scale,XCoords,YCoords,ZCoords);
axis(h1, new_scale);
  
% Make scatter plot
scatter3(XCoords,YCoords,ZCoords,30,[0 0 1],'filled');

h2 = subplot(1,2,2);
new_scale2 = scale_axis_check2d(scale2,XCoords,YCoords, XVeloci, YVeloci);
axis(h2, new_scale2);
%mag = sqrt( XVeloci.^2 + YVeloci.^2);
quiver(h2, XCoords, YCoords, XVeloci, YVeloci,1,'b');

h1 = subplot(1,2,1);
