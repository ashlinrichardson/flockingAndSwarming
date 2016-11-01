function scaleOut = scale_axis_check2d(scale,XCoords,YCoords, XVel, YVel)

 % Fix scaling if necessary
 
 % Get xmin xmax ymin ymax
 xmin = scale(1);
 xmax = scale(2);
 ymin = scale(3);
 ymax = scale(4);
% zmin = scale(5);
% zmax = scale(6);

 %===============================================================

   xMx = max( [ max(XCoords+XVel) max(XCoords) ] );
   xMn = min( [ min(XCoords+XVel) min(XCoords) ] );
   yMx = max( [ max(YCoords+YVel) max(YCoords) ] ); 
   yMn = min( [ min(YCoords+YVel) min(YCoords) ] ); 
  
  if (xMx>xmax || xMn<xmin || yMx>ymax || yMn<ymin)
% || max(ZCoords)>zmax || min(ZCoords)<zmin)
   xmax = xMx + 100;
   xmin = xMn - 100;
   ymax = yMx + 100; 
   ymin = yMn - 100; 
 %  zmax = max(ZCoords) + 100; 
 %  zmin = min(ZCoords) - 100; 
  end
 
 minv = [xmin ymin];% zmin];
 maxv = [xmax ymax];% zmax];

 axlen = abs(maxv-minv);
 topv = minv + (1+1/10.)*max(axlen);

 scaleOut = [xmin-max(axlen)/10. topv(1) ymin-max(axlen)/10. topv(2)];% zmin topv(3)];

 return;
