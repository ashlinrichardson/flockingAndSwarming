function scaleOut = scale_axis_check(scale,XCoords,YCoords,ZCoords)

 % Fix scaling if necessary
 
 % Get xmin xmax ymin ymax
 xmin = scale(1);
 xmax = scale(2);
 ymin = scale(3);
 ymax = scale(4);
 zmin = scale(5);
 zmax = scale(6);
 
 %===============================================================
  
  if (max(XCoords)>xmax || min(XCoords)<xmin || max(YCoords)>ymax || min(YCoords)<ymin || max(ZCoords)>zmax || min(ZCoords)<zmin)
   xmax = max(XCoords) + 100;
   xmin = min(XCoords) - 100;
   ymax = max(YCoords) + 100; 
   ymin = min(YCoords) - 100; 
   zmax = max(ZCoords) + 100; 
   zmin = min(ZCoords) - 100; 
  end
 
 minv = [xmin ymin zmin];
 maxv = [xmax ymax zmax];

 axlen = abs(maxv-minv);
 topv = minv + max(axlen);

 scaleOut = [xmin topv(1) ymin topv(2) zmin topv(3)];

 return;
  %scaleOut = [xmin xmax ymin ymax zmin zmax];
