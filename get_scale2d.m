% From the initial conditions, get the maximum/minimum x and y
% coordinates so that the graph scale can be adjusted accordingly.

function scale = get_scale2d(ic)

% Put all possible initial x-positions in one vector
xpos = zeros(1,size(ic,2));
ypos = zeros(1,size(ic,2));
%zpos = zeros(1,size(ic,2));
xvel = zeros(1,size(ic,2));
yvel = zeros(1,size(ic,2));

for i=1:size(ic,2),
  xpos(i) = ic(1,i);
  ypos(i) = ic(3,i);
 % zpos(i) = ic(5,i);
  
  xvel(i) = ic(2,i);
  yvel(i) = ic(4,i);
end


dl = 100;
xmin = min([xpos+xvel xpos]) - dl;
xmax = max([xpos+xvel xpos]) + dl;
ymin = min([ypos+yvel ypos]) - dl;
ymax = max([ypos+yvel ypos]) + dl;

minv = [xmin ymin];% zmin];
maxv = [xmax ymax];% zmax];

axlen = abs(maxv-minv);
topv = minv + 2*max(axlen);

scale = [xmin-max(axlen) topv(1) ymin-max(axlen) topv(2)];% zmin topv(3)];
