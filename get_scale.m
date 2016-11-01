% From the initial conditions, get the maximum/minimum x and y
% coordinates so that the graph scale can be adjusted accordingly.

function scale = get_scale(ic)

% Put all possible initial x-positions in one vector
xpos = zeros(1,size(ic,2));
ypos = zeros(1,size(ic,2));
zpos = zeros(1,size(ic,2));

for i=1:size(ic,2),
  xpos(i) = ic(1,i);
  ypos(i) = ic(3,i);
  zpos(i) = ic(5,i);
end

xmin = min(xpos) - 10;
xmax = max(xpos) + 10;
ymin = min(ypos) - 10;
ymax = max(ypos) + 10;
zmin = min(zpos) - 10;
zmax = max(zpos) + 10;

minv = [xmin ymin zmin];
maxv = [xmax ymax zmax];

axlen = abs(maxv-minv);
topv = minv + max(axlen);

scale = [xmin topv(1) ymin topv(2) zmin topv(3)];
