% Important: initial conditions ic must be passed in the form of a
% matrix with each column representing a different bird and each
% bird having it's column vector of the form (x,xdot,y,ydot) where
% x and y are the coordinates and xdot and ydot are the time
% derivatives. Therefore this matrix will be of dimensions 4xn
% where n is the number of birds.

function nothing = n_interac_3D_mov_1(ic,tf)

% Set constants
par = getparameters();

% Number of birds
global n;
n = size(ic,2);

% Make number of time points where ODE will be evaluated
time_stamps = 500;
times = zeros(1,time_stamps);
for i=1:time_stamps,
  times(i)=(i-1)*(tf/(time_stamps-1));
end
%'solving ode...'
% Solve ODE in inertial frame
[T, W] = ode45(@(times,ic) n_interac_3D_1(times,ic,par), times, ic);
%'done solving ode...'

% ------------------ PLAY MOVIE ------------------ %
% Prepare movie
numframes = size(T,1);
movie_inertial=moviein(numframes);

% Setup the movies and adjust the plotting scale
figure(1)
set(gca,'NextPlot','replacechildren')
scaleI = get_scale(ic);
axis(scaleI)
title('Flocking Simulation')

% Plot and run movie, loop over all values of t
for i=1:numframes,
  %i
  
  % Create and display frame of movie for CM frame
  scaleI = show_frame(0,i,n,W,scaleI);
   
  % Get frame for the movie 
  movie_inertial(:,i)=getframe;

end

cla;
movie(movie_inertial,1,10)
save birdsmov.mat movie_inertial
% ------------------ PLAY MOVIE ------------------ %

