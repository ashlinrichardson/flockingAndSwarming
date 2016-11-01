% Important: initial conditions ic must be passed in the form of a
% matrix with each column representing a different bird and each
% bird having it's column vector of the form (x,xdot,y,ydot) where
% x and y are the coordinates and xdot and ydot are the time
% derivatives. Therefore this matrix will be of dimensions 4xn
% where n is the number of birds.

function nothing = n_interac_3D_mov(ic,par,tf,ts)

% Set constants
%par = getparameters();

% Number of birds
global n;
n = size(ic,2);

% Make number of time points where ODE will be evaluated
time_stamps = ts;%200;
times = zeros(1,time_stamps);
for i=1:time_stamps,
  times(i)=(i-1)*(tf/(time_stamps-1));
end

% Solve ODE in inertial frame
[T, W] = ode45(@(times,ic) n_interac_3D(times,ic,par), times, ic);


% ------------------ PLAY MOVIE ------------------ %
% Prepare movie


numframes = size(T,1);
movie_inertial=moviein(numframes);

% Setup the movies and adjust the plotting scale
h_fig = figure(1)
scrsz = get(0,'ScreenSize');
figure('Position',[1 scrsz(4)/2 scrsz(3) scrsz(3)])
h1 = subplot(1,2,1)
%[left, bottom, width, height]

set(h1,'NextPlot','replacechildren')
scaleI = get_scale(ic);
axis(h1, scaleI)
title('Flocking Simulation')

h2 = subplot(1,2,2)
scaleI2 = get_scale2d(ic);
axis(h2, scaleI2)
set(h2,'NextPlot','replacechildren')

% Plot and run movie, loop over all values of t
for i=1:numframes,
  
  % Create and display frame of movie for CM frame
  scaleI = show_frame(0,i,n,W,scaleI, scaleI2);
   
  % Get frame for the movie 
  movie_inertial(:,i)=getframe(gcf);
  

end

%un comment this for AVI file creation...........

'writing avi file...'
aviobj = avifile('mymovie.avi','FPS',10);
for i=1:numframes,
    i
    aviobj = addframe(aviobj,movie_inertial(:,i));
end

aviobj = close(aviobj);
'done writing avi file...'

cla;
%movie2avi( movie_inertial, 'movie_inertial.avi');
movie(movie_inertial);%,1,10)
save birdsmov.mat movie_inertial

% ------------------ PLAY MOVIE ------------------ %
