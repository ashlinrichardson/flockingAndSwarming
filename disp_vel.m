function final_vels = disp_vel(W,n) 
% Show final velocities for all birds (Column: Vx,Vy,norm(V))
tfinal = size(W,1);
final_vels = zeros(3,n);

for i=1:n,
  
  % X-comp of velocity
  final_vels(1,i) = W(tfinal,((i-1)*4)+2);
  
  % Y- comp of velocity
  final_vels(2,i) = W(tfinal,((i-1)*4)+4);
  
  % Speed of bird
  final_vels(3,i) = ( (final_vels(2,i)*final_vels(2,i)) + (final_vels(1, ...
						  i)*final_vels(1,i)) ).^(1/2);   
end
