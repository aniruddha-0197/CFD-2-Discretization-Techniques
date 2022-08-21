function [V_norm, V_tan] = boundary(u_L,u_R,u_B,u_T,v_L,v_R,v_B,v_T,N,h_d)
%% Construct Boundary Vectors

% Normal Velocity
u_n 	= repmat([u_L,u_R],1,N);
v_n 	= repelem([v_B,v_T],N);
V_norm 	= [u_n,v_n];

% Tangential Velocity
u_t = repmat(u_B,1,N+1);
u_t = [u_t,u_T*h_d];

v_t   = repelem([v_L,v_R],N+1);
V_tan = [u_t,v_t];

end

