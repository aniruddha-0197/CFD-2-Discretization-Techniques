function [Vel,up,vp] = vel_int(u,hd,N)
%% Interpolate velocity from dual edges onto inner dual points and compute magnitude

% Initialize
up = zeros(N^2,1); vp = up;

% Interpolate
for i = 1:N
    for j = 1:N
        
        k = j + (i-1)*N;

        up(k) = (u(j + (i-1)*(N+1))/hd(j) + u(j+1 + (i-1)*(N+1))/hd(j+1))/2;
        vp(k) = (u(N*(N+1) +j + (i-1)*N)/hd(i) + u(N*(N+1) + j + i*N)/hd(i+1))/2;
     
    end
end

% Velocity Magnitude
Vel = sqrt(up.^2 + vp.^2);

end

