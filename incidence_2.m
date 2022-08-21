function [E21, tE10, E_tan] = incidence_2(N)
%% Construct E21,tE10 and E_tan

% E21 and tE10
B1_1 = spdiags(ones((N+1)^2,1)*[-1,1],[0,-(N+1)],(N+1)^2,N*(N+1));
D    = spdiags(ones(N+1,1)*[1,-1],[0,-1],(N+1),N);
B1_2 = kron(speye(N+1),D);

E21   = [B1_1,B1_2];		% Dual Curl
tE10  = E21.';              % Primal Curl

% Boundary Matrix (Tangential)
D2 	  = sparse([1,N+1],[1,2],[1,-1],N+1,2);
B2_1  = kron(D2,speye(N+1));
B2_2  = kron(speye(N+1),-D2);
E_tan = [B2_1,B2_2];

end

