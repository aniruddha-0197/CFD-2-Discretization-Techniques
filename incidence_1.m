function [tE21, E10, E_norm] = incidence_1(N)
%% Construct tE21,E10 and E_norm

% First Block of tE21
D1 	 = spdiags(ones(N,1)*[-1,1],[0,1],N,N+1);
I 	 = speye(N);
B1_1 = kron(I,D1);
B1_2 = spdiags(ones(N^2,1)*[-1,1],[0,N],N^2,N*(N+1));
B1 	 = [B1_1 B1_2];

% Second BLock of tE21
D2 	 = sparse([1,2],[1,N+1],[1,-1],2,N+1);
B2_1 = kron(I,D2);
B2_2 = kron(D2,I);
B2   = blkdiag(B2_1,B2_2);

tE21 = [B1;B2];	     % Primal Divergence
E10  = -tE21.';      % Dual Gradient 

% Boundary Matrix (Normal)
diag   = [repmat([-1,1],1,N),repelem([-1,1],N)];
E_norm = [sparse(N^2,4*N);spdiags(diag.',0,length(diag),length(diag))];

end

