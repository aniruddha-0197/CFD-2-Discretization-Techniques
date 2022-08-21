function [Ht1_1, H1_t1, Ht0_2] = hodge(N,h_p,h_d)
%% Construct Hodge Matrices

% H1_t1 and Ht1_1
d = zeros(1,2*N*(N+1));
k = 1;

for i = 1:N

    for j = 1:N+1
        d(k) = h_p(i)/h_d(j);
        k    = k+1;
    end

end

for i = 1:N+1
    for j = 1:N
	    d(k) = h_p(j)/h_d(i);
	    k    = k+1;
    end
end

di  = 1./d;
len = length(d); 

Ht1_1 = spdiags(d.',0,len,len);	    % Dual to Primal
H1_t1 = spdiags(di.',0,len,len);	% Primal to Dual

% Ht0_2
S_d   = (h_d.'*h_d).';
S_d   = S_d(:);
Ht0_2 = spdiags(1./S_d,0,length(S_d),length(S_d));

end

