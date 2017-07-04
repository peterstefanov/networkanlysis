%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Compute normalized mutual information %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% normalisation - I(x,y)/ (1/2(H(x)+H(y))) %%%%%%%%%%%%%%%%%%%%%%%%%%%

function NMI = nmi(x, y)

n = numel(x);
x = reshape(x, 1, n);
y = reshape(y, 1, n);

min = min([x, y]);
x = x - min + 1;
y = y - min + 1;
max = max([x, y]);

indexVector = 1:n;
Mx = sparse(indexVector, x, 1, n, max, n);
My = sparse(indexVector, y, 1, n, max, n);
Pxy = nonzeros(Mx'*My/n);
Hxy = -(Pxy'*log2(Pxy));

%avoid the 0log0 
Px = nonzeros(mean(Mx, 1));
Py = nonzeros(mean(My, 1));

% entropy of Px  - dot product of Px'*(logPx)
Hx = -(Px'*log2(Px));
% entropy of Py  - dot product of Py'*(logPy)
Hy = -(Py'*log2(Py));

% mutual information
MI = Hx + Hy - Hxy;

% normalized mutual information
NMI = MI/((Hx + Hy)/2);