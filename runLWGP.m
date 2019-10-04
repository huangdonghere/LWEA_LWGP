%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                   %
% This is a demo for the LWEA and LWGP algorithms. If you find this %
% code useful for your research, please cite the paper below.       %
%                                                                   %
% Dong Huang, Chang-Dong Wang, and Jian-Huang Lai.                  %
% "Locally weighted ensemble clustering."                           %
% IEEE Transactions on Cybernetics, 2018, 48(5), pp.1460-1473.      %
%                                                                   %
% The code has been tested in Matlab R2014a and Matlab R2015a on a  %
% workstation with Windows Server 2008 R2 64-bit.                   %
%                                                                   %
% https://www.researchgate.net/publication/316681928                %
%                                                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function labels = runLWGP(bcs,baseClsSegs, ECI, clsNum)

% Build the locally weighted bipartite graph
lwB = bsxfun(@times, baseClsSegs, ECI);

% Partition the graph
labels = zeros(size(bcs,1),numel(clsNum));
for i = 1:numel(clsNum) % clustering number. 
    disp(['Obtain ',num2str(clsNum(i)),' clusters by LWGP.']); 
    tic;
    labels(:,i) = bipartiteGraphPartitioning(lwB',clsNum(i));
    toc;
end 
disp('.');

function labels = bipartiteGraphPartitioning(B,Nseg)

% B - |X|-by-|Y|, cross-affinity-matrix

[Nx,Ny] = size(B);
if Ny < Nseg
    error('The cluster number is too large!');
end

dx = sum(B,2);
dx(dx==0) = 1e-10;
Dx = sparse(1:Nx,1:Nx,1./dx); clear dx
Wy = B'*Dx*B;

%%% compute Ncut eigenvectors
% normalized affinity matrix
d = sum(Wy,2);
D = sparse(1:Ny,1:Ny,1./sqrt(d)); clear d
nWy = D*Wy*D; clear Wy
nWy = (nWy+nWy')/2;

% computer eigenvectors
[evec,eval] = eig(full(nWy)); clear nWy % use eigs for large superpixel graphs  
[~,idx] = sort(diag(eval),'descend');
Ncut_evec = D*evec(:,idx(1:Nseg)); clear D

%%% compute the Ncut eigenvectors on the entire bipartite graph (transfer!)
evec = Dx * B * Ncut_evec; clear B Dx Ncut_evec

% normalize each row to unit norm
evec = bsxfun( @rdivide, evec, sqrt(sum(evec.*evec,2)) + 1e-10 );

% k-means
labels = kmeans(evec,Nseg,'MaxIter',100,'Replicates',3);