function [apen, result] = ApEn( dim, r, signal)

%   dim : embedded dimension
%   r : tolerance (typically 0.2 * std)
%   signal : time-series signal
%   tau : delay time for downsampling

%   Changes in version 1
%       Ver 0 had a minor error in the final step of calculating ApEn
%       because it took logarithm after summation of phi's.
%       In Ver 1, I restored the definition according to original paper's
%       definition, to be consistent with most of the work in the
%       literature. Note that this definition won't work for Sample
%       Entropy which doesn't count self-matching case, because the count 
%       can be zero and logarithm can fail.
%
%       A new parameter tau is added in the input argument list, so the users
%       can apply ApEn on downsampled signal by skipping by tau. 
%---------------------------------------------------------------------
% coded by Kijoon Lee,  kjlee@ntu.edu.sg
% Ver 0 : Aug 4th, 2011
% Ver 1 : Mar 21st, 2012
%---------------------------------------------------------------------
N = size(signal,2); %number of columns
result = zeros(1,2);

for j = 1:2
    m = dim+j-1;
    phi = zeros(1,N-m+1);
    signal_emb = zeros(m,N-m+1);
    
    % setting up signal matrix
    for i = 1:m
        signal_emb(i,:) = signal(i:N-m+i); 
    end
    
    % counting similar patterns using distance calculation
    for i = 1:N-m+1
        tempMat = abs(signal_emb - repmat(signal_emb(:,i),1,N-m+1));
        boolMat = any( (tempMat > r),1);
        phi(i) = sum(~boolMat)/(N-m+1);
    end
    
    % summing over the counts
    result(j) = sum(log(phi))/(N-m+1);
end

apen = result(1)-result(2);

end
    
    