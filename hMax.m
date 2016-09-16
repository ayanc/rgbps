% function [q,lq,rho] = hMax(hin,ropts)
%
% Find peaks in albedo histogram.
%
% Copyright (C) 2016, Ayan Chakrabarti <ayanc@ttic.edu>
function [q,lq,rho] = hMax(hin,ropts)

LMAX=ropts.LMAX;
num = ropts.nms_num;

Q=sqrt(size(hin,1));
LQ=size(hin,2);

% Consider immediate neighbors
NGB=[-1:1]; sm=3;

% Reshape histogram into cube, and smooth slightly
hin = reshape(hin,[Q Q LQ]);
hin = (hin+convn(hin,ones([sm sm 1])/sm/sm,'same'))/2;


h = hin;
for i = NGB
  for j = NGB
    for k = NGB
      
      if abs(i+j+k) == 0
	continue;
      end;

      h(max(1,1-i):min(end,end-i), ...
	max(1,1-j):min(end,end-j), ...
	max(1,1-k):min(end,end-k))  = ...
	  h(max(1,1-i):min(end,end-i), ...
	    max(1,1-j):min(end,end-j), ...
	    max(1,1-k):min(end,end-k)) .* ...
	  ( hin(max(1,1+i):min(end,end+i), ...
		max(1,1+j):min(end,end+j), ...
		max(1,1+k):min(end,end+k)) < ...
	    hin(max(1,1-i):min(end,end-i), ...
		max(1,1-j):min(end,end-j), ...
		max(1,1-k):min(end,end-k)) );
      
    end;
  end;
end;


[~,idx] = sort(-h(:)); 
idx(h(idx) == 0) = [];
idx = idx(1:min(num,end));

[q,lq] = ind2sub([Q^2 LQ],idx);

if nargout == 3
  rho = qChrom(Q); rho = rho(q,:);
  rho = bsxfun(@times,rho, (lq-1)/(LQ-1)*LMAX );
  rho = reshape(rho,[1 size(rho)]);
end;