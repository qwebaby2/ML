%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Anum Ali                                           
% Last Modified: Nov, 2017
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the paper: Anum Ali, Nuria González-Prelcic and Robert W. Heath Jr., 
% "Millimeter Wave Beam-Selection Using Out-of-Band Spatial Information", 
% IEEE Transactions on Wireless Communications.
%
% Contact person email: anumali@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input Arguments:
% y: Measurements
% A: The sensing matrix
% p: the weights
% sigma: noise variance
% Output Arguments:
% OMP_index: The index selected by the OMP algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SOMP_index]=LWSOMP(y,A,p,sigma)
if exist('p','var')
    d=size(y,2);
    Bj_SOMP=A'*y;
    [~,SOMP_index]=max(sum(abs(Bj_SOMP),2)+(1/2)*(1+2*sigma^2)*log2(p./(1-p)));
    % Based on: Compressed Sensing With Prior Information Information-Theoretic Limits and Practical Decoders
    % Equation (18) in our manuscript
else
    Bj_SOMP=A'*y;
    [~,SOMP_index]=max(sum(abs(Bj_SOMP),2));
    % Simultaneous OMP algorithm
end
end
