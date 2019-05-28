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
function [OMP_index]=LWOMP(y,A,p,sigma)
if exist('p','var')
    Bj_OMP=A'*y;
    [~,OMP_index]=max(abs(Bj_OMP)+(1/2)*(1+2*sigma^2)*log2(p./(1-p)));
    % Compressed Sensing With Prior Information Information-Theoretic Limits and Practical Decoders
    % Equation (16) in our manuscript
else
    Bj_OMP=A'*y;
    [~,OMP_index]=max(abs(Bj_OMP));
    % Equation (12) in our manuscript
end
end
