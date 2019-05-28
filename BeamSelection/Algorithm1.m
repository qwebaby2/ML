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
% D: The number of phase quantization bit
% Phase_Shifts: The actual phase shift values
% N_bar: The number of codewords in the supercodebook F_bar or Q_bar
% N: The number of desired codewords
% M: The number of antennas
% WZ: The precoder or combiner matrix 
% Output Arguments:
% FQ_structured: The structured random training codebook F or Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [FQ_structured]=Algorithm1(D,Phase_Shifts,N_bar,N,M,WZ)
% FQ implies F or Q, similarly WZ implies W or Z
j=sqrt(-1);
%
FQ_bar=Random_Codebook(D,Phase_Shifts,N_bar,M);
%
Res=FQ_bar'*WZ;
norm_res=sqrt(diag(Res*Res'));
[~,ind]=sort(norm_res,'descend');
FQ_structured=FQ_bar(:,ind(1:N));