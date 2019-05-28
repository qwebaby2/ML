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
% N: The number of desired codewords
% M: The number of antennas
% Output Arguments:
% FQ: The random training codebook F or Q
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function FQ=Random_Codebook(D,Phase_Shifts,N,M)
j=sqrt(-1);
Integers_for_Random_Phases = randi(2^D,[M,N]);
FQ=sqrt(1/M)*exp(j*Phase_Shifts(Integers_for_Random_Phases));
