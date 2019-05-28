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
% sector_start: The starting point of the sector under consideration
% sector_end: The ending point of the sector under consideration
% M: The number of antennas
% no_codewords: The desired number of codewords (typically same as the number of codewords
% at mmWave)
% Output Arguments:
% par1: The Sub-6 GHz dictionary underline(Z) or underline(W)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dictionary]=Dictionary_Generation_s(sector_start,sector_end,M,no_codewords)
%% Generating equispaced sin(theta)
j=sqrt(-1);
ax=@(x) 1/sqrt(M)*exp(j*pi*[0:M-1]'*x); % Assuming lam/2 distance, x is sin(theta)
range=(sin(sector_end)-sin(sector_start))/no_codewords;
sampling_points=linspace(sin(sector_start)+range/2,sin(sector_end)-range/2,no_codewords);
Dictionary=ax(sampling_points);