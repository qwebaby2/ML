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
% Phase_Shifts: The actual phase shift values
% Output Arguments:
% par1: The mmWave dictionary Z or W
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Dictionary]=Dictionary_Generation_m(sector_start,sector_end,...
    M,Phase_Shifts)
%% Generating equispaced sin(theta)
j=sqrt(-1);
ax=@(x) 1/sqrt(M)*exp(j*pi*[0:M-1]'*x); % Assuming lam/2 distance, x is sin(theta)
range=(sin(sector_end)-sin(sector_start))/M;
sampling_points=linspace(sin(sector_start)+...
    range/2,sin(sector_end)-range/2,M);
A=ax(sampling_points);
%% Quantizing based on phase shifters' resolution
AngA=angle(A);
AngAwr=wrapTo2Pi(AngA);
given_phase=[Phase_Shifts Phase_Shifts(end) + ...
    (Phase_Shifts(end)-Phase_Shifts(end-1))]; 
AngAwr_quantized = interp1(given_phase,given_phase,AngAwr,'nearest');
Dictionary=1/sqrt(M)*exp(j*AngAwr_quantized);