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
% C1: The number of clusters in channel 1
% C1: The number of clusters in channel 2
% tau_sigma1: The RMS delay spread of channel 1
% tau_sigma2: The RMS delay spread of channel 2
% gamma: frequency separation between Sub-6 GHz and mmWave
% sector_start: The starting point of the sector under consideration
% sector_end: The ending point of the sector under consideration
% Output Arguments:
% par1: The channel parameters of channel 1
% par2: The channel parameters of channel 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [par1,par2]=Algorithm2(C1,C2,tau_sigma1,tau_sigma2,gamma,sector_start,sector_end)
% _m becomes 1
% _s becomes 2
%% Generating both channel parameters independently
j=sqrt(-1);
% Channel 1
tau_c1 = tau_sigma1*log(rand(C1,1));
tau_c1 = sort(tau_c1-min(tau_c1));
%
theta_c1 = sector_start + rand(C1,1)*(sector_end - sector_start);
phi_c1 = sector_start + rand(C1,1)*(sector_end - sector_start);
% Channel 2
tau_c2 = tau_sigma2*log(rand(C2,1));
tau_c2 = sort(tau_c2-min(tau_c2));
%
theta_c2 = sector_start + rand(C2,1)*(sector_end - sector_start);
phi_c2 = sector_start + rand(C2,1)*(sector_end - sector_start);
%% Replacing clusters 
tau_max1=max(tau_c1);
tau_max2=max(tau_c2);
%
R1=[rand(C1,1)>gamma*(tau_c1/tau_max1)];
R2=[rand(C2,1)>gamma*(tau_c2/tau_max2)];
replacement_indices=intersect(find(R1==1),find(R2==1));
%
tau_c1(replacement_indices)=tau_c2(replacement_indices);
theta_c1(replacement_indices)=theta_c2(replacement_indices);
phi_c1(replacement_indices)=phi_c2(replacement_indices);
%% Adding perturbation
perturbation = rand(C1,1);
perturbation_time = gamma.*tau_c1.*perturbation;
perturbation_angle = gamma.*(tau_c1/tau_max1).*perturbation;
%
tau_ind_max = (tau_c1+perturbation_time > tau_max1);
tau_ind_min = (tau_c1-perturbation_time < 0);
tau_ind_minmax = (tau_c1+perturbation_time <= tau_max1 & tau_c1-perturbation_time >= 0);
tau_c1=tau_ind_max.*(tau_c1-perturbation_time)+tau_ind_min.*(tau_c1+perturbation_time)+tau_ind_minmax.*(tau_c1+sign(randn(C1,1)).*perturbation_time);
%
theta_ind_max = (theta_c1+perturbation_angle > sector_end);
theta_ind_min = (theta_c1-perturbation_angle < sector_start);
theta_ind_minmax = (theta_c1+perturbation_angle <= sector_end & theta_c1-perturbation_angle >= sector_start);
theta_c1=theta_ind_max.*(theta_c1-perturbation_angle)+theta_ind_min.*(theta_c1+perturbation_angle)+theta_ind_minmax.*(theta_c1+sign(randn(C1,1)).*perturbation_angle);
%
phi_ind_max = (phi_c1+perturbation_angle > sector_end);
phi_ind_min = (phi_c1-perturbation_angle < sector_start);
phi_ind_minmax = (phi_c1+perturbation_angle <= sector_end & phi_c1-perturbation_angle >= sector_start);
phi_c1=phi_ind_max.*(phi_c1-perturbation_angle)+phi_ind_min.*(phi_c1+perturbation_angle)+phi_ind_minmax.*(phi_c1+sign(randn(C1,1)).*perturbation_angle);
%
par1=[tau_c1 theta_c1 phi_c1];
par2=[tau_c2 theta_c2 phi_c2];
%
end

