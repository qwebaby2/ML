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
% par: The parameters of the channel (generated via Algorithm 2)
% tau_sigma: The RMS delay spread of the channel
% mu: The exponential PDF parameter
% C: The number of clusters
% Rc: The number of rays in a cluster
% tau_sigma_clust: The RMS delay spread of a cluster
% sigma_theta: The RMS AoA spread
% sigma_phi: The RMS AoD spread
% Output Arguments:
% homni_par: The parameters of the omni directional channel (eq. 19)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [homni_par]=omni_directional_channel_par...
    (par,tau_sigma,mu,C,Rc,tau_sigma_clust,sigma_theta,sigma_phi)
j=sqrt(-1);
%
tau_c=par(:,1);
theta_c=par(:,2);
phi_c=par(:,3);
%
PDP=exppdf(tau_c/tau_sigma,mu); % Exponential PDP
PDP=PDP/sum(PDP); % Normalizing to unit power PDP
%
tau_rc = tau_sigma_clust*randn(Rc,C);
tau_c=ones(Rc,1)*tau_c';
taus=tau_rc+tau_c;
taus=taus(:);
taus=taus-min(taus);
%
vartheta_rc=angle(exp(j*sigma_theta*randn(Rc,C)));
theta_c=ones(Rc,1)*theta_c';
thetas=vartheta_rc+theta_c;
thetas=wrapToPi(thetas(:));
%
varphi_rc=angle(exp(j*sigma_phi*randn(Rc,C)));
phi_c=ones(Rc,1)*phi_c';
phis=varphi_rc+phi_c;
phis=wrapToPi(phis(:));
%
alphas=zeros(Rc,C);
%
for ii=1:C
    alphas(:,ii)=sqrt(PDP(ii)/2/Rc)*(randn(Rc,1)+j*randn(Rc,1));
end
%
alphas=alphas(:);
%
homni_par=[alphas taus thetas phis];