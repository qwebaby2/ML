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
% Mrx: Number of receive antennas
% Mtx: Number of transmit antennas
% L: The number of channel taps
% eval_points: The points where the continous channel is to be evaluated
% C: The number of clusters
% Rc: The number of rays in a cluster
% alpha: The complex path coefficients of the rays
% sintheta: The AoAs (sin operated)
% sinphi: The AoDs (sin operated)
% BW: The bandwidth
% tau: The delays
% beta: the roll-off factor of the raised cosine filter
% K: The number of subcarriers
% Output Arguments:
% H: The time domain MIMO channel taps
% H_freq: The frequency domain MIMO channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H,H_freq]=MIMO_channel(Mrx,Mtx,L,eval_points,...
    C,Rc,alpha,sintheta,sinphi,BW,tau,beta,K)
%% Defining the raised cosine filter
prc=@(x,beta) nansum([double((abs(x)==1/2/beta).* pi/4*sinc(1/2/beta))'...
    double(~(abs(x)==1/2/beta).*cos(pi*beta*x)./...
    (1-(2*beta*x).^2).*sinc(x))'],2); % https://en.wikipedia.org/wiki/Raised-cosine_filter
%
j=sqrt(-1);
%
aTx=@(x) 1/sqrt(Mtx)*exp(j*pi*[0:Mtx-1]'*x); % Assuming lam/2 distance, x is sin(theta)
aRx=@(x) 1/sqrt(Mrx)*exp(j*pi*[0:Mrx-1]'*x);
%
H=zeros(Mrx,Mtx,L);
h=zeros(1,L);
%
for ii=1:L
    for jj=1:C*Rc
        H(:,:,ii)=H(:,:,ii)+sqrt(Mrx*Mtx)*alpha(jj)*aRx(sintheta(jj))*aTx(sinphi(jj))'...
            *prc(BW*(tau(jj)-eval_points(ii)/(1+beta)/BW),beta);
        h(ii)=h(ii)+alpha(jj)*prc(BW*(tau(jj)-eval_points(ii)/(1+beta)/BW),beta);
        % The argument of pulse shape is BW(tau_k - n x Ts), where n is the tap
        % delay index, Ts is the sampling time, now for raised cosine with
        % roll off beta, in BW B, Ts=1/((1+beta)xB), so the argument becomes
        % BW(tau_k - n / (1+beta)/B)
        % "A MATLAB - based Object-Oriented Approach to Multipath Fading
        % Channel Simulation"
    end
end
%
H=H/norm(h)*norm(alpha);
%
H_freq=zeros(Mrx,Mtx,K);
for ii=1:K
    for jj=1:L
        H_freq(:,:,ii)=H_freq(:,:,ii)+H(:,:,jj)*exp(-j*2*pi/K*(ii-1)*(jj-1));
        % ii-1, meaning the subcarriers are from 0 to N-1
        % Fourier transform formula mentioned in equation (4) in
        % Frequency Selective Hybrid Precoding for Limited Feedback Millimeter Wave Systems
    end
end
end