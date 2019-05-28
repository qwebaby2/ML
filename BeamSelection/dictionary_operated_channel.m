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
% H_freq: The frequency domain channel
% v: noise
% W: precoding dictionary
% Z: combining dictionary
% K: number of subcarriers
% subcarrier_index: The index of the subcarrier to execute the function (if 0, execute 
% on all subcarriers)
% noise_type_ind: noise type indicator (i) 0 means Z'HW+Z'v, (ii) means Z'(H+v)W
% Output Arguments:
% ZHW: The dictionary operated channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ZHW]=dictionary_operated_channel(H_freq,v,W,Z,K,subcarrier_index,noise_type_ind)

if subcarrier_index~=0
    ii=subcarrier_index;
    if(noise_type_ind==0)
        ZHW=Z'*H_freq(:,:,ii)*W+Z'*v(:,:,ii);
    elseif (noise_type_ind==1)
        ZHW=Z'*(H_freq(:,:,ii)+v(:,:,ii))*W;
    end
else
    ZHW=zeros(size(Z',1),size(W,2),K);
    for ii=1:K
        if(noise_type_ind==0)
            ZHW(:,:,ii)=Z'*H_freq(:,:,ii)*W+Z'*v(:,:,ii);
        elseif (noise_type_ind==1)
            ZHW(:,:,ii)=Z'*(H_freq(:,:,ii)+v(:,:,ii))*W;
        end
    end
end
