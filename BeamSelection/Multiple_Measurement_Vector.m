%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Anum Ali                                           
% Last Modified: Nov, 2017
%
% If you use this code or any (modified) part of it in any publication, please cite 
% the paper: Anum Ali, Nuria Gonz�lez-Prelcic and Robert W. Heath Jr., 
% "Millimeter Wave Beam-Selection Using Out-of-Band Spatial Information", 
% IEEE Transactions on Wireless Communications.
%
% Contact person email: anumali@utexas.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This scripts generates the results presented in Fig. 14 in the manuscript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the code, variables ending with _m are for MmWave and _s are for Sub-6 GHz
clc;
clear all;
close all;
%
j=sqrt(-1);
rng(1); % fixing the random seed
%% Simulation Parameters
E=1e1;2e3; % Number of trials/runs/iterations
%
Nrx_vec=[8:4:24]; % The number of measurements at the RX
Ntx_vec=[8:4:24]/4; % The number of measurements at the TX
%
Tx_Rx_dist=60; % The distance between TX and RX in meters
%% System Input Parameters
K_m=128; % Number of sub-carriers in OFDM system at mmWave
Lc_m=K_m/4; % Cyclic prefix length for mmWave
%
K_s=32; % Number of sub-carriers in OFDM system at sub-6 GHz
Lc_s=K_s/4; % Cyclic prefix length for sub-6 GHz
%
f_s=3.5e9; % Sub-6 GHz frequency
BW_s=150e6; % Sub-6 GHz bandwidth
%
f_m=28e9; % MmWave frequency
BW_m=850e6; % MmWave bandwidth
%
Mrx_s=8; % Number of Rx antennas at Sub-6 GHz
Mtx_s=2; % Number of Tx antennas at Sub-6 GHz
%
Mrx_m=64; % Number of Rx antennas at MmWave
Mtx_m=16; % Number of Tx antennas at MmWave
%
Drx=2;  % Number of phase shifters quantization bits at the Rx
Dtx=2;  % Number of phase shifters quantization bits at the Tx
%
sector_start=-pi/3; % Sectoar start point
sector_end=pi/3; % Sector end point
%
beta_s=1; % roll-off factor of the pulse shaping filter at Sub-6 GHz
beta_m=1; % roll-off factor of the pulse shaping filter at MmWave
%
Pt_dBm_m=43-10*log10(Mtx_m); % Total Sub-6 GHz transmit power
Pt_dBm_s=30+10*log10(BW_s/25e6)-10*log10(Mtx_s); % Total MmWave transmit power
%% Channel Input Parameters
C_s=10; % Number of clusters at Sub-6 GHz
C_m=5; % Number of clusters at MmWave
%
Rc_s=20; % Number of rays in a cluster at Sub-6 GHz
Rc_m=20; % Number of rays in a cluster at MmWave
%
sigma_theta_m=2*pi/180; % RMS arrival angle spread at MmWave
sigma_phi_m=2*pi/180; % RMS departure angle spread at MmWave 
%
sigma_theta_s=4*pi/180; % RMS arrival angle spread at Sub-6 GHz
sigma_phi_s=4*pi/180; % RMS departure angle spread at Sub-6 GHz
%
pl_coeff_s=3; % Path-loss coefficient for Sub-6 GHz
pl_coeff_m=3; % Path-loss coefficient for MmWave
%
gamma=abs(f_s-f_m)/max(f_s,f_m); % Frequency separation between Sub-6 GHz and mmWave
% The frequency separation is used in perturbation generation
%
mu_s=0.1; % exponential distribution parameter at Sub-6 GHz
mu_m=0.2; % exponential distribution parameter at MmWaves
% Used in power-delay-profile (PDP) generation
%% Miscellaneous Parameters Calculated based on input
L_m=Lc_m+1; %  Number of taps in the wideband mmWave channel
eval_points_m=0:L_m-1;
%
L_s=Lc_s+1; %  Number of taps in the wideband Sub-6 GHz channel
eval_points_s=0:L_s-1;
%
Phase_Shifts_tx=2*pi*[0:(2^Dtx)-1]/(2^Dtx); % The possible phase-shifts for mmWave TX
Phase_Shifts_rx=2*pi*[0:(2^Drx)-1]/(2^Drx); % The possible phase-shifts for mmWave RX
%
tau_max_m=eval_points_m(end)/(1+beta_m)/BW_m; % Maximum delay-spread at MmWave
tau_max_s=eval_points_s(end)/(1+beta_s)/BW_s; % Maximum delay-spread at Sub-6 GHz
%
tau_sigma_m=1/(log(1-0.999)/(-tau_max_m)); % RMS delay-spread at MmWave
tau_sigma_s=1/(log(1-0.999)/(-tau_max_s)); % RMS delay-spread at Sub-6 GHz
%
tau_sigma_clust_m=tau_sigma_m/10; % RMS delay-spread of a cluster at MmWave
tau_sigma_clust_s=tau_sigma_s/10; % RMS delay-spread of a cluster at Sub-6 GHz
%% SNR Calculations
Pt_avg_m=10^((Pt_dBm_m-30)/10); % MmWave transmit power in watts
Pt_avg_s=10^((Pt_dBm_s-30)/10); % Sub-6 GHz transmit power in watts
%
lam_m=3e8/f_m; % MmWave wavelength
lam_s=3e8/f_s; % Sub-6 GHz wavelength
%
d0=5; % Reference distance
%
pl_m=(4*pi*d0/lam_m)^2*(Tx_Rx_dist/d0)^pl_coeff_m; % MmWave path-loss
pl_s=(4*pi*d0/lam_s)^2*(Tx_Rx_dist/d0)^pl_coeff_s; % Sub-6 GHz path-loss
%
Pr_avg_dBm_m=Pt_dBm_m-10*log10(pl_m); % Average received power MmWave
Pr_avg_dBm_s=Pt_dBm_m-10*log10(pl_s); % Average received power Sub-6 GHz
% linear version of alpha+10*beta*log10(d/d0), where alpha=10log10(4*pi*d0/lambda)^2
%
N0_dBm_m=-173.8+10*log10(BW_m); % MmWave Noise power in dBm
% 10*log10(kT)+30=-173.8, where k=1.381*10^(-23) is boltzman constant, and
% T=300 is temperature in kelvins, +30 is used to convert into dBm
% Total N0_dBm=10*log10(kTB)+30=10*log10(kT)+30*10*log10(B)
N0_dBm_s=-173.8+10*log10(BW_s); % Sub-6 GHz Noise power in dBm
%
SNR_dB_m=Pr_avg_dBm_m-N0_dBm_m; % MmWave SNR in dB
SNR_dB_s=Pr_avg_dBm_s-N0_dBm_s; % Sub-6 GHz SNR in dB
%
SNR_m=10^(SNR_dB_m/10); % MmWave SNR in linear scale
SNR_s=10^(SNR_dB_s/10); % Sub-6 GHz SNR in linear scale
%% Array Related Parameters (Used for Calculating Rate)
[W_m]=Dictionary_Generation_m(sector_start,sector_end,Mtx_m,Phase_Shifts_tx); % Transmit dictionary
[Z_m]=Dictionary_Generation_m(sector_start,sector_end,Mrx_m,Phase_Shifts_rx); % Receive dictionary
% This TX/RX dictionary corresponds to exhaustive search (eq. 7 in the manuscript)
[W_s]=Dictionary_Generation_s(sector_start,sector_end,Mtx_s,Mtx_m);
[Z_s]=Dictionary_Generation_s(sector_start,sector_end,Mrx_s,Mrx_m);
% The Sub-6 GHz dictionary is used in information retrieval (eq. 15 in the manuscript)
WconjkronZ_m=sqrt(Mtx_m*Mrx_m)*kron(transpose(W_m),Z_m')';
%% Variables to store results
Rate_exh_sear=zeros(length(Nrx_vec),1); % Rate for exhaustive-search
Rate_SOMP=zeros(length(Nrx_vec),1); % Rate for Simultaneous OMP
Rate_StrucLWSOMP=zeros(length(Nrx_vec),1); % Rate for Simultaneous structured LW-OMP
%% Loop over Measurements
for ee=1:E
    %% Channel's Parameter Generation
    % Generating the mean time of arrival and AoA/AoD of the clusters at
    % sub-6 GHz based on Algorithm 2 in the manuscript
    [par_m,par_s]=Algorithm2(C_m,C_s,tau_sigma_m,tau_sigma_s,gamma,sector_start,sector_end);
    %% Omni-directional channel impulse response par
    % Getting the omni-directional channel impulse response parametres
    % (used in eq. 19 in the manuscript)
    % For MmWave
    [homni_par_m]=omni_directional_channel_par(par_m,tau_sigma_m,mu_m,C_m,Rc_m,tau_sigma_clust_m,sigma_theta_m,sigma_phi_m);
    alpha_m=homni_par_m(:,1);tau_m=homni_par_m(:,2);sintheta_m=sin(homni_par_m(:,3));sinphi_m=sin(homni_par_m(:,4));
    % For Sub-6 GHz
    [homni_par_s]=omni_directional_channel_par(par_s,tau_sigma_s,mu_s,C_s,Rc_s,tau_sigma_clust_s,sigma_theta_s,sigma_phi_s);
    alpha_s=homni_par_s(:,1);tau_s=homni_par_s(:,2);sintheta_s=sin(homni_par_s(:,3));sinphi_s=sin(homni_par_s(:,4));
    %
    %% MIMO Channel Generation
    [H_m,H_freq_m]=MIMO_channel(Mrx_m,Mtx_m,L_m,eval_points_m,C_m,Rc_m,alpha_m,sintheta_m,sinphi_m,BW_m,tau_m,beta_m,K_m);
    [H_s,H_freq_s]=MIMO_channel(Mrx_s,Mtx_s,L_s,eval_points_s,C_s,Rc_s,alpha_s,sintheta_s,sinphi_s,BW_s,tau_s,beta_s,K_s);
    %% Noise Generation
    v_s=sqrt(K_s)/sqrt(2*SNR_s)*(randn(Mrx_s,Mtx_s,K_s)+j*randn(Mrx_s,Mtx_s,K_s));
    v_m=sqrt(K_m)/sqrt(2*SNR_m)*(randn(Mrx_m,Mtx_m,K_m)+j*randn(Mrx_m,Mtx_m,K_m));
    % Noise Power is further multiplied by sqrt(K), because the overall
    % signal power is divided into K subcarriers, hence the SNR on
    % each tone is K times weaker compared to the total symbol SNR
    %% Random Training Codebooks
    F=Random_Codebook(Dtx,Phase_Shifts_tx,Mtx_m,Mtx_m);
    Q=Random_Codebook(Drx,Phase_Shifts_rx,Mrx_m,Mrx_m);
    %% Fetching out-of-band information
    %% The dominant directions/index from Sub-6 GHz
    [ZHW_s]=dictionary_operated_channel(H_freq_s,v_s,W_s,Z_s,K_s,0,1);
    ZHW_s_abs=mean(abs(ZHW_s).^2,3);
    %
    [~,Res_s]=sort(ZHW_s_abs(:),'descend');
    Res_s=Res_s(1:Mrx_m/Mrx_s); % O=Mrx_m/Mrx_s dominant precoders are selected
    [I_s,J_s] = ind2sub(size(ZHW_s_abs),Res_s); % Indices of dominant precoders and combiners
    %% Structured Random Training Codebooks based on out-of-band information (Algorithm 1 in the manuscript)
    N_bar_tx=1e2*Mtx_m;
    N_bar_rx=1e2*Mrx_m;
    %
    F_structured=Algorithm1(Dtx,Phase_Shifts_tx,N_bar_tx,Mtx_m,Mtx_m,W_m(:,J_s));
    Q_structured=Algorithm1(Drx,Phase_Shifts_rx,N_bar_rx,Mrx_m,Mrx_m,Z_m(:,I_s));
    %% Weight calculation for structured LW-OMP (based on eq. 17 in the manuscript)
    J_p=0.1;
    p=ZHW_s_abs-min(min(ZHW_s_abs));
    p=p/max(max(ZHW_s_abs));
    p=J_p*p;
    p=p(:);
    %% Randomly chosen subcarrier for sensing
    subcarrier_index=0;
    %% Main Loop
    for mm=1:length(Nrx_vec)
        [ee mm]
        %
        Nrx=Nrx_vec(mm);
        Ntx=Ntx_vec(mm);
        %% Exhaustive search at infinite SNR (no noise) to get the best beam-direction
        [ZHW_exh_sear_no_noise_m]=dictionary_operated_channel(H_freq_m,zeros(Mrx_m,Mtx_m,K_m),W_m,Z_m,K_m,subcarrier_index,0);
        ZHW_exh_sear_no_noise_abs_m=mean(abs(ZHW_exh_sear_no_noise_m).^2,3);
        [~,Beam_indices_sorted]=sort(abs(ZHW_exh_sear_no_noise_abs_m(:)),'descend');
        %
        Best_beam_index=Beam_indices_sorted(1);
        Best_five_beams_index=Beam_indices_sorted(1:5);
        %% Exhaustive Search
        [ZHW_exh_sear_m]=dictionary_operated_channel(H_freq_m,v_m,W_m,Z_m,K_m,subcarrier_index,0);
        ZHW_exh_sear_abs_m=mean(abs(ZHW_exh_sear_m).^2,3);
        Result_exh_sear=find(ZHW_exh_sear_abs_m(:)==max(max(ZHW_exh_sear_abs_m)));
        [I_exh_sear,J_exh_sear] = ind2sub(size(ZHW_exh_sear_abs_m),Result_exh_sear);
        %
        for ii=1:K_m
            Rate_exh_sear(mm)=Rate_exh_sear(mm)+1/E*1/K_m*log2(1+abs(Z_m(:,I_exh_sear)'*H_freq_m(:,:,ii)*W_m(:,J_exh_sear))^2*SNR_m/K_m);
        end
        %% SOMP
        [ZHW_SOMP_m]=dictionary_operated_channel(H_freq_m,v_m(:,1:Ntx,:),F(:,1:Ntx),Q(:,1:Nrx),K_m,subcarrier_index,0);
        ZHW_SOMP_m=reshape(ZHW_SOMP_m(:),[Ntx*Nrx K_m]);
        FtranskronQherm=kron(transpose(F(:,1:Ntx)),Q(:,1:Nrx)');
        Psi=FtranskronQherm*WconjkronZ_m;
        onebyPsi_column_norms=sqrt(diag(1./diag(Psi'*Psi)));
        Psi=Psi*onebyPsi_column_norms;
        %
        Result_SOMP=LWSOMP(ZHW_SOMP_m,Psi);
        [I_SOMP,J_SOMP] = ind2sub(size(ZHW_exh_sear_no_noise_abs_m),Result_SOMP);
        %
        for ii=1:K_m
            Rate_SOMP(mm)=Rate_SOMP(mm)+1/E*1/K_m*log2(1+abs(Z_m(:,I_SOMP)'*H_freq_m(:,:,ii)*W_m(:,J_SOMP))^2*SNR_m/K_m);
        end
        %% SLW-OMP
        [ZHW_StrucLWSOMP_m]=dictionary_operated_channel(H_freq_m,v_m(:,1:Ntx,:),F_structured(:,1:Ntx),Q_structured(:,1:Nrx),K_m,subcarrier_index,0);
        ZHW_StrucLWSOMP_m=reshape(ZHW_StrucLWSOMP_m(:),[Ntx*Nrx K_m]);
        FstructranskronQstrucherm=kron(transpose(F_structured(:,1:Ntx)),Q_structured(:,1:Nrx)');
        Psistruc=FstructranskronQstrucherm*WconjkronZ_m;
        onebyPsistruc_col_norms=sqrt(diag(1./diag(Psistruc'*Psistruc)));
        Psistruc=Psistruc*onebyPsistruc_col_norms;
        %
        Result_StrucLWSOMP=LWSOMP(ZHW_StrucLWSOMP_m,Psistruc,p,sqrt(K_m/SNR_m));
        %
        [I_StrucLWSOMP,J_StrucLWSOMP] = ind2sub(size(ZHW_exh_sear_no_noise_abs_m),Result_StrucLWSOMP);
        %
        for ii=1:K_m
            Rate_StrucLWSOMP(mm)=Rate_StrucLWSOMP(mm)+1/E*1/K_m*log2(1+abs(Z_m(:,I_StrucLWSOMP)'*H_freq_m(:,:,ii)*W_m(:,J_StrucLWSOMP))^2*SNR_m/K_m);
        end
    end
end
%% Figure 14
Tc=Mrx_m*Mtx_m*4;
Effective_Rate_exh_sear=(1-(1./Tc')*(Mrx_m*Mtx_m*ones(1,length(Nrx_vec)))).*(ones(length(Tc),1)*mean(Rate_exh_sear)*ones(1,length(Nrx_vec)));
Effective_Rate_OMP=(1-(1./Tc')*(Nrx_vec.*Ntx_vec)).*(ones(length(Tc),1)*Rate_SOMP');
Effective_Rate_SLWOMP=(1-(1./Tc')*(Nrx_vec.*Ntx_vec)).*(ones(length(Tc),1)*Rate_StrucLWSOMP');
%
figure,plot(Nrx_vec.*Ntx_vec,Effective_Rate_exh_sear(1,:),'k-','linewidth',2);
hold on,plot(Nrx_vec.*Ntx_vec,Effective_Rate_OMP(1,:),'bs-','linewidth',2);
hold on,plot(Nrx_vec.*Ntx_vec,Effective_Rate_SLWOMP(1,:),'rd-','linewidth',2);
%
le=legend({'exhaustive-search','SOMP','structured LW-SOMP'});le.Interpreter='latex';le.Location='southeast';le.FontSize=14;
xlabel('Measurements ($N_{\mathrm{RX}}\times N_{\mathrm{TX}}$)','Interpreter','latex'),ylabel('$R_{\textrm{eff}}$ (b/s/Hz)','Interpreter','latex')
ax=gca;ax.TickLabelInterpreter='latex';ax.FontSize=14;
box off
grid on
saveas(gcf, '../results/fig14.png')
%
% save Res_MMV
