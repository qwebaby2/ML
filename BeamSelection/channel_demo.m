clear
close
clc
%% Parameters
iter = 100000;                                        % iterations of main loop
Nt = 64;                                            % transmit antennas in mmBS
Nr = 8;                                             % receive antennas of train
D = 500;                                            % distances between train and mmBS
beta = 3.76;                                        % path loss exponent
rho =D^(-beta/2);                                   % path loss between train and mmBS
L = 4;                                              % scattering paths
SNR_m = 6;                                          % SNR/dB
Drx=2;  % Number of phase shifters quantization bits at the Rx
Dtx=2;  % Number of phase shifters quantization bits at the Tx
v_m=1/sqrt(2*SNR_m)*(randn(Nr,Nt)+1j*randn(Nr,Nt));
sector_start = -pi/3;
sector_end = pi/3;
Phase_Shifts_tx=2*pi*[0:(2^Dtx)-1]/(2^Dtx); % The possible phase-shifts for mmWave TX
Phase_Shifts_rx=2*pi*[0:(2^Drx)-1]/(2^Drx); % The possible phase-shifts for mmWave RX
[W_m]=Dictionary_transmit_m(sector_start,sector_end,Nt,Phase_Shifts_tx); % Transmit dictionary
[Z_m]=Dictionary_receive_m(sector_start,sector_end,Nr,Phase_Shifts_rx); % Receive dictionary
WconjkronZ_m=sqrt(Nt*Nr)*kron(transpose(W_m),Z_m')';
F=Random_Codebook(Dtx,Phase_Shifts_tx,5,Nt);
Q=Random_Codebook(Drx,Phase_Shifts_rx,2,Nr);
size_recover = zeros(10,4);
data = zeros(iter, 4*L+1);
Beam_Index = zeros(iter, 2);
%% Main Loop
for kk = 1:iter
    %% Channel Simulation
    H = 0;                                              % MIMO channel
    alpha = 1/sqrt(2)*(randn(1,L)+1j*randn(1,L));       % complex path gain
    theta = 2*pi*rand(1,L);                             % AoDs
    phi = 2*pi*rand(1,L);                                % AoAs
    sintheta = sin(theta);
    sinphi = sin(phi);
    At = @(x) 1/sqrt(Nt)*exp(1j*pi*[0:Nt-1]'*x);
    Ar = @(x) 1/sqrt(Nr)*exp(1j*pi*[0:Nr-1]'*x);
    %
    for jj = 1:L
        H = H + alpha(jj) * At(sintheta(jj)) * Ar(sinphi(jj))';
    end
    %
    H = sqrt(Nt*Nr)*H;
    %% Exhaustive search at infinite SNR (no noise) to get the best beam-direction
    ZHW_exh_sear_no_noise_m=W_m'*H*Z_m;
    ZHW_exh_sear_no_noise_abs_m=mean(abs(ZHW_exh_sear_no_noise_m).^2,3);
    % the sparsity of ZHW
    % bar3(ZHW_exh_sear_no_noise_abs_m);
    [~,Beam_indices_sorted]=sort(abs(ZHW_exh_sear_no_noise_abs_m(:)),'descend');
    %
    Best_beam_index=Beam_indices_sorted(1);
    %     %% LW-OMP
    %     ZHW_LWOMP_m = F'*H*Q;
    %     FtranskronQherm=kron(transpose(F),Q');
    %     Psi=FtranskronQherm*WconjkronZ_m;
    %     onebyPsi_col_norms=sqrt(diag(1./diag(Psi'*Psi)));
    %     Psi=Psi*onebyPsi_col_norms;
    %     normlize Psi
    %     Result_LWOMP=LWOMP(ZHW_LWOMP_m(:),Psi);
    %
    [I_Exha,J_Exha] = ind2sub(size(size_recover),Best_beam_index);
    %
    %% Generate Training Data
    data(kk,:) = [real(alpha), imag(alpha), theta, phi, Best_beam_index];
    Beam_Index(kk,:) = [I_Exha,J_Exha];
end
csvwrite('data.csv',data)