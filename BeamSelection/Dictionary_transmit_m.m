function [Dictionary]=Dictionary_transmit_m(sector_start,sector_end,...
    M,Phase_Shifts)
%% Generating equispaced sin(theta)
j=sqrt(-1);
ax=@(x) 1/sqrt(M)*exp(j*pi*[0:M-1]'*x); % Assuming lam/2 distance, x is sin(theta)
range=(sin(sector_end)-sin(sector_start))/M;
sampling_points=linspace(sin(sector_start)+...
    range/2,sin(sector_end)-range/2,10);
A=ax(sampling_points);
%% Quantizing based on phase shifters' resolution
AngA=angle(A);
AngAwr=wrapTo2Pi(AngA);
given_phase=[Phase_Shifts Phase_Shifts(end) + ...
    (Phase_Shifts(end)-Phase_Shifts(end-1))]; 
AngAwr_quantized = interp1(given_phase,given_phase,AngAwr,'nearest');
Dictionary=1/sqrt(M)*exp(j*AngAwr_quantized);