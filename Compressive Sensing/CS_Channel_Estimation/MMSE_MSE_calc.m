%Function Declaration:
function ms_error=MMSE_MSE_calc(X,H,Y,Rgg,variance,N,L);
%This function generates mean squared error for the the MMSE estimator..
%EVALUATION OF Hmmse
%Hmmse=F*Rgg*inv(Rgy)*Y;
F=fft(eye(N,L));
I=eye(N,N);
Rgy=Rgg * F'* X';
Ryy=X * F * Rgg * F' *X' + variance * I;
for i=1:N
    yy(i,i)=Y(i);
end
Gmmse=Rgy * inv(Ryy)* Y;
Hmmse=F*Gmmse;

ms_error_mat=mean(((abs(H)-abs(Hmmse))/abs(H)).^2);
for i=1:N
    if(ms_error_mat(i)~=0)
        ms_error=ms_error_mat(i);
    end
end