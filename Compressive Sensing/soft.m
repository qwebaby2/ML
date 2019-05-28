

function [ x ] = soft (b,T)
    x = sign(b).*max(abs(b) - T,0);
end


% r = 1;
% x = -3:0.1:3;
% S1 = (max(abs(x) - r,0) ./ (max(abs(x) - r,0) + r)) .* x;
% S2 = sign(x) .* max(abs(x) - r,0);
% subplot(1,2,1)
% plot(x,S1,'b');
% subplot(1,2,2)
% plot(x,S2,'r');
