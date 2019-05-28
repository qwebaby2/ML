clear
close
clc
%% main function
M = 64;                                         % nums of observation point
N = 250;                                        % nums of frequency point
% creat pilot matrix
toeplitz_rows = randn(1, M);
toeplitz_columns = randn(1, N);
toeplitz_columns(1) = toeplitz_rows(1);
X = toeplitz(toeplitz_rows, toeplitz_columns);  % Toeplitz matrix

data = xlsread('C:\Users\dell\Desktop\Data.xlsx');
h = data(:, 1) + 1j * data(:, 2);               % channel value
y = X * h;                                      % observation signal
