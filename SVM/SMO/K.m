function f = K(x1, x2)
sigma = 3;
f = exp(-1/(2*sigma^2)*norm(x1-x2)^2);