% x = soft(x + Phi'*(y-Phi*x),lambda); update x


function [ x ] = IST_Basic( y,Phi,lambda,epsilon,loopmax )


    if nargin < 5    
        loopmax = 200;    
    end    
    if nargin < 4      
        epsilon = 1e-20;      
    end     
    if nargin < 3      
        lambda = 0.1*max(abs(Phi'*y));     
    end     
    [y_rows,y_columns] = size(y);      
    if y_rows<y_columns      
        y = y';% y should be a column vector      
    end 
    %soft = @(x,T) sign(x).*max(abs(x) - T,0);      %%%%%%Òþ²Øº¯Êý
    n = size(Phi,2);
    x = zeros(n,1);% Initialize x=0
    f = 0.5*(y-Phi*x)'*(y-Phi*x)+lambda*sum(abs(x));
    loop = 0;
    fprintf('\n');
    while 1
        x_pre = x;
        x = soft(x_pre + Phi'*(y-Phi*x),lambda);% update x
        loop = loop + 1;
        f_pre = f;
        f = 0.5*(y-Phi*x)'*(y-Phi*x)+lambda*sum(abs(x));
        if abs(f-f_pre)/f_pre<epsilon
            fprintf('abs(f-f_pre)/f_pre<%f\n',epsilon);
            fprintf('IST loop is %d\n',loop);
            break;
        end
        if loop >= loopmax
            fprintf('loop > %d\n',loopmax);
            break;
        end
        if norm(x-x_pre)<epsilon
            fprintf('norm(x-x_pre)<%f\n',epsilon);
            fprintf('IST loop is %d\n',loop);
            break;
        end
    end  
end   
