%SMO algorithm solve the SVM optimization question;
%build the inner kernel at first, since the optimization algorithm is the
%easiest.
%then develop to the selection method;
%-------------------20140204----------------------------------------------
%revise the stop judgement:whether the alpha parameters and related points
%meet the KKT conditions. no obvious improvement.
%17:16 try another method: C-SVM;
%at first, C=100;
%and at the same time, I need efxt;
%-------------------------------------------------------------------------

%{

上面是对这个问题的大致说明，我的方法是由计算内核（核函数）开始，逐步往上构建整个算法，由于运算涉及到多次使用核函数，因此如何简化核函数、避免重复运算，成为提升效率、优化代码的关键；

当然，在这之前，我也写下了整个程序流程图，这让我知道每一步的意义，也知道在编程的过程中，代码处于什么位置，又将向何处去。

%}


clear all;                      
close all;


%关闭图形和清除内存中各参数；

st = cputime;                   %record current time; 这个是为了统计整体计算时间；


%-------------------------------------------------
%    set the initial value;  
%-------------------------------------------------


N=300;                          %the training data set size;                       300个训练点；
NN=2000;                        %the number of validation data set;     2000个验证点；
NS=N;                           %number of SVs.                                         支持向量的数目，一开始假设所有数据点都是支持向量；


%------------------------------------------------
%produce the training set                                
%------------------------------------------------                                


    fprintf('Producing the training set...\n');
    fprintf('-----------------------------------------\n');

%运算中我们需要知道程序进行到哪一步了，这就需要在代码中嵌入这些指示性文字；



r=10;                               %半径=10；
wth=6;                           %宽度=6；
dis=-6.5;                             %距离=-6；
si=0.085;


%构建双月形状，设定半径、宽度、距离、以及核函数的方差，这个方差不是一蹴而就的，可以在第一次运算后，通过统计支持向量之间的最大距离和支持向量数目运算得出；

for i=1:N+NN;                       %at the beginning of every round, to produce the training set in random; 看了Yanbo Xue的程序，他是一次性产生3000个数据点，
                                    %然后用这三千个点中的一千个反复训练，训练50个回合；先试试他的方法；
                                    %make the rand vectors X;
    if  rand>0.5;                   %make sure the points occur in random. and half in half for each class;
        theta=pi*rand;              %set the angle;
        lenth=r+wth*(rand-0.5);     %the range;
        h(i)=lenth*cos(theta);      %the X axis;
        ve(i)=lenth*sin(theta);     %the Y axis;
        y(i)=-1;                     %the expected response; class 1;
    else
        beta=-pi*rand;
        len_2=r+wth*(rand-0.5);
        h(i)=len_2*cos(beta)+r;
        ve(i)=len_2*sin(beta)-dis;
        y(i)=1;                   %class 2;
    end
      
end



%这里是产生双月形状的代码，主要到数量是N+NN，也就是说训练集和验证集产生于同一个模型；
 
    fprintf('Normalization...\n');
    fprintf('-----------------------------------------\n');
    %the input data need to deal with;去均值和正则化？
    % and use the random data.
    %首先去掉均值；


    miu_x=mean(h);                 %取得x1的均值；
    miu_y=mean(ve);                %取得x2的均值；                
    h=h-miu_x*ones(1,NN+N);           %减去均值；相当于减去了直流分量；
    ve=ve-miu_y*ones(1,NN+N);         %减去均值；
    %正则化？高斯化？让取值落在[0,1]区间?
    

    max_h=max(abs(h));             %取得绝对值的最大值；
    max_ve=max(abs(ve));        
    h=h./max_h;                    %除以这个数；
    ve=ve./max_ve;
    tr_set=[h;ve;y];               %组合成一个训练集，响应没有变化；
    %X=[h;ve];                      %只取数据点，不取标号，但不打乱顺序，这样可以恢复出响应来；
    pos = find(y == 1);
    neg = find(y == -1);
    hold on;
    plot(tr_set(2,pos),tr_set(1,pos),'k.');
    plot(tr_set(2,neg),tr_set(1,neg),'b.');
    seq_tr_set=randperm(N + NN);         %重新洗牌；                     
    seq_tr_set = seq_tr_set(1,1:N);
    plot(tr_set(2,seq_tr_set),tr_set(1,seq_tr_set),'r.');
    hold off;
%---------------------------------------------------------------------------------------------------------------

%注意到这里，randperm的意思是打乱次序，也就是随机取点，不受产生顺序的限制；

%---------------------------------------------------------------------------------------------------------------
    mix_tr_set=tr_set(:,seq_tr_set);  %注意到从2300个点中取得300个；
    h=mix_tr_set(1,:);                %第一行是x轴数据，或x1；
    ve=mix_tr_set(2,:);              
    y=mix_tr_set(3,:);                %注意到响应是跟着数据点走的；
    %tr_set_1=[h;ve;y];       
    X=[h;ve];                    %某些时候只需要用到两列，但是响应还是随着数据点走的；
    



%以上步骤是把产生的数据归一化，这样便于处理；这样数据已经产生完毕；
    
    fprintf('Now we start the SVM by SMO algorithm ...\n');
    fprintf('-----------------------------------------\n');
    

%从这里开始，我们准备构建SMO算法；
    %------------------------------------------------------
    %set the SVM initial value;
    %------------------------------------------------------


    a=zeros(1,N);                   %the first initial value; 这就是著名的拉格朗日算子，每个数据点都有一个；
    b=0;                            %bias.                                   偏置；
    C=100;                          %the boundry;                  惩罚算子；  
    FX=zeros(1,N);                  %f(x(i)) i=1toN, since {a(i)}=0, so f(x)=0; f(x)=w'*x-b; or f(x)=SUM ai*yi*K(xi,xj)-b=0; f(xi)函数，一般有f(xi)=w'×xi-b；这里是f(xi)=sum ai*yi*K(xi,xj)-b;
    E=FX-y;                         %y: the expected response;      这是函数f（xi）与期望响应（或者叫做标号）的差异；
                                    %??E should be followed by training set;
    efxt=a./C;                      %efxt is the gap distance related parameters;                                这个是对于那些算子=C的点对应的间隔；
    
    eps=0.007;                      %threshold for stopping judgement;                                 阈值，（原目标-对偶目标）/（原目标+1）的比值对应停机条件，小于这个阈值，运算停止；
    times=0;                        %at the start, times=0;the number of                                   外循环运算多少次了；
                                    %external loops;
    presv=0;                        %at the beginning, suppose there is no                              相当于flag，表示是否找到不满足KKT条件的支持向量；
                                    %support vector found;
    Gram=eye(N,N);                 %build the gram matrix, for recording the calculation results;  

                                                   %Gram矩阵，我们不需要运算矩阵中所有的点，但是对于已经计算了的点，我们无需重复计算；注意到采用的是eye矩阵，这样对角数据全为1；

                                                   %这样有两个用意，第一，数据确实为1，第二，其他数据为0，这样可以判断这个位置是否运算过了；
                                   
                                   
    ot=0;                                      %指示是否运算数目太多，迭代时间太长了；                    
    totaltimes=0;                   %how many times of calculation operated?   一共运行多少次两算子计算了；
    in=1;                           %set the initial value of in;in start from 1;                为了区分违反KKT条件的三种算子，有三种指示，a=0；a=C；C>a>0;
    ic=1;
    i0=1;
    
while(1)                           %the first loop, when ratio<eps, loop stop;                  开始外循环，采用启发式方法选择第一个点；
%select the first point;
%--------------------------------------------------------------------------
if (times==0)                      %first selection or calculation;                                     如果算法初次运行，随意选择一个点，这里选择点1；
    i1=1;                          %select the point in random, so chose 1;  
   
else                               %if this is not the first selection:
    %----------------------------------------------------------------------
    %when other selection: 
    %1, choose the (0,C) alpha which break the KKT conditions;
    %2, then alpha=0 and alpha=C points which break the KKT conditions;
    %3, if alpha choosen is in the last of the squene, re-start from 1;
    %----------------------------------------------------------------------
    
    while (in<=N)                  %search in all (0,C)alpha;                                          以后的运算，首先考虑那些界内的、不满足KKT条件的算子；
        if (a(in)>0) && (a(in)<C)  %here, C is 100.
            if (y(in)*FX(in)>1.01)||(y(in)*FX(in)<0.99)   %meet the KKT condition or not:               这里不能用不等于1来表示违反KKT条件，而是用一个容许范围，取0.99～1.01；
                i1=in;             %if break the KKT rule, i1=i;  
                presv=1;           %we found a SV already;                                                                         如果满足，记录下当前序号，并改变flag标志；
            end
        end        
        in=in+1;                   %if we don't found the break KKT condition SV, continue...               循环进行；
        if (presv)                 %we've found SV and not exceed the NS times;                                   这里采用的是do while结构，所以采用了一个判决语句，满足条件时退出；
            break;
        end                
    end
    
    if presv==0                    %if we don't found SV which break KKT condition, or, the times_i1 out of NS times;        

                                             %如果界内的算子都满足KKT条件，考虑界上的算子（a=0或a=C）        
        in=1;                            %及时对界内算子的序号复位，这里复位为1； 
        while (ic<=N)
                                   %if the alpha on the boundry,i.e. a=0 or C;
            if a(ic)==C                 %这里首先考虑a=C的算子，因为这也是支持向量；
                if (y(ic)*FX(ic)>1.01)                        %注意到这里也采用了容许范围；
                    i1=ic;         %if break the KKT rule, i1=i;                    
                    presv=1;       %we found a SV already;                
                end
            end
            ic=ic+1;
            if (presv)
                break;
            end            
        end
    end
    
    if presv==0
        in=1;
        
        while (i0<=N)
                                   %if the alpha on the boundry,i.e. a=0 or C;
            if a(i0)==0     %现在考虑a=0的算子，这些算子是非常多的；这些就不是支持向量了；
                if (y(i0)*FX(i0)<0.99)
                    i1=i0;         %if break the KKT rule, i1=i;                    
                    presv=1;       %we found a SV already;                
                end
            end
            i0=i0+1;
            if (presv)
                break;
            end           
        end
    end
      
    
    if (presv==0)                  %if we didn't found a SV which break KKT;        
        in=1;
        ic=1;
        i0=1;
        i1=floor(rand*N)+1;      %如果所有的支持向量都满足KKT条件，那么随意选择一个算子；
    end
    
    presv=0;                       %back to the initial value;i.e.no sv found;           标志位清零；
    
end
%--------------------------------------------------------------------------       
    
  
%{
  --------------------------------------------------------
  now I have the first point, and search the second point.

第一个算子找到之后，开始第二个算子；
  --------------------------------------------------------
%}
                           
max_i=1;                 
times_2=0;                            %how many times of internal loop?
i2old=i1;                             %the initial value of i2;
while(1)                              %the important question is how to
                                      %stop the loop and this is internal
                                      %loop, too.


if (times==0)                         %at the first time, most points have                初次运算，所有的算子都要更新一下；
                                      %a(i)=0;
   if i2old<N                        %i2old from i1:N; i.e. 1:N.
        i2=i2old+1;                   %i2 from 2:N+1;
        i2old=i2;                     %remember the choice;
   end
   
else                                  %at the other loops,there should be              在之后的运算中，我们要寻找那些|E1-E2|最大的点；
                                      %some points a(i)>0;
    min_E=E(i1);                          
    max_E=E(i1);
  
    if (E(i1)>=0)
        for i=1:N
            if (a(i)>0) && (a(i)<C)
                if (E(i)<min_E)%&&(i~=i2old_x)
                    min_E=E(i);
                    i2=i; 
                    %i2old_x=i2;
                end
            end            
        end        
    else
        for i=1:N
            if (a(i)>0) && (a(i)<C)
                if (E(i)>max_E)%&&(i~=i2old_m)
                    max_E=E(i);
                    i2=i;
                    %i2old_m=i2;
                end
            end            
        end
    end
  end
%end
    times_2=times_2+1;      %totally N for times=0; 
    if times_2 == 300
        N = N;
    end
if times_2>=NS            %finish one complete loop;         
   break;
end
    totaltimes=totaltimes+1;         %total times of calculation;
    
    %---------------------------------------------------------------
    %now I have two points, start the calculation;
    %---------------------------------------------------------------

%now we get the i2, so we could start the optimization;
       
%select the second point;


%SMO optimization algorithm;            从现在开始就是SMO的运算了，每次更新两个点，然后将更新后的Fx1, Fx2上传；

%{
possible parameters;
SV(i); collect all support vectors into one group and calculate the E1&E2
using these elements;
x1,x2;
x_sv(i) and y_sv(i) compare to the SV(i); 
---------------------------------------------
maybe we don't need to collect the SV group.
we just choose the {a(i)>0}is okay,
when the algorithm stops, we could update the
SVs and give the final decision hyper plane.
---------------------------------------------
a1,a2;
a2new,a2new-unc,a2old;a1new,a1old;
y1,y2;
E1,E2;
L,H;
K11,K22,K12;
%}


%and we need the kernel function;
%{
    K(x,z)=exp(-1/(2*sigma^2)*norm(x-z)^2);
%}
%--------------------------------------------------------------------------
%
%RUN the SMO algorithm
%
%--------------------------------------------------------------------------
x1=X(:,i1);                               %这里都比较容易看懂了，只解释比较隐晦的部分；
x2=X(:,i2);
y1=y(i1);
y2=y(i2);
a1old=a(i1);
a2old=a(i2);
%when C=inf, the limits should be adjusted as below.
if(y1~=y2)
    L=max(0,a2old-a1old);
    H=min(C,C-a1old+a2old);
else
    L=max(0,a1old+a2old-C);
    H=min(C,a1old+a2old);
end


if times~=0 && L>=H                         %first time, we should let L=H=0;
    break;
end


K11=1;%K(x1,x1); for Guass kernel function;
K22=1;%K(x2,x2); ditto;
if Gram(i1,i2)~=0                         %注意到这里的选择机制，一旦Gram中的元素不等于0，就说明这两个点已经运算过，可以跳过了；
    K12=Gram(i1,i2);
else
    K12=K(x1,x2,si);
    Gram(i1,i2)=K12;    
end


%we need two loop for E1&E2 calculation;
%l, the number of support vectors; at the first, this value is 0, then grow
%slowly.
%I need the mapping between support vectors and original points. 
%set N(i)=number of the original points, that means the pointer. such as
%N(1)=30, means the first support vector is the 30th points.
%so the pointer of original points is important.
%and I want to know the verse # of the SV, i.e.O(30)=1, means the 30th
%point in the original points, is the first SV.




E1=E(i1);                        %this value got from memory;
E2=E(i2);                        %ditto;


k=K11+K22-2*K12;                 %parameter couldn't be 0 when it worked as divider.


if k==0                          %the possibility is rare.
    k=0.01;
end
%?? how about when k=0???!!!!!----------------------
%---------------------------------------------------
a2new=a2old+y2*(E1-E2)/k;                    %这就是算子的更新运算了；


if (a2new>H)
    a2new=H;
else
    if(a2new<L)
        a2new=L;
    end    
end


a1new=a1old+y1*y2*(a2old-a2new);
a1new=max(0,a1new);                               %算子a1的更新；?
%{
if abs(a2new-a2old)<(eps*(a2new+a2old-eps))          %if the difference is little, jump out of the loop;
    break;
end
%}
a(i1)=a1new;                           %now we could update the{a(i)}i=1toN
a(i2)=a2new;                           %update the {a(i)};
%now we starts update the bias;
%--------------------------------------------------------------------------
bold=b;
a1e=y1*(a1new-a1old);
a2e=y2*(a2new-a2old);
b1new=E1+a1e+a2e*K12+bold;                        %偏置的更新；
b2new=E2+a1e*K12+a2e+bold;


if a1new>0 && a1new<C                  %if a1new is in the bounds;
    b=b1new;
else
    if a2new>0 && a2new<C
        b=b2new;
    else 
        if (a1new==0||a1new==C)&&(a2new==0||a2new==C)&&(L~=H)
            b=(b1new+b2new)/2;
        end
    end
end
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%**************************************************************************
%--------------------------------------------------------------------------
%更新F(x1),F(x2);
FX1=0;
FX2=0;
for (i=1:N)
    if a(i)>0
        if Gram(i,i1)~=0;
             Ki1=Gram(i,i1);
        else
            Ki1=K(X(:,i),x1,si);
            Gram(i,i1)=Ki1;
        end
        FX1=FX1+a(i)*y(i)*Ki1;
        if Gram(i,i2)~=0;
            Ki2=Gram(i,i2);
        else
            Ki2=K(X(:,i),x2,si);
            Gram(i,i2)=Ki2;
        end
        FX2=FX2+a(i)*y(i)*Ki2;
    end
end
FX(i1)=FX1-b;                           %FX=SUM ai*yi*K(xi,xj)-b;
FX(i2)=FX2-b;
E(i1)=FX(i1)-y1;                        %store the E(i) into the E matrix;
E(i2)=FX(i2)-y2;


%N, the number of training set;


%now we calculate the stop formula and judge whether the algorithm should
%stop.
end
    %---------------------------------------------------------------
    %one internal loops complete, times++ for external loop;
    %---------------------------------------------------------------
%C=max(a)+1;
%--------------------------------------------------------------------------
% I hope that calculation could be complete less than (n-1)*(n-2)/2 times;
%--------------------------------------------------------------------------
    
    if totaltimes>100*N                                           %运算时间过长、迭代次数太多，算法自动终止；
        fprintf('how many knives?\n');
        fprintf('-----------------------------------------\n');
        break;
    end
    
times=times+1;                        %times++, i.e. outer loop ++;  

%do we need update the E(i) and FX(i)? I think so. FX(i) is neccesary, E(i)
%don't need.
%now I hesitate.


%这里开始记录所有的支持向量，并更新所有的F(xi)，i=1:N;
i=1;
l=0;
while i<=N
    if i == 300
        N = N;
    end
    if a(i)>0
        l=l+1;
        SV(l)=a(i);
        y_sv(l)=y(i);
        x_sv(:,l)=X(:,i);
        ptr(l)=i;                   %remember the pointer;
    end
    i=i+1;
end
NS=l;                               %the number of support vectors.
%{
if NS<=(N/10)
    p=0;
    for j=1:NS;
        for i=1:NS;
            d_max=norm(x_sv(:,i)-x_sv(:,j));
            if d_max>p;
                p=d_max;
            end
        end
    end
    d_max=p;
    si=d_max^2/(2*NS);
end
%}
%{
lold=1;
while (a(lold)==C)
    lold=lold+1;
end
FV=0;
for l=1:NS
    FV=FV+SV(l)*y_sv(l)*K(x_sv(:,l),x_sv(:,lold));
end
b=FV-y_sv(lold);
lold=lold+1;
if lold>NS
    lold=1;
end
%}
FX=zeros(1,N);
for i=1:N;
    for l=1:NS
        if Gram(ptr(l),i)~=0;
            Kli=Gram(ptr(l),i);
        else
            Kli=K(x_sv(:,l),X(:,i),si);
            Gram(ptr(l),i)=Kli;
        end
        FX(i)=FX(i)+SV(l)*y_sv(l)*Kli;                          
    end
end
sv_x=x_sv;                         %this vector for demostation.
x_sv=zeros(2,NS);                  %clear past data for preparation.
FX=FX-b*ones(1,N);
E=FX-y;


for i=1:N    
    efxt(i)=max(0,1-y(i)*FX(i));
end


%-------------------------------------
%Now we revise the stop criteria.这里是关于停机条件的运算；
%-------------------------------------


W=0;
AAYYK=0;
for i=1:N
    if (a(i)>0)
        for j=1:N
            if (a(j)>0)
                if Gram(i,j)~=0
                    Kij=Gram(i,j);
                else
                    Kij=K(X(:,i),X(:,j),si);
                    Gram(i,j)=Kij;
                end
                AAYYK=AAYYK+a(i)*a(j)*y(i)*y(j)*Kij;
            end
        end
    end
end


suma=sum(a);
sume=C*sum(efxt);
W=suma-0.5*AAYYK;


ratio=((suma-2*W+sume)/(suma-W+sume+1));             %比值接近于0才行的；


fprintf('Now we can see the difference between original objective and dual one...\n');
fprintf('Ratio = %f\n',ratio);


if (ratio<eps)
    break;                          %if ratio meet the stop threshold, loop
                                    %stops.
end


end                                 %this end compare to while(1);


%now the SV calculation is complete.
%and there should be some SVs occured, now I can collect the SV group, 
%check the validation data set and draw the decision hyper plane.

%运算结束后，开始记录支持向量；

i=1;
l=0;
while i<=N
    if a(i)>0
        l=l+1;
        SV(l)=a(i);
        y_sv(l)=y(i);
        x_sv(:,l)=X(:,i);   
        ptr(l)=i;
    end
    if i == 299
        N = N;
    end
    i=i+1;
end
NS=l;                               %the number of support vectors.


%now we should update the value of b.!!!!!!更新偏置；


lold=floor(rand*NS)+1;
while (a(lold)==C)
    lold=lold+1;
end
FV=0;
for l=1:NS
    if Gram(ptr(l),ptr(lold))~=0
        Kllo=Gram(ptr(l),ptr(lold));
    else
        Kllo=K(x_sv(:,l),x_sv(:,lold),si);
        Gram(ptr(l),ptr(lold))=Kllo;
    end
    FV=FV+SV(l)*y_sv(l)*Kllo;
end
b=FV-y_sv(lold);


%now I collect all support vectors and related expected response, input;


%now I can check the validation set.验证；


% this module for producing the validation set.
%------------------------------------------------------------------
    seq_tr_set=randperm(NN);         %重新洗牌；
    mix_tr_set=tr_set(:,seq_tr_set);  %注意到从3000个点中取得1000个；
    h=mix_tr_set(1,:);                %第一行是x轴数据，或x1；
    ve=mix_tr_set(2,:);              
    y=mix_tr_set(3,:);                %注意到响应是跟着数据点走的；
    tr_set=[h;ve;y];       
    X=[h;ve];                    %某些时候只需要用到两列，但是响应还是随着数据点走的；
%------------------------------------------------------------------
% and do some normlization for this set.


err=0;


for j=1:NN
    FV=0;
    for l=1:NS        
        FV=FV+SV(l)*y_sv(l)*K(X(:,j),x_sv(:,l),si);
    end
    FV=FV-b;
    if (y(j)*FV)>=1
        FV=0;
    else
        if (y(j)*FV)<0
            err=err+1;
            FV=0;
        end
    end
end


%now the validation is complete.


%start the decision curve drawing;画出决策超平面；


fprintf('Draw the judgement curve ...\n');
fprintf('------------------------------------\n');


figure;
plot(h,ve,'.');
hold on;
plot(x_sv(1,:),x_sv(2,:),'r*');


for il=1:400;    
    xl(il)=-1+il/200; 
    ul=10;
    pl=2;
    for jl=1:400;
        yl(jl)=-1+jl/200;
        Cur=[xl(il) yl(jl)]';
        oo=0;
        for l=1:NS
            oo=oo+SV(l)*y_sv(l)*K(Cur,x_sv(:,l),si);
        end          
        oo=oo-b;
        
        zl=abs(oo);
        if  zl<ul;
            ul=zl;
            pl=yl(jl);
        end
    end
    pl_l(il)=pl;
end
hold on;
plot(xl,pl_l,'r.');
fprintf('run time = %4.2f seconds\n',cputime-st); %统计时间，判断运算效率；

%{
l=0;
for i=1:N
    if y(i)*FX(i)==1
        l=l+1;
    end
end
%}




