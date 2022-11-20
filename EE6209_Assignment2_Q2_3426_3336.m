%% Question 2
%% Projections for Numerical Method
clear;
close all;
clc;

figure(1) 
o_f = randi([1 10],5,5);
subplot(2,5,1);imshow(o_f,[0 max(o_f,[],'all')]);
title('Orginal image');



proj_angles = [0,-45,45,90];
proj = radon(o_f,proj_angles);
subplot(2,5,2);imshow(proj',[0 max(proj,[],'all')]);
title('4 Forward projections');


rec = ones(size(o_f));%intial recontructed image
subplot(2,5,4);imshow(rec,[0 max(rec,[],'all')]);
title('initial image with 0 iterations');
image_ones = ones(size(proj));
sense = iradon(image_ones,proj_angles,'none',5);
col_o_f = o_f(:);
for i=1:6
    proj_r = radon(rec,proj_angles);
    ratio = proj./(proj_r + 0.0001); %to avoid division by zero +0.0001
    subplot(2,5,3);imshow(ratio',[0 max(ratio,[],'all')]);
    title('ratio of real vs reconstructed');
    bp_ratio = iradon(ratio,proj_angles,'none',5);
    rec = (rec.*bp_ratio)./sense;
    subplot(2,5,i+4);imshow(rec,[0 max(rec,[],'all')]);
    title("Iteration no."+i+"i")

    col_rec = rec(:);
    error = sum ((col_o_f-col_rec).^2);
    rmse1(i) = sqrt(error/25);
end  
figure(2)
i = 1:1:6;
plot(i,rmse1,'-o','MarkerEdgeColor','red')
ylim([1.5 2.6])
xlim([0 7])
title('NRMS Error of 6 iterations');

ylabel('NRMSE');
xlabel('Iterations');


%% Projections for Linear Algebra


azi_angles =[0,-45,45,90];
Projection1=radon(o_f,azi_angles);
imshow(Projection1');



o_f = randi([1 10],5,5);
imagesc(o_f);
P1 =zeros(1,11);
P2 =zeros(1,11);
P1 = sum(o_f,1);
P2 = sum(o_f,2)';
% P3=zeros(1,9);
summ = zeros(1,1);
for i=1:5
    for j=1:i
        k=j;
        m=5-i+j;
        summ = summ + o_f(k,m);
    end
    P3(i)=summ;
end    

for i=1:5
    for j=1:i
        m=j;
        k=5-i+j;
        summ = summ + o_f(k,m);
    end
    P3(10-i)=summ;
    summ=0;
end 

P4=zeros(1,9);
for i=1:5
    for j=1:i
        k=j;
        m=i+1-j;
        summ = summ + o_f(k,m);
    end
    P4(i)=summ;
    summ=0;
end    

for i=1:5
    for j=1:i
        k=j;
        m=5-i+j;
        summ = summ + o_f(k,m);
    end
    P4(10-i)=summ;
    summ=0;
end 

%% Back Projection using Linear Algebra

syms X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 ;
X = [X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25 ];
clear X1 X2 X3 X4 X5 X6 X7 X8 X9 X10 X11 X12 X13 X14 X15 X16 X17 X18 X19 X20 X21 X22 X23 X24 X25;

%from 1st projection
for i= 1:5
    eqn(i) = X(i)+X(i+5)+X(i+10)+X(i+15)+X(i+20)-P1(i);
end

%from 2nd projection
for i= 1:5
    eqn(i+5) = X(5*(i-1)+1)+X(5*(i-1)+2)+X(5*(i-1)+3)+X(5*(i-1)+4)+X(5*(i-1)+5)-P2(i);
end
% from 3rd projection
for i= 1:5
    maxx=i*5;
    minn=6-i;
    now=minn;
    summ=X(now);
    for j=1:i-1
    if i>1    
    now = minn +j*(maxx-minn)/(i-1);   
    summ = summ+X(now);
    end
    end
    eqn3(1,i) = summ -P3(i) ;
end
for i= 1:4
    maxx=25-i;
    minn=5*i+1;
    now=minn;
    summ=X(now);
    for j=1:4-i
    if 4-i>0    
    now = minn +j*(maxx-minn)/(4-i);   
    summ = summ+X(now);
    end
    end
    eqn3(1,i+5) = summ-P4(i) ;
end


% from 4th projection
for i= 1:5
    maxx=5*(i-1)+1;
    minn=i;
    now=minn;
    summ=X(now);
    for j=1:i-1
    if i>1    
    now = minn +j*(maxx-minn)/(i-1);   
    summ = summ+X(now);
    end
    end
    eqn4(1,i) = summ -P4(i) ;
end
for i= 1:4
    maxx=21+i;
    minn=5*(i+1);
    now=minn;
    summ=X(now);
    for j=1:4-i
    if 4-i>0    
    now = minn +j*(maxx-minn)/(4-i);   
    summ = summ+X(now);
    end
    end
    eqn4(1,i+5) = summ -P4(i) ;
end



equations = horzcat(eqn,eqn3,eqn4);


F=@(X)[double(equations(:))];
fun=@(u) F(u(1),u(2),u(3));
x0=0;
[val1,fval1]=fsolve(F,x0),
[u2,~,fval2]=lsqnonlin(fun,x0)
    
[A,B] = equationsToMatrix(equations, X);%(1) X(2) X(3) X(4) X(5) X(6) X(7) X(8) X(9) X(10) X(11) X(12) X(13) X(14) X(15) X(16) X(17) X(18) X(19) X(20) X(21) X(22) X(23) X(24) X(25)])
Xvalues = mldivide(A,B);

sol = solve(equations,X);
X1val = sol.X1;