%% 清空变量
clc,clear;
close all;
%% 数据处理
ab00=xlsread('c1.xlsx');
ab0=ab00(1:36,:);%选择36组数据作为训练数据
mu=mean(ab0);sig=std(ab0); %求均值和标准差
rr=corrcoef(ab0);   %求相关系数矩阵
ab=zscore(ab0); %数据标准化
a=ab(:,[1:30]);    %提出标准化后的自变量数据
b=ab(:,[31:end]);  %提出标准化后的因变量数据

%% 判断提出成分对的个数
[XL,YL,XS,YS,BETA,PCTVAR,MSE,stats] =plsregress(a,b);
xw=a\XS;  %求自变量提出成分的系数，每列对应一个成分，这里xw等于stats.W
yw=b\YS;  %求因变量提出成分的系数
a_0=PCTVAR(1,:);b_0=PCTVAR(2,:);% 自变量和因变量提出成分的贡献率
a_1=cumsum(a_0);b_1=cumsum(b_0);% 计算累计贡献率
i=1;%赋初始值
while ((a_1(i)<0.9)&(a_0(i)>0.1)&(b_1(i)<0.9)&(b_0(i)>0.1)) % 提取主成分的条件
    i=i+1;
end
ncomp=i;% 选取的主成分对的个数
fprintf('%d对成分分别为：\n',ncomp);% 打印主成分的信息
for i=1:ncomp
    fprintf('第%d对成分：\n',i);
    fprintf('u%d=',i);
    for k=1:size(a,2)%此处为变量x的个数
        fprintf('+(%f*x_%d)',xw(k,i),k);
    end
    fprintf('\n');
        fprintf('v%d=',i);
    for k=1:size(b,2)%此处为变量y的个数
        fprintf('+(%f*y_%d)',yw(k,i),k);
    end
    fprintf('\n');
end
%% 确定主成分后的回归分析
[XL2,YL2,XS2,YS2,BETA2,PCTVAR2,MSE2,stats2] =plsregress(a,b,ncomp);
n=size(a,2); m=size(b,2);%n是自变量的个数,m是因变量的个数
beta3(1,:)=mu(n+1:end)-mu(1:n)./sig(1:n)*BETA2([2:end],:).*sig(n+1:end); %原始数据回归方程的常数项
beta3([2:n+1],:)=(1./sig(1:n))'*sig(n+1:end).*BETA2([2:end],:); %计算原始变量x1,...,xn的系数，每一列是一个回归方程
fprintf('最后得出如下回归方程：\n')
for i=1:size(b,2)%此处为变量y的个数
    fprintf('y%d=%f',i,beta3(1,i));
    for j=1:size(a,2)%此处为变量x的个数
        fprintf('+(%f*x%d)',beta3(j+1,i),j);
    end
    fprintf('\n');
end
%% 求预测值
ab1=ab00(37:48,:);%选择12组数据作为验证数据
y1 = repmat(beta3(1,:),[12,1])+ab1(:,[1:n])*beta3([2:end],:);  %求y1,..,ym的预测值
y0 = ab1(:,end-size(y1,2)+1:end);  % 真实值

%% 贡献率画图
figure 
percent_explained1 = 100 * PCTVAR(1,:) / sum(PCTVAR(1,:));
pareto(percent_explained1);
xlabel('主成分')
ylabel('贡献率(%)')
title('主成分对自变量的贡献率')

figure
percent_explained = 100 * PCTVAR(2,:) / sum(PCTVAR(2,:));
pareto(percent_explained);
xlabel('主成分')
ylabel('贡献率(%)')
title('主成分对因变量的贡献率')
%% 绘制预测结果和真实值的对比
N = 12;% 样本个数
for i =1:size(b,2)
    yz = y0(:,i);% 真实值
    yc = y1(:,i);% 预测值
    R2 = (N*sum(yc.*yz)-sum(yc)*sum(yz))^2/((N*sum((yc).^2)-(sum(yc))^2)*(N*sum((yz).^2)-(sum(yz))^2)); %计算R方
end

[total, percentage] =clarke(yz*18,yc*18);