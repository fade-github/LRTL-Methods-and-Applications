% 主程序
clear all;
rng('shuffle') % 控制随机数生成
KnownPercentage = 50; % 图像已知比例
Image = imread('good_brother.jpg');

% 返回一个缺失数组
Image = double(Image); 
imSize = size(Image);
known = randperm(prod(imSize)/imSize(3),floor(KnownPercentage/100*(prod(imSize)/imSize(3)))); % randperm, 随机置换函数，随机生成缺失数据的下标，返回行向量
X = zeros(imSize);
X = double(ReplaceInd(X,known,Image,imSize));
imwrite(uint8(X),'Simulated_Data.jpg')

% SiLRTC算法结果
X1 = SiLRTC(Image,100,known); % 人为设置迭代次数,运行255秒 
imwrite(uint8(X1),'SiLRTC_Result.jpg');

% HaLRTC算法结果
X1 = HaLRTC(Image,50,known); % 人为设置迭代次数,运行108秒
imwrite(uint8(X1),'HaLRTC_Result.jpg');

% FaLRTC算法结果
X1 = FaLRTC(Image,50,known);  % 人为设置迭代次数,运行166秒
imwrite(uint8(X1),'FaLRTC_Result.jpg');

subplot(2,2,1)
imshow('good_brother.jpg');
title('Original Image');
subplot(2,2,2)
imshow('Simulated_Data.jpg');
title('Simulated Image');
subplot(2,2,3)
imshow('FaLRTC_Result.jpg');
title('FaLRTC Image');
subplot(2,2,4)
imshow('FaError.jpg');
title('Error Image');

            

