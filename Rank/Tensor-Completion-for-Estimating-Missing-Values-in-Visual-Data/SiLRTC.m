% 原图像 Image,人为设置迭代次数 K,缺失机制 known
function X2 = SiLRTC(Image,K,known)
tic; % 计时器
imSize = size(Image); % 计算数组各维度大小，返回行向量

% Corrupting Image
X = zeros(imSize); % 创建全零数组，全部缺失
X = double(ReplaceInd(X,known,Image,imSize)); % ReplaceInd，用真实数值替换全零矩阵，保证Omega指标下数组数值相同，构造缺失数组X 

% Create random alphas and betas
b = abs(randn(3,1))/200; % 算法中的超参数设置,randn函数生成随机矩阵
a = abs(randn(3,1));
a = a./sum(a); % 标准化

% Block Coordinate Decent Algorithm
for k = 1:K % Iteration step
    M = zeros(imSize); % initialize Mi's
    for i = 1:3 
        % Create Mi's multiply by appropriate beta value fold them and add into tensor to get Mi values 
        % M = M + Mfold(b(i)*shrinkage(unfold(X,i),a(i)/b(i)),i,imSize);
        M = M + Mfold(b(i)*shrinkage(ten2mat(X,i),a(i)/b(i)),i,imSize); % The closed form of convex problem,是缺失数据张量X的函数
    end
    % Divide by sum of betas
    M = M./sum(b); % 点除，矩阵元素对应相除，对于分母是数值来说与除法没有区别，更新缺失值下的张量X

    % Update indices that we know from Image into M and set X equal to M
    M = ReplaceInd(M,known,Image,imSize); % Image 对应的是真实数据，数值保持不变，M是不断进行迭代更新出来的(Convex closed form)
    X = M; % BCD 分块算法的更新过程
end
% Output Absolute Value Error Image
imwrite(uint8(abs(double(Image)-double(X))),'SiError.jpg')
X2 = X;
toc;
end