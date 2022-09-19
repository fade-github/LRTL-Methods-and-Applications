% 原图像 Image,人为设置迭代次数 K,缺失机制 known
function X2 = FaLRTC(Image,K,known)
tic;

imSize = size(Image);
X = zeros(imSize);
X = double(ReplaceInd(X,known,Image,imSize));

% Create needed constants and holder variables
a = abs(randn(3,1));
a = a./sum(a);
C = 0.5;
u = 10^5;
ui = a./u;
Z = X;
W = X;
B = 0;
L = sum(ui);

for k = 0:K
    while true
        Theta = (1+sqrt(1+4*L*B)) / (2*L);
        W = Theta/(B+Theta) * Z + B/(B+Theta) * X;

        % Compute f(x), f(w) and f'(w).
        dfw = zeros(imSize);
        fx = 0;
        fw = 0;

        for i = 1:3
            [Sig,~] = Trimming(ten2mat(X,i),ui(i)/a(i));
            fx = fx + sum(Sig(:));
            [Sig,T] = Trimming(ten2mat(W,i),ui(i)/a(i));
            fw = fw + sum(Sig(:));
            dfw = dfw + Mfold(a(i)^2/ui(i)*T,i,imSize); % 求导运算规则
        end
        dfw = ReplaceInd(dfw,known,zeros(imSize),imSize); % 每次循环梯度的更新，缺失数据补全的核心步骤

        % If step size is too low then break
        if fx <= fw - sum(dfw(:).^2)/L
            break;
        end

        % Solve for new f(X) using modified W matrix
        Xp = W - dfw/L; % 这三者都会进行参数更新，使得张量Xp更新 (一次循环)
        fxp = 0;
        for i = 1:3
            [Sig,~] = Trimming(ten2mat(Xp,i),ui(i)/a(i));
            fxp = fxp + sum(Sig(:));
        end

        % If step size is not too large update and break
        if fxp <= fw - sum(dfw(:).^2)/2/L
            X = Xp; % 本次循环结束，返回补全的张量X
            break;
        end

        % Update L and rerun 循环结束的最终条件
        L = L/C; 
        if L == Inf
            disp('Error')
            imwrite(uint8(X),'FaFinal.jpg') % If L is too large throw error and exit
            exit
        end
        
    end
    
    % Update Z and B 下一次循环开始
    Z = Z - Theta*dfw;
    B = B + Theta;
    
end
% Output Image
% imwrite(uint8(X),'FaFinal.jpg')
%Output Absolute Value Error Image
imwrite(uint8(abs(double(Image)-double(X))),'FaError.jpg')
X2 = X;

toc;
end
