% 矩阵 X, 截断参数 t,返回截断算子T和特征根
function [Sig,T] = Trimming(X,t) % 对SVD分解结果进行限制
[U,Sig,V] = svd(X,'econ');
SigT = Sig;
for i = 1:size(SigT,1)
    SigT(i,i) = min(SigT(i,i),t);
end
T = U*SigT*V';
end
