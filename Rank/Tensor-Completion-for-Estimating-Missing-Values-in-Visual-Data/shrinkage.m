% 矩阵 X, 截断参数 t,返回压缩算子
function D = shrinkage(X,t) % 对SVD分解结果进行限制，构造新的算子D
[U,Sig,V] = svd(X,'econ'); % 矩阵X SVD分解
for i = 1:size(Sig,1)
    Sig(i,i) = max(Sig(i,i)-t,0); % 生成对角阵
end
D = U*Sig*V';
end