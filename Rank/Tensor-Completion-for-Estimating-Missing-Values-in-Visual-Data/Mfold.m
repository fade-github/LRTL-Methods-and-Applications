% 重组矩阵 a，任意张量维度 i,目标张量的维数 OrDim,返回一个重塑的张量(folding)
function c = Mfold(a,i,OrDim) % 对矩阵 a 按照 i 维度组合成三维张量
c = [];
for j = 1:OrDim(i) 
    switch i % switch-case 语句，表示按照不同维度重塑
        case 1
            c(j,:,:) = reshape(a(j,:),OrDim(2),OrDim(3));
        case 2
            c(:,j,:) = reshape(a(j,:),OrDim(1),OrDim(3));
        case 3
            c(:,:,j) = reshape(a(j,:),OrDim(1),OrDim(2));
    end
end
end