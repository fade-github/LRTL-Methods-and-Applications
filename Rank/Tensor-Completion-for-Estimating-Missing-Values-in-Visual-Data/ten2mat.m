% 可展目标张量 ten,按照维度 k 展开,返回矩阵
function mat = ten2mat(ten,k) % 张量展开为矩阵
dim = size(ten);
% Tensor unfolding for tensor with size of n1*n2*...*nd
% dim = (n1,n2,...,nd), k=1,2,...,d
redim = permute(ten,[k,1:k-1,k+1:length(dim)]); % 重新定义张量ten的维度，张量数值并未发生改变
mat = reshape(redim,dim(k),[]);