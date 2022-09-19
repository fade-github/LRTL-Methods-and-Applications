% num 一维位置坐标,mSize 张量大小,返回矩阵行列信息
function [a,b] = FindInd(num,mSize) % 给定Folding数组一维index找到原先二维数组的index
a = ceil(num/mSize(2)); % 向上取整
b = mSize(1) - mod(num,mSize(1)); % mod, 取余运算
end