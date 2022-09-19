% 张量 X,随机缺失机制 Known,缺失数据张量 Image,张量大小 imSize,返回X更新后的数组
function mOut = ReplaceInd(X,known,Image,imSize) % 对数组X进行随机替换(数组更新)，随机性由参数向量known决定
for i = 1:length(known)
    [in1,in2] = FindInd(known(i),imSize);
    X(in1,in2,:) = double(Image(in1,in2,:)); % 人为缺失，
end
mOut = X;
end