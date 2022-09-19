function Output = MultiGrNMF(Input, Options)
X=Input.X;
H=Input.HInit;
W=Input.WInit;

GraphNum=size(Input.A,2);

tau=ones(1,GraphNum)/GraphNum;

for k=1:GraphNum
    L{k}=diag(sum(Input.A{k}))-Input.A{k};
end

   figure;
for t=1:Options.MaxIter
    A=zeros(size(Input.A{1}));
    for k=1:GraphNum
        A=A+tau(k)*Input.A{k};
    end
    U=diag(sum(A));

    H=(X*W')./(H*W*W').*H;
    W=(H'*X + Options.alpha*W*A)./(H'*H*W + Options.alpha*W*U).*W;
    
    for k=1:GraphNum
        s(k)=trace(W*L{k}*W');
    end
    
    H1=eye(GraphNum)*Options.beta/Options.alpha;
    f=s';
    Aeq=ones(1,GraphNum);beq=1;    
    LB=zeros(GraphNum,1);
    UB=ones(GraphNum,1);
    tau = quadprog(H1,f,[],[],Aeq,beq,LB,UB);
    
    clf;
    StrTitle=sprintf('t: %d',t);
    bar(tau); hold on;
    title(StrTitle);hold on;
    ylabel('tau');
    drawnow;    
end
Output.H=H;
Output.W=W;
Output.tau=tau;
end