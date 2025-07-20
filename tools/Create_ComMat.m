function K = Create_ComMat(p,q)
K=[];
for i=1:p
    Ki=[];
    for j=1:q
        Mat=zeros(q,p);
        Mat(j,i)=1;
        Ki=[Ki Mat];
    end
    K=[K; Ki];
end
end
