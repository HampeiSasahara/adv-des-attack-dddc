function dF = AutoDiff(F,X0)

h0=1e-8;
F0=F(X0);
nF0=norm(F0,'fro');
if nF0<1e-0
    h=h0;
else
    h=h0*nF0;
end
[n,m]=size(F0);
[p,q]=size(X0);
dF=zeros(n*m,p*q);

for j=1:q
    for k=1:p
        Ejk=zeros(p,q); Ejk(k,j)=h;
        DXjkp=X0+Ejk; DXjkn=X0-Ejk;
        DFjk=(F(DXjkp)-F(DXjkn))/(2*h);
        dF(:,(j-1)*p+k)=DFjk(:);
    end
end

%{
for i=1:m % f_i
    DFi=zeros(n,p*q);
    %fi0=F0(:,i);
    ei=zeros(m,1);
    ei(i)=1;
    fi=@(X) F(X)*ei;
    for j=1:q % x_j
        DFij=zeros(n,p);
        for k=1:p
            ek=zeros(p,q);
            ek(k,j)=1;
            DXjk=X0+h*ek; DXjkn=X0-h*ek;
            %dfidxjk=(fi(DXjk)-fi0)/h;
            dfidxjk=(fi(DXjk)-fi(DXjkn))/(2*h);
            DFij(:,k)=dfidxjk;
        end
        DFi(:,(j-1)*p+1:j*p)=DFij;
    end
    DF((i-1)*n+1:i*n,:)=DFi;
end
%}

%{
DF=[];
for i=1:m % f_i
    DFi=[];
    fi0=F0(:,i);
    ei=zeros(m,1);
    ei(i)=1;
    fi=@(X) F(X)*ei;
    for j=1:q % x_j
        DFij=zeros(n,p);
        for k=1:p
            ek=zeros(p,q);
            ek(k,j)=1;
            DXjk=X0+h*ek;
            DXjkn=X0-h*ek;
            %dfidxjk=(fi(DXjk)-fi0)/h;
            dfidxjk=(fi(DXjk)-fi(DXjkn))/(2*h);
            DFij(:,k)=dfidxjk;
        end
        DFi=[DFi DFij];
    end
    DF=[DF;DFi];
end
%}

end

