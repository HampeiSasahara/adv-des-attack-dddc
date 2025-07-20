function out=proj_D(Del,D_norms,pert_size,pert_norm,n,m)

[Del_Z,Del_X,Del_U]=sep_ZXU(Del,n,m);
Z_nor=D_norms(1); X_nor=D_norms(2); U_nor=D_norms(3);
Del_Z=proj_block(Del_Z,Z_nor,pert_size,pert_norm);
Del_X=proj_block(Del_X,X_nor,pert_size,pert_norm);
Del_U=proj_block(Del_U,U_nor,pert_size,pert_norm);
out=[Del_Z;Del_X;Del_U];

end

function pD=proj_block(Del_V,V_nor,pert_size,pert_norm)
tmp=pert_size*V_nor;
if strcmp(pert_norm,'Fro')
    if norm(Del_V,'fro')>=tmp
        pD=tmp*Del_V/norm(Del_V,'fro');
    else
        pD=Del_V;
    end
elseif strcmp(pert_norm,'max')
    pD=Del_V;
    pD(pD>tmp)=tmp;
    pD(pD<-tmp)=-tmp;
end
end