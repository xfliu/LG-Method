function diag_degree(p)
% Double-precision LG cluster bound on the 8-element regular tet at degree p.
% A2 midpoint is identical to the old thin-wrap code (b'LMA b == w1'A_rt w1 in
% double), so this isolates the DEGREE effect on the cluster bound from rigor.
vfem_root='/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root,'VFEM_Lib')));
addpath(fullfile(fileparts(vfem_root),'vfem2d','veigs'));
global INTERVAL_MODE; INTERVAL_MODE=0;
n_eig=4; lhat=200; rho2=264.5223894731; rho5=386.4419682205;
NodeList=[0,0,0;1,0,0;0.5,sqrt(3)/2,0;0.5,sqrt(3)/6,sqrt(6)/3];
mesh0=make_one_tet_mesh(NodeList,[1,2,3,4]); mesh=refine_tet_mesh(mesh0);
NumElt=mesh.NumElt; NL=mesh.NodeList; EL=mesh.ElementList;
DegK=get_DOF(3,p);
[L2G,DOF_BD,DimCG]=local_register_cg_dof(mesh,p);
interior=setdiff(1:DimCG,DOF_BD);
Amm=double(getInnerProdMatrix_Reference(p,p)); Amm1=double(getInnerProdMatrix_Reference(p-1,p-1));
iK=zeros(NumElt*DegK^2,1);jK=iK;vK=iK;vM=iK;ptr=0;
for e=1:NumElt
  LN=NL(EL(e,:),:); vol=get_volume(LN); dofs=L2G(e,:);
  MG=get_GradMat(p,LN); LK=zeros(DegK);
  for d=1:3, LK=LK+double(MG(:,:,d))'*Amm1*double(MG(:,:,d)); end
  LK=LK*vol; LM=Amm*vol; [jl,il]=meshgrid(1:DegK,1:DegK); rr=ptr+(1:DegK^2);
  iK(rr)=dofs(il(:)); jK(rr)=dofs(jl(:)); vK(rr)=LK(:); vM(rr)=LM(:); ptr=ptr+DegK^2;
end
K=sparse(iK,jK,vK,DimCG,DimCG); Mm=sparse(iK,jK,vM,DimCG,DimCG);
K_int=(K(interior,interior)+K(interior,interior)')/2; M_int=(Mm(interior,interior)+Mm(interior,interior)')/2;
[Ve,De]=eigs(K_int,M_int,n_eig+2,1.0,struct('disp',0));
[lr,si]=sort(real(diag(De))); Ve=real(Ve(:,si(1:n_eig))); lr=lr(1:n_eig);
v_all=zeros(DimCG,n_eig);
for k=1:n_eig, vf=zeros(DimCG,1); vf(interior)=Ve(:,k); vf=vf/sqrt(vf(interior)'*M_int*vf(interior)); v_all(:,k)=vf; end
rt=build_scalar_rt_matrices(mesh,p); A_rt=rt.A_rt;B_rt=rt.B_rt;M_dg=rt.M_dg;
DimRT=rt.DimRT;DimDG=rt.DimDG;DegKd=rt.DegK;
SP=[A_rt,B_rt';B_rt,-lhat*M_dg]; [Ls,Us,Ps,Qs]=lu(SP);
w1=zeros(DimRT,n_eig); w2=zeros(DimDG,n_eig); Vdg=zeros(DimDG,n_eig);
for i=1:n_eig
  for e=1:NumElt, Vdg((e-1)*DegKd+(1:DegKd),i)=v_all(L2G(e,:),i); end
  rhs=[zeros(DimRT,1);-M_dg*Vdg(:,i)]; sol=Qs*(Us\(Ls\(Ps*rhs)));
  sol=sol+Qs*(Us\(Ls\(Ps*(rhs-SP*sol)))); w1(:,i)=sol(1:DimRT); w2(:,i)=sol(DimRT+1:end);
end
A0=v_all(interior,:)'*((K_int+lhat*M_int)*v_all(interior,:));
A1=v_all(interior,:)'*(M_int*v_all(interior,:));
A2=zeros(n_eig);
for i=1:n_eig, for j=i:n_eig
  A2(i,j)=w1(:,i)'*A_rt*w1(:,j)+lhat*(w2(:,i)'*M_dg*w2(:,j)); A2(j,i)=A2(i,j);
end, end
A0=(A0+A0')/2; A1=(A1+A1')/2; A2=(A2+A2')/2;
fprintf('p=%d 8elt: ub=%.6f %.6f %.6f %.6f\n', p, lr);
fprintf('  A2(k,k)*(lr+lhat): '); for k=1:n_eig, fprintf('%.8f ', A2(k,k)*(lr(k)+lhat)); end; fprintf('\n');
% lambda1
r2=rho2+lhat; a=A0(1,1)-r2*A1(1,1); b=A0(1,1)-2*r2*A1(1,1)+r2^2*A2(1,1); nu=a/b;
lb1=r2-r2/(1-nu)-lhat;
% cluster
idx=2:4; r5=rho5+lhat; Acl=A0(idx,idx)-r5*A1(idx,idx); Bcl=A0(idx,idx)-2*r5*A1(idx,idx)+r5^2*A2(idx,idx);
nus=sort(real(eig(Acl,Bcl)));
fprintf('  cluster nu: %.6f %.6f %.6f\n', nus);
lbc=zeros(3,1); for j=1:3, nuk=nus(3+1-j); lbc(j)=r5-r5/(1-nuk)-lhat; end
fprintf('  lb: lam1=%.6f (gap %.4f%%) | cluster=%.6f %.6f %.6f (gap %.4f%%)\n', ...
  lb1,(lr(1)-lb1)/lr(1)*100, lbc, (lr(2)-min(lbc))/lr(2)*100);
fprintf('DIAGDEG DONE p=%d\n', p);
end

function mesh=make_one_tet_mesh(NodeList,ElementList)
ElementList=sort(ElementList,2); FacetList=get_FacetList(ElementList);
mesh.ElementList=ElementList; mesh.NodeList=NodeList; mesh.FacetList=FacetList;
mesh.NumF=size(FacetList,1); mesh.NumNode=size(NodeList,1); mesh.NumElt=size(ElementList,1);
edges=[ElementList(:,[1 2]);ElementList(:,[1 3]);ElementList(:,[1 4]);ElementList(:,[2 3]);ElementList(:,[2 4]);ElementList(:,[3 4])];
mesh.EdgeList=unique(sort(edges,2),'rows'); mesh.NumEdge=size(mesh.EdgeList,1);
lf=[2 3 4;1 3 4;1 2 4;1 2 3]; NumElt=size(ElementList,1); NumF=size(FacetList,1);
E2F=zeros(NumElt,4); F2E=zeros(NumF,2);
for e=1:NumElt, for k=1:4, f=sort(ElementList(e,lf(k,:)),2); ix=find(ismember(FacetList,f,'rows'),1); E2F(e,k)=ix; if F2E(ix,1)==0,F2E(ix,1)=e;else,F2E(ix,2)=e;end; end, end
mesh.Facet2Element=F2E; mesh.Element2Facet=E2F;
end
