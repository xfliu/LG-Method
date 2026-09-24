vfem_root='/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root,'VFEM_Lib')));
global INTERVAL_MODE; INTERVAL_MODE=0;
p=6; lhat=200; n_eig=4;
NodeList=[0,0,0;1,0,0;0.5,sqrt(3)/2,0;0.5,sqrt(3)/6,sqrt(6)/3];
EL=sort([1,2,3,4],2); FacetList=get_FacetList(EL);
mesh.ElementList=EL; mesh.NodeList=NodeList; mesh.FacetList=FacetList;
mesh.NumF=size(FacetList,1); mesh.NumNode=4; mesh.NumElt=1;
lf=[2 3 4;1 3 4;1 2 4;1 2 3]; E2F=zeros(1,4); F2E=zeros(mesh.NumF,2);
for k=1:4, f=sort(EL(lf(k,:)),2); idx=find(ismember(FacetList,f,'rows'),1); E2F(k)=idx; if F2E(idx,1)==0,F2E(idx,1)=1;else,F2E(idx,2)=1;end; end
mesh.Facet2Element=F2E; mesh.Element2Facet=E2F;
mesh.EdgeList=unique(sort([EL([1 2]);EL([1 3]);EL([1 4]);EL([2 3]);EL([2 4]);EL([3 4])],2),'rows');
mesh.NumEdge=size(mesh.EdgeList,1);

% bubble trial functions
DegP=get_DOF(3,p); DegQ=get_DOF(3,p-4); ijp=get_IJKL(p); ijq=get_IJKL(p-4);
C=zeros(DegP,DegQ); for j=1:DegQ, a=ijq(j,:)+1; r=find(ismember(ijp,a,'rows'),1); C(r,j)=1; end
LN=NodeList; vol=get_volume(LN);
Amm=double(getInnerProdMatrix_Reference(p,p)); Amm1=double(getInnerProdMatrix_Reference(p-1,p-1));
MG=double(get_GradMat(p,LN)); Kp=zeros(DegP); for d=1:3, Kp=Kp+MG(:,:,d)'*Amm1*MG(:,:,d); end
Kp=vol*Kp; Mp=vol*Amm; Kb=C'*Kp*C; Kb=(Kb+Kb')/2; Mb=C'*Mp*C; Mb=(Mb+Mb')/2;
sc=1./sqrt(diag(Mb)); Sd=diag(sc); Ks=Sd*Kb*Sd; Ks=(Ks+Ks')/2; Ms=Sd*Mb*Sd; Ms=(Ms+Ms')/2;
[Y,D]=eig(Ks,Ms); [~,si]=sort(real(diag(D))); Z=Sd*Y(:,si(1:n_eig));
for k=1:n_eig, Z(:,k)=Z(:,k)/sqrt(Z(:,k)'*Mb*Z(:,k)); end
Vdg=C*Z;

rt_d=build_scalar_rt_matrices(mesh,p); A_rt=rt_d.A_rt; B_rt=rt_d.B_rt; M_dg=rt_d.M_dg;
DimRT=rt_d.DimRT; DimDG=rt_d.DimDG;
SP=[A_rt,B_rt';B_rt,-lhat*M_dg]; [Ls,Us,Ps,Qs]=lu(SP);
w1=zeros(DimRT,n_eig);
for i=1:n_eig
  rhs=[zeros(DimRT,1);-M_dg*Vdg(:,i)]; sol=Qs*(Us\(Ls\(Ps*rhs)));
  sol=sol+Qs*(Us\(Ls\(Ps*(rhs-SP*sol)))); w1(:,i)=sol(1:DimRT);
end
fprintf('||w1|| per col: '); fprintf('%.3e ', sqrt(sum(w1.^2,1))); fprintf('\n');
fprintf('w1''*A_rt*w1 diag (double): '); fprintf('%.6e ', diag(w1'*A_rt*w1)); fprintf('\n');

rt_i=build_interval_rt_matrices(mesh,p);
INTERVAL_MODE=1;
Bi=rt_i.B_rt; Ai=rt_i.A_rt;
fprintf('B_rt: max|mid-dbl|=%.3e  maxrad=%.3e\n', full(max(abs(mid(Bi(:))-B_rt(:)))), full(max(rad(Bi(:)))));
fprintf('A_rt: max|mid-dbl|=%.3e  maxrad=%.3e\n', full(max(abs(mid(Ai(:))-A_rt(:)))), full(max(rad(Ai(:)))));
w1i=intval(w1);
q = w1i'*(Ai*w1i);
fprintf('w1''*A_rt_iv*w1 diag: mid / rad:\n');
for k=1:n_eig, fprintf('  k=%d mid=%.6e rad=%.6e\n', k, mid(q(k,k)), rad(q(k,k))); end
INTERVAL_MODE=0;
fprintf('DIAGP6 DONE\n');
