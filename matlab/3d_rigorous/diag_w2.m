% Isolate the w2 verified solve on the REAL rhs (one-element regular tet, p=12).
vfem_root = '/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root,'VFEM_Lib')));
global INTERVAL_MODE; INTERVAL_MODE=0;
p=12; lhat=200; n_eig=4;
NodeList=[0,0,0;1,0,0;0.5,sqrt(3)/2,0;0.5,sqrt(3)/6,sqrt(6)/3];
EL=sort([1,2,3,4],2); FacetList=get_FacetList(EL);
mesh.ElementList=EL; mesh.NodeList=NodeList; mesh.FacetList=FacetList;
mesh.NumF=size(FacetList,1); mesh.NumNode=4; mesh.NumElt=1;
lf=[2 3 4;1 3 4;1 2 4;1 2 3]; E2F=zeros(1,4); F2E=zeros(mesh.NumF,2);
for k=1:4, f=sort(EL(lf(k,:)),2); idx=find(ismember(FacetList,f,'rows'),1); E2F(k)=idx; if F2E(idx,1)==0,F2E(idx,1)=1;else,F2E(idx,2)=1;end; end
mesh.Facet2Element=F2E; mesh.Element2Facet=E2F;
mesh.EdgeList=unique(sort([EL([1 2]);EL([1 3]);EL([1 4]);EL([2 3]);EL([2 4]);EL([3 4])],2),'rows');
mesh.NumEdge=size(mesh.EdgeList,1);

% bubble trial functions (as in driver)
DegP=get_DOF(3,p); DegQ=get_DOF(3,p-4); ijp=get_IJKL(p); ijq=get_IJKL(p-4);
C=zeros(DegP,DegQ);
for j=1:DegQ, a=ijq(j,:)+1; r=find(ismember(ijp,a,'rows'),1); C(r,j)=1; end
LN=NodeList; vol=get_volume(LN);
Amm=double(getInnerProdMatrix_Reference(p,p)); Amm1=double(getInnerProdMatrix_Reference(p-1,p-1));
MG=double(get_GradMat(p,LN)); Kp=zeros(DegP);
for d=1:3, Kp=Kp+MG(:,:,d)'*Amm1*MG(:,:,d); end
Kp=vol*Kp; Mp=vol*Amm;
Kb=C'*Kp*C; Kb=(Kb+Kb')/2; Mb=C'*Mp*C; Mb=(Mb+Mb')/2;
sc=1./sqrt(diag(Mb)); S=diag(sc); Ks=S*Kb*S; Ks=(Ks+Ks')/2; Ms=S*Mb*S; Ms=(Ms+Ms')/2;
[Y,D]=eig(Ks,Ms); [~,si]=sort(real(diag(D))); Z=S*Y(:,si(1:n_eig));
for k=1:n_eig, Z(:,k)=Z(:,k)/sqrt(Z(:,k)'*Mb*Z(:,k)); end
Vdg=C*Z;

rt=build_scalar_rt_matrices(mesh,p); A_rt=rt.A_rt; B_rt=rt.B_rt; M_dg=rt.M_dg;
DimRT=rt.DimRT; DimDG=rt.DimDG;
SP=[A_rt,B_rt';B_rt,-lhat*M_dg]; [Ls,Us,Ps,Qs]=lu(SP);
w1=zeros(DimRT,n_eig);
for i=1:n_eig
  rhs=[zeros(DimRT,1);-M_dg*Vdg(:,i)]; sol=Qs*(Us\(Ls\(Ps*rhs)));
  sol=sol+Qs*(Us\(Ls\(Ps*(rhs-SP*sol)))); w1(:,i)=sol(1:DimRT);
end

% Build the interval DG block and real rhs (element 1), test verifylss strategies.
INTERVAL_MODE=1;
Kvol=get_volume(intval(LN));
A_M_M = getInnerProdMatrix(p,p,Kvol);          % tight DG mass block
B_rt_iv=intval(B_rt); M_dg_iv=A_M_M;            % one element: M_dg block == A_M_M
w1_iv=intval(w1); Vdg_iv=intval(Vdg);
rhs = B_rt_iv*w1_iv + M_dg_iv*Vdg_iv;           % (DegP x n)
DV = lhat*A_M_M;
fprintf('||mid(rhs)||=%.3e  cond(mid DV)~%.3e\n', norm(mid(rhs),'fro'), cond(mid(DV)));
test = @(tag,x) fprintf('%-10s rad(w2)=%.3e  residual rad=%.3e\n', tag, ...
    full(max(rad(x(:)))), full(max(rad(reshape(DV*x-rhs,[],1)))));
try, x=verifylss(full(DV),rhs); test('plain',x); catch e, fprintf('plain threw %s\n',e.message); end
try, x=verifylss(full(DV),rhs,'normal'); test('normal',x); catch e, fprintf('normal threw %s\n',e.message); end
try, x=verifylss(full(DV),rhs,'illco'); test('illco',x); catch e, fprintf('illco threw %s\n',e.message); end
INTERVAL_MODE=0;
fprintf('DIAGW2 DONE\n');
