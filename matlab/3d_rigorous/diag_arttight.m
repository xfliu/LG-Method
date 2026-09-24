function diag_arttight(p)
% Test: tight (w1,w1)_{A_rt} via b = MFB\w1 (verified solve) + raw Bernstein mass LMA,
% instead of the ||w1||^2-amplified w1'*A_rt_iv*w1.  One regular tetrahedron.
vfem_root='/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root,'VFEM_Lib')));
global INTERVAL_MODE; INTERVAL_MODE=0;
lhat=200; n_eig=4;
NodeList=[0,0,0;1,0,0;0.5,sqrt(3)/2,0;0.5,sqrt(3)/6,sqrt(6)/3];
EL=sort([1,2,3,4],2); FacetList=get_FacetList(EL);
mesh.ElementList=EL; mesh.NodeList=NodeList; mesh.FacetList=FacetList;
mesh.NumF=size(FacetList,1); mesh.NumNode=4; mesh.NumElt=1;
lf=[2 3 4;1 3 4;1 2 4;1 2 3]; E2F=zeros(1,4); F2E=zeros(mesh.NumF,2);
for k=1:4, f=sort(EL(lf(k,:)),2); idx=find(ismember(FacetList,f,'rows'),1); E2F(k)=idx; if F2E(idx,1)==0,F2E(idx,1)=1;else,F2E(idx,2)=1;end; end
mesh.Facet2Element=F2E; mesh.Element2Facet=E2F;
mesh.EdgeList=unique(sort([EL([1 2]);EL([1 3]);EL([1 4]);EL([2 3]);EL([2 4]);EL([3 4])],2),'rows');
mesh.NumEdge=size(mesh.EdgeList,1);

% bubble trial + SP solve for w1 (double)
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
rt=build_scalar_rt_matrices(mesh,p); A_rt=rt.A_rt; B_rt=rt.B_rt; M_dg=rt.M_dg;
DimRT=rt.DimRT;
SP=[A_rt,B_rt';B_rt,-lhat*M_dg]; [Ls,Us,Ps,Qs]=lu(SP);
w1=zeros(DimRT,n_eig);
for i=1:n_eig
  rhs=[zeros(DimRT,1);-M_dg*Vdg(:,i)]; sol=Qs*(Us\(Ls\(Ps*rhs)));
  sol=sol+Qs*(Us\(Ls\(Ps*(rhs-SP*sol)))); w1(:,i)=sol(1:DimRT);
end
fprintf('p=%d ||w1||=%.2e  double diag(w1''A_rt w1): ', p, max(sqrt(sum(w1.^2,1))));
fprintf('%.6e ', diag(w1'*A_rt*w1)); fprintf('\n');

% --- local RT DOF map and signs (one element) ---
[~,E2Fl,sgn]=mesh_get_Facet2Element_with_sign_fast(mesh);
M=p; DegF=get_DOF(2,M); DegK=get_DOF(3,M); DegRTInner=get_DOF(3,M-1)*3; DegRTElt=DegF*4+DegRTInner;
L2GMapping=zeros(DegRTElt,1); negmask=zeros(1,DegF*4);
for k=1:4
  IdxF=E2Fl(1,k); LocalDOF=(1:DegF)+(k-1)*DegF;
  L2GMapping(LocalDOF)=(1:DegF)+(IdxF-1)*DegF;
  if sgn(1,k)<0, negmask(LocalDOF)=1; end
end
negp=find(negmask>0);
DOF_Facet=mesh.NumF*DegF;
L2GMapping((1:DegRTInner)+DegF*4)=(1:DegRTInner)+DOF_Facet;

% raw Bernstein RT mass LMA (interval) and LMA (double) for check
INTERVAL_MODE=1;
LNi=intval(LN); Kvol=get_volume(LNi);
AMMr=getInnerProdMatrix_Reference(M,M); AMMp=getInnerProdMatrix_Reference(M,M+1); APP=getInnerProdMatrix_Reference(M+1,M+1);
ijl=get_IJKL(M); hidx=find(ijl(:,4)==0);
A_M_M=AMMr*Kvol; AMMplus=AMMp*Kvol; APPv=APP*Kvol;
Up=get_MatDegreeUpByX(M,LNi); UpH=Up(:,hidx,:);
LMA=intval(zeros(DegRTElt)); ia=1:DegK; ib=DegK+(1:DegK); ic=2*DegK+(1:DegK); id=(3*DegK+1):DegRTElt;
LMA(ia,ia)=A_M_M; LMA(ib,ib)=A_M_M; LMA(ic,ic)=A_M_M;
LMd=intval(zeros(numel(id))); for k=1:3, LMd=LMd+UpH(:,:,k)'*APPv*UpH(:,:,k); end
LMA(id,id)=LMd; LMA(ia,id)=AMMplus*UpH(:,:,1); LMA(ib,id)=AMMplus*UpH(:,:,2); LMA(ic,id)=AMMplus*UpH(:,:,3);
LMA(id,ia)=LMA(ia,id)'; LMA(id,ib)=LMA(ib,id)'; LMA(id,ic)=LMA(ic,id)';

MFB = get_rt_funcmat(LN, Kvol, M);   % interval, scaled functional matrix

% per trial: w1_signed -> b = MFB \ w1_signed (verified) -> A2rt = b'LMA b
fprintf('   tight (w1,w1)_A via b=MFB\\w1:  value / rad / (rad/val)\n');
for i=1:n_eig
  w1l = w1(L2GMapping, i); w1l(negp) = -w1l(negp);
  b = verifylss(MFB, intval(w1l));
  q = b' * (LMA * b);
  fprintf('     i=%d  ||b||=%.3e  val=%.6e  rad=%.3e  rad/val=%.2e\n', ...
      i, norm(mid(b)), mid(q), rad(q), rad(q)/abs(mid(q)));
end
INTERVAL_MODE=0;
fprintf('DIAGART DONE p=%d\n', p);
end
