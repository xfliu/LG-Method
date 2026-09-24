vfem_root='/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(vfem_root); addpath(genpath(fullfile(vfem_root,'VFEM_Lib')));
global INTERVAL_MODE; INTERVAL_MODE=0;
p=12;
NodeList=[0,0,0;1,0,0;0.5,sqrt(3)/2,0;0.5,sqrt(3)/6,sqrt(6)/3];
EL=sort([1,2,3,4],2); FacetList=get_FacetList(EL);
mesh.ElementList=EL; mesh.NodeList=NodeList; mesh.FacetList=FacetList;
mesh.NumF=size(FacetList,1); mesh.NumNode=4; mesh.NumElt=1;
lf=[2 3 4;1 3 4;1 2 4;1 2 3]; E2F=zeros(1,4); F2E=zeros(mesh.NumF,2);
for k=1:4, f=sort(EL(lf(k,:)),2); idx=find(ismember(FacetList,f,'rows'),1); E2F(k)=idx; if F2E(idx,1)==0,F2E(idx,1)=1;else,F2E(idx,2)=1;end; end
mesh.Facet2Element=F2E; mesh.Element2Facet=E2F;
mesh.EdgeList=unique(sort([EL([1 2]);EL([1 3]);EL([1 4]);EL([2 3]);EL([2 4]);EL([3 4])],2),'rows');
mesh.NumEdge=size(mesh.EdgeList,1);

rt_d = build_scalar_rt_matrices(mesh, p);     % double, proven
rt_i = build_interval_rt_matrices(mesh, p);   % interval, new

INTERVAL_MODE=1;
Bd = rt_d.B_rt; Bi = rt_i.B_rt;
Ad = rt_d.A_rt; Ai = rt_i.A_rt;
Md = rt_d.M_dg; Mi = rt_i.M_dg;
fprintf('B_rt: max|mid(Bi)-Bd| = %.3e   max rad = %.3e   (max|Bd|=%.3e)\n', ...
    full(max(abs(mid(Bi(:))-Bd(:)))), full(max(rad(Bi(:)))), full(max(abs(Bd(:)))));
fprintf('A_rt: max|mid(Ai)-Ad| = %.3e   max rad = %.3e   (max|Ad|=%.3e)\n', ...
    full(max(abs(mid(Ai(:))-Ad(:)))), full(max(rad(Ai(:)))), full(max(abs(Ad(:)))));
fprintf('M_dg: max|mid(Mi)-Md| = %.3e   max rad = %.3e\n', ...
    full(max(abs(mid(Mi(:))-Md(:)))), full(max(rad(Mi(:)))));
fprintf('nnz: Bd=%d Bi=%d  Ad=%d Ai=%d\n', nnz(Bd), nnz(mid(Bi)), nnz(Ad), nnz(mid(Ai)));
INTERVAL_MODE=0;
fprintf('DIAGBRT DONE\n');
