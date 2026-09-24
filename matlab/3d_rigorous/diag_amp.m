function diag_amp()
% Why is the interval cluster bound loose while the double one is tight?
% Load the saved p10 8-element interval matrices and inspect the cancellation
% B = A0 - 2 r5 A1 + r5^2 A2 and the veig nu enclosure widths.
vfem_root='/home/xfliu/Workspace/New_Project_SchrodingerEVP/VFEM_LIB/VFEM3D';
old=pwd; cd(vfem_root); my_intlab_mode_config; cd(old);
addpath(fullfile(fileparts(vfem_root),'vfem2d','veigs'));
S=load('regular_tet_ref1_p10_lg_rigorous_results.mat'); R=S.R;
A0=R.A0_iv; A1=R.A1_iv; A2=R.A2_iv; lhat=R.lhat; rho5=R.rho5;
fprintf('max rad: A0=%.2e A1=%.2e A2=%.2e\n', max(rad(A0(:))),max(rad(A1(:))),max(rad(A2(:))));
idx=2:4; r5=intval(rho5+lhat);
A=A0(idx,idx)-r5*A1(idx,idx);
B=A0(idx,idx)-2*r5*A1(idx,idx)+r5^2*A2(idx,idx);
A=hull(A,A'); B=hull(B,B');
fprintf('r5=%.2f  r5^2=%.3e\n', mid(r5), mid(r5)^2);
fprintf('B midpoint diag: %.4f %.4f %.4f\n', mid(diag(B)));
A2c = A2(idx,idx);
fprintf('B max rad: %.3e   (r5^2 * rad(A2) = %.3e)\n', max(rad(B(:))), mid(r5)^2*max(rad(A2c(:))));
[nu,~]=veig(A,B,1:3);
fprintf('veig nu enclosures:\n');
for j=1:3, fprintf('  nu%d in [%.5f, %.5f]  width %.3e\n', j, inf(nu(j)), sup(nu(j)), sup(nu(j))-inf(nu(j))); end
for j=1:3
  nuk=nu(3+1-j); lb=inf(r5 - r5/(1-nuk) - intval(lhat));
  fprintf('  cluster lb(%d) = %.6f\n', j, lb);
end
fprintf('DIAGAMP DONE\n');
end
