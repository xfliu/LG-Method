function res = lg_certify(A0_iv, A1_iv, A2_iv, idx, rho, lhat)
%LG_CERTIFY  Certified Lehmann-Goerisch lower bounds for one eigenvalue group.
%
% Given the interval Goerisch matrices A0, A1, A2 (assembled rigorously) and a
% group of trial indices idx that bound the eigenvalues lambda_{m-q+1..m} using
% a rigorous lower bound rho <= lambda_{m+1}, returns certified lower bounds.
%
% Implements Theorem (Lehmann-Goerisch) with the shift technique:
%   A := A0 - rho_sh*A1,   B := A0 - 2*rho_sh*A1 + rho_sh^2*A2,   rho_sh = rho+lhat
%   nu_1<=...<=nu_q eigenvalues of A z = nu B z (B verified SPD)
%   lambda_{m+1-k} >= rho - rho/(1-nu_k),  here returned ascending within idx:
%   lb(j) corresponds to nu_{q+1-j}  (j=1 smallest .. q largest in the group)
% All steps are interval; the guaranteed lower bound is inf(.) of each enclosure.
%
% INPUT
%   A0_iv,A1_iv,A2_iv : (N x N) interval Goerisch matrices
%   idx               : indices of this group (subset of 1:N)
%   rho               : rigorous lower bound of lambda_{m+1} (double)
%   lhat              : shift (double)
% OUTPUT (struct)
%   res.lb            : (q x 1) certified lower bounds (ascending within group)
%   res.nu            : verified nu enclosures
%   res.B_spd         : true if B certified positive definite
%   res.Lambda_max    : enclosure of the largest eig of A0 x = Lambda A1 x (group)
%   res.rho_sh        : rho + lhat
%   res.B_pd_margin   : inf(rho_sh) - sup(Lambda_max)  (>0 => B PD by the Lemma)

    lhat_iv = intval(lhat);
    rho_sh  = intval(rho) + lhat_iv;

    A0s = A0_iv(idx, idx);
    A1s = A1_iv(idx, idx);
    A2s = A2_iv(idx, idx);
    A0s = hull(A0s, A0s'); A1s = hull(A1s, A1s'); A2s = hull(A2s, A2s');
    q = numel(idx);

    A_lg = A0s - rho_sh*A1s;
    B_lg = A0s - 2*rho_sh*A1s + rho_sh^2*A2s;
    A_lg = hull(A_lg, A_lg');
    B_lg = hull(B_lg, B_lg');

    % --- certify B positive definite (verified) ---
    B_spd = false;
    try
        B_spd = isspd(B_lg);   % INTLAB verified SPD test
    catch
        B_spd = false;
    end

    % --- Lemma route: rho_sh > Lambda_n (max eig of A0 x = Lambda A1 x) => B PD ---
    Lambda_max = intval(NaN);
    B_pd_margin = NaN;
    try
        [lam_enc, ~] = veig(A0s, A1s, 1:q);
        Lambda_max = lam_enc(q);
        B_pd_margin = inf(rho_sh) - sup(Lambda_max);
    catch
    end

    % --- nu and lower bounds ---
    res.error = '';
    if q == 1
        nu = A_lg / B_lg;
        lb = inf(rho_sh - rho_sh/(1 - nu) - lhat_iv);
        res.nu = nu;
    else
        try
            [nu, ~] = veig(A_lg, B_lg, 1:q);   % ascending enclosures
            lb = nan(q, 1);
            for j = 1:q
                nuk = nu(q + 1 - j);
                lb(j) = inf(rho_sh - rho_sh/(1 - nuk) - lhat_iv);
            end
            res.nu = nu;
        catch ME
            % veig could not certify B PD / the enclosure is too wide.
            res.error = ME.message;
            res.nu = [];
            lb = nan(q, 1);
        end
    end

    res.lb = lb;
    res.B_spd = B_spd;
    res.Lambda_max = Lambda_max;
    res.rho_sh = rho_sh;
    res.B_pd_margin = B_pd_margin;
    res.B_lg = B_lg;
    res.A_lg = A_lg;
end
