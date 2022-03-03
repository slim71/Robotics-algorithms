function qdot = rp(tt, qq, sym_Js, sym_twistes, lambda)
    % initialization
    qdot_lp1 = zeros(max(size(qq)),1);

    % init
    qdot_k = qdot_lp1;
    JRA_k = [];

    for k = [max(size(sym_twistes)):-1:1]
        Jk = sym_Js{k}; % function handle
        twistk = sym_twistes{k}; % function handle

        Jk_qq = Jk(qq);
        JRA_k = [Jk_qq; JRA_k];

        [JRAPseudo, Tk, Hk] = LSEstimateSequential(JRA_k, size(Jk_qq, 1));

        qdot_k = qdot_k + Tk * pinv(Jk_qq*Tk) * (twistk(tt)' - Jk_qq * qdot_k);
    end

    qdot = zeros(max(size(qq)),1);
    for k = [1:max(size(qq))]
        qdot(k,1) = qdot_k(k);
    end
end