function [ETs, Trs, Ess, Eis, Esbs, SM, RF, GW, snp] = Cal_SiTHM(Rni, Tai, Precii, Pai, Gi, LAIii, Top, s_VODi, ...
        soilpar, pftpar, wa, zgw, snp, optpara, spinfg)

    if spinfg == 1

        ETs = zeros(size(Rni, 1), 1);
        Trs = zeros(size(Rni, 1), 1);
        Ess = zeros(size(Rni, 1), 1);
        Eis = zeros(size(Rni, 1), 1);
        Esbs = zeros(size(Rni, 1), 1);
        SM = zeros(size(Rni, 1), 3);
        RF = zeros(size(Rni, 1), 1);
        GW = zeros(size(Rni, 1), 1);

        for k = 1 : 500 % set the spin-up time (500 years)
            for i = 1:size(Rni, 1)

                Rn = Rni(i, 1);
                Ta = Tai(i, 1);
                Pe = Precii(i, 1);
                Pa = Pai(i, 1);
                G = Gi(i, 1);
                LAI = LAIii(i, 1);
                s_VOD = s_VODi(i, 1);

                [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, ~, ~, ~] = SiTHM(Rn, Ta, Top, ...
                    Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp, optpara);

                ETs(i, 1) = Et;
                Trs(i, 1) = Tr;
                Ess(i, 1) = Es;
                Eis(i, 1) = Ei;
                Esbs(i, 1) = Esb;
                SM(i, :) = wa;
                RF(i, 1) = srf;
                GW(i, 1) = zgw;

            end
        end

    else

        ETs = zeros(size(Rni, 1), 1);
        Trs = zeros(size(Rni, 1), 1);
        Ess = zeros(size(Rni, 1), 1);
        Eis = zeros(size(Rni, 1), 1);
        Esbs = zeros(size(Rni, 1), 1);
        SM = zeros(size(Rni, 1), 3);
        RF = zeros(size(Rni, 1), 1);
        GW = zeros(size(Rni, 1), 1);

        for i = 1:size(Rni, 1)

            Rn = Rni(i, 1);
            Ta = Tai(i, 1);
            Pe = Precii(i, 1);
            Pa = Pai(i, 1);
            G = Gi(i, 1);
            LAI = LAIii(i, 1);
            s_VOD = s_VODi(i, 1);

            [Et, Tr, Es, Ei, Esb, wa, srf, zgw, snp, ~, ~, ~] = SiTHM(Rn, Ta, Top, ...
                Pe, Pa, s_VOD, G, LAI, soilpar, pftpar, wa, zgw, snp, optpara);

            ETs(i, 1) = Et;
            Trs(i, 1) = Tr;
            Ess(i, 1) = Es;
            Eis(i, 1) = Ei;
            Esbs(i, 1) = Esb;
            SM(i, :) = wa;
            RF(i, 1) = srf;
            GW(i, 1) = zgw;

        end

    end
end
