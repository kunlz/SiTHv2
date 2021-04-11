% ----------------------- %
%    Soil Water Balance   %
% ----------------------- %
function [Tr, Es, wa, theta_new, DeepDrain] = sw_balance(Topt, Ta, IWS, pTr_ly, pEs, wo, wet, soilpar, soildeep, d)
    % ---------- function input -------
    % soilpar   :   Soil parameters
    % wo        :   Previous step soil water content
    % IWS       :   Infiltrated water to soil
    % soildepth :   The depths of different soil layers
    % wet       :   Canopy wetness fraction
    % pTr_ly    :   The potential Transpiration from different soil layers
    % pSoilE    :   The potential soil Evaporation from first soil layer
    % Ta        :   Air temperature
    % ---------- function output -------
    % Tr        :   Actual Plant Transpiration, mm/day
    % Es        :   Actual Soil Evaporation, mm/day
    % wa_new    :   Actual Soil moisture for three soil layers, %
    % theta_new :
    % DeepDrain :
    % ----------
    % Constrains of temperature
    S_tem = StressTemp(Topt, Ta);

    R_sb_max = 39; % mm day-1
    f = 1.25e-3; % mm-1

    % water holding capacity (%)
    whc = wfc - wp; % field capacity (wfc, %) minus wilting point (wp, %)

    d = [d1 d2 d3]
    z_gw = sum(d);
    zm_s = [soildeep(1), soildeep(1) + soildeep(2), sum(soildeep)];

    % Case 1: groundwater table above layer 1
    if z_gw <= zm_s(1)

        % layer #1
        % existed water column in the unsat-zone
        wc_s1 = d1 * wa1_unsat;

        % maximum water column
        wc_m1 = d1 * theta_sat;

        if wc_s1 + IWS >= wc_m1

            wa1 = theta_sat; % current soil water content
            vw1 = wc_s1 + IWS - wc_m1; % exceeded water

        else

            % soil water content in unsaturated zone
            wa1_unsat = wa1_unsat + IWS / d1;

            % calculate the adjusted swc#1 with considering the groundwater depth
            wa1 = (wa1_unsat * d1 + theta_sat * (z1 - d1)) / z1;

            vw1 = 0; % no exceeded water

        end

        % layer #2
        wa2 = theta_sat;

        % layer #3
        wa3 = theta_sat;

        % -------------------------- evapotration section

        % the maximum transpiration rate of each soil layer, Tr_p
        % Get all avaliable water contents through root distribution
        wr = r1 * (wa1 / theta_sat)^b + r2 * (wa2 / theta_sat)^b + r3 * (wa3 / theta_sat)^b;

        % Root distribution adjusted by soil water content
        beta1 = r1 * (wa1 / theta_sat)^b / wr;
        beta2 = r2 * (wa2 / theta_sat)^b / wr;
        beta3 = r3 * (wa3 / theta_sat)^b / wr;

        % potentail transpiration rate for different layers
        Tr_p1 = beta1 * pEc;
        Tr_p2 = beta2 * pEc;
        Tr_p3 = beta3 * pEc;

        % divide Tr_p1 into unsaturated zone and saturated zone
        Tr_p1_u = Tr_p1 * (d1 * wa1_unsat) / (d1 * wa1_unsat + (z1 - d1) * theta_sat);
        Tr_p1_g = Tr_p1 * ((z1 - d1) * theta_sat) / (d1 * wa1_unsat + (z1 - d1) * theta_sat);

        % Calculate the moisture constrains for plant and soil in unsaturated zone
        [S_plant, S_soil] = swc_stress([wa1, wa2, wa3], soilpar);

        % Actrual Transpiration
        Tr1 = (1 - wetness) .* S_plant1 .* S_tem .* Tr_p1;
        Tr2 = (1 - wetness) .* S_plant2 .* S_tem .* Tr_p2;
        Tr3 = (1 - wetness) .* S_plant3 .* S_tem .* Tr_p3;

        Tr1_g = Tr1 * ((z1 - d1) * theta_sat) / (d1 * wa1_unsat + (z1 - d1) * theta_sat);

        % Soil Evaporation
        Es = S_soil1 .* pEs; % Only consider about the first layer

        % ------------------------------------------------------------
        % The Percolation Process -- ref.to :  Ducharne & Polcher 1998
        % ------------------------------------------------------------
        % layer #1
        f1 = min(ks, d1 * (wa1_unsat - theta_fc)); % drainage from unsaturated zone

        if f1 < 0
            f1 = 0;
        end

        F1 = f1 + vw1; % total water recharge to groundwater

        % total transpiration from groundwater
        Tr_g = Tr1_g + Tr2 + Tr3;

        % R_sb groundwater discaharge
        R_sb_max = 39; % mm day-1
        f = 1.25e-3; % mm-1
        R_sb = R_sb_max * exp(-f * z_gw); 
        
        % variation of water storaged in the saturated zone
        delta_w = F1 - Tr_g - R_sb;

        % changes in groundwater table depth
        delata_z_gw = delta_w / (theta_sat - wa1); 

    % Case 2: groundwater table between layer 1 and layer 2
    elseif z_gw > zm_s(1) && z_gw <= zm_s(2)
        
        % layer #1
        % existed water column in the unsat-zone
        wc_s1 = d1 * wa1_unsat;

        % maximum water column
        wc_m1 = d1 * theta_sat;

        if wc_s1 + IWS >= wc_m1

            wa1 = theta_sat; % current soil water content
            vw1 = wc_s1 + IWS - wc_m1; % exceeded water

        else

            % soil water content in unsaturated zone
            wa1 = wa1_unsat + IWS / d1;
            vw1 = 0; % no exceeded water

        end

        % layer #2
        % existed water column in the unsat-zone
        wc_s2 = d2 * wa2_unsat; 

        % maximum water column in d2
        wc_m2 = d2 * theta_sat;

        if wc_s2 + vw1 >= wc_m2
            wa2 = theta_sat;
            vw2 = wc_s2 + vw1 - wc_m2;
        else
            % soil water content in unsaturated zone
            wa2_unsat = wa2_unsat + vw1 / d2;
            % calculate the adjusted swc#1 with considering the groundwater depth
            wa2 = (wa2_unsat * d2 + theta_sat * (z2 - d2)) / z2;

            vw2 = 0; % no exceeded water
        end

        % layer #3
        wa3 = theta_sat;

        % -------------------------- evapotration section
        % distributed the potential T to different layers
        [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, D50, c, b, theta_sat, wet, zm); 

        % divide Tr_p2 into unsaturated zone and saturated zone
        Tr_p2_u = Tr_p2 * (d2 * wa2_unsat) / (d2 * wa2_unsat + (z2 - d2) * theta_sat);
        Tr_p2_g = Tr_p2 * ((z2 - d2) * theta_sat) / (d2 * wa2_unsat + (z2 - d2) * theta_sat);

        % Calculate the moisture constrains for plant and soil in unsaturated zone
        [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar); 
        [f_sm2, ~] = swc_stress(wa2_unsat, soilpar); 

        % Actrual Transpiration
        Tr1 = f_sm1 .* S_tem .* Tr_p1;
        Tr2_u = f_sm2 .* S_tem .* Tr_p2_u;
        Tr2_g = S_tem .* Tr_p2_g;
        Tr2 = Tr2_u + Tr2_g;
        Tr3 = S_tem .* Tr_p3; 
        Tr = Tr1 + Tr2 + Tr3;

        % Soil Evaporation
        Es = f_sm_s1 .* pEs; % Only consider about the first layer

        % ------------------------------
        % The groundwater table changes
        % ------------------------------
        % layer #1
        f1 = min(ks, d1 * (wa1 - theta_fc)); % drainage from unsaturated zone

        if f1 < 0
            f1 = 0;
        end

        F1 = f1 + vw1 + vw2; % total water recharge to groundwater

        % total transpiration from groundwater
        Tr_g = Tr2_g + Tr3;

        % R_sb groundwater discaharge
        R_sb = R_sb_max * exp(-f * z_gw); 
        
        % variation of water storaged in the saturated zone
        delta_w = F1 - Tr_g - R_sb;

        % changes in groundwater table depth
        delata_z_gw = delta_w / (theta_sat - wa2_unsat); 

    % Case 3: groundwater table between layer 2 and layer 3
    elseif z_gw > zm_s(2) && z_gw < zm_s(3)

        % layer #1
        % existed water column in the unsat-zone
        wa1_unsat = wa1; 
        wc_s1 = d1 * wa1_unsat;

        % maximum water column
        wc_m1 = d1 * theta_sat;

        if wc_s1 + IWS >= wc_m1

            wa1 = theta_sat; % current soil water content
            vw1 = wc_s1 + IWS - wc_m1; % exceeded water

        else

            % soil water content in unsaturated zone
            wa1 = wa1_unsat + IWS / d1;
            vw1 = 0; % no exceeded water

        end

        % layer #2
        % existed water column in the unsat-zone
        wa2_unsat = wa2;
        wc_s2 = d2 * wa2_unsat; 

        % maximum water column in d2
        wc_m2 = d2 * theta_sat;

        if wc_s2 + vw1 >= wc_m2
            wa2 = theta_sat;
            vw2 = wc_s2 + vw1 - wc_m2;
        else
            % soil water content in unsaturated zone
            wa2 = wa2_unsat + vw1 / d2; 
            vw2 = 0; % no exceeded water
        end

        % layer #3
        % existed water column in the unsat-zone
        wc_s3 = d3 * wa3_unsat; 

        % maximum water column in d2
        wc_m3 = d3 * theta_sat;

        if wc_s3 + vw2 >= wc_m3
            wa3 = theta_sat;
            vw3 = wc_s3 + vw2 - wc_m3;
        else
            % soil water content in unsaturated zone
            wa3_unsat = wa3_unsat + vw2 / d3; 
            % calculate the adjusted swc#1 with considering the groundwater depth
            wa3 = (wa3_unsat * d3 + theta_sat * (z3 - d3)) / z3;

            vw3 = 0; % no exceeded water
        end

        % -------------------------- evapotration section
        % distributed the potential T to different layers
        [Tr_p1, Tr_p2, Tr_p3] = pTr_partition(pEc, wa1, wa2, wa3, D50, c, b, theta_sat, wet, zm); 

        % divide Tr_p3 into unsaturated zone and saturated zone
        Tr_p3_u = Tr_p3 * (d3 * wa3_unsat) / (d3 * wa3_unsat + (z3 - d3) * theta_sat);
        Tr_p3_g = Tr_p3 * ((z3 - d3) * theta_sat) / (d3 * wa3_unsat + (z3 - d3) * theta_sat);

        % Calculate the moisture constrains for plant and soil in unsaturated zone
        [f_sm1, f_sm_s1] = swc_stress(wa1, soilpar); 
        [f_sm2, ~] = swc_stress(wa2, soilpar); 
        [f_sm3, ~] = swc_stress(wa3_unsat, soilpar); 

        % Actrual Transpiration
        Tr1 = f_sm1 .* S_tem .* Tr_p1;
        Tr2 = f_sm2 .* S_tem .* Tr_p2;
        Tr3_g = S_tem .* Tr_p3_g;
        Tr2 = Tr2_u + Tr2_g;
        Tr3 = S_tem .* Tr_p3; 
        Tr = Tr1 + Tr2 + Tr3;

        % Soil Evaporation
        Es = f_sm_s1 .* pEs; % Only consider about the first layer

        % ------------------------------
        % The groundwater table changes
        % ------------------------------
        % layer #1
        f1 = min(ks, d1 * (wa1 - theta_fc)); % drainage from unsaturated zone

        if f1 < 0
            f1 = 0;
        end

        F1 = f1 + vw1 + vw2; % total water recharge to groundwater

        % total transpiration from groundwater
        Tr_g = Tr2_g + Tr3;

        % R_sb groundwater discaharge
        R_sb = R_sb_max * exp(-f * z_gw); 
        
        % variation of water storaged in the saturated zone
        delta_w = F1 - Tr_g - R_sb;

        % changes in groundwater table depth
        delata_z_gw = delta_w / (theta_sat - wa2_unsat); 


    % Case 4: groundwater table below layer 3
    else

    end

    % unsatturated soil water infiltration for different layers
    % #1 Fist soil layer water content

    % Calculate the moisture constrains for plant and soil (3 layers)
    [S_plant1, S_soil1] = swc_stress(w1_adj, soilpar);

    % % Cal Actrual Transpiration and Soil Evaporation
    Tr1 = S_plant1 .* S_tem .* pTr_ly(1);
    Es1 = (wet + (1 - wet) .* S_soil1) .* pEs; % Only consider about the first layer

    d = 1.5;
    % Percolation from the first layer
    Dmin = 0.048; % mm/day
    Dmax = 4.8; % mm/day
    Perc1 = Dmin .* w1_unsat;
    Perc_x = Dmin .* w1_unsat + (Dmax - Dmin) .* w1_unsat.^d;
    Perc1(w1_unsat >= 0.75) = Perc_x(w1_unsat >= 0.75);

    % Drainage from the secdond layer
    Dmin = 0.0005 .* 24; % mm/day
    Dmax = 0.05 .* 24; % mm/day
    Perc_2 = Dmin .* w_2;
    Perc_x = Dmin .* w_2 + (Dmax - Dmin) .* w_2.^d;
    Perc_2(w_2 >= 0.75) = Perc_x(w_2 >= 0.75);

    % Drainage from the third layer
    Dmin = 0.0005 .* 24; % mm/day
    Dmax = 0.05 .* 24; % mm/day
    Perc_3 = Dmin .* w_3;
    Perc_x = Dmin .* w_3 + (Dmax - Dmin) .* w_3.^d;
    Perc_3(w_3 >= 0.75) = Perc_x(w_3 >= 0.75);

    % Second layer
    wa_2 = w_2 + (Perc1 - Tr_2 - Perc_2) ./ (soildeep(2) .* soilpar(:, :, 4));
    theta_2 = soilpar(:, :, 6) + soilpar(:, :, 4) .* wa_2;
    theta_2(wa_2 < 0) = S6(wa_2 < 0);
    wa_2(wa_2 < 0) = 0;

    VI2 = (wa_2 - 1) .* (soildeep(2) .* soilpar(:, :, 4));
    VI2(wa_2 <= 1) = 0;
    theta_2(wa_2 > 1) = S9(wa_2 > 1);
    wa_2(wa_2 > 1) = 1;

    % Third layer
    wa_3 = w_3 + (Perc_2 + VI2 - Tr_3 - Perc_3) ./ (soildeep(3) .* soilpar(:, :, 5));
    theta_3 = soilpar(:, :, 6) + soilpar(:, :, 5) .* wa_3;
    theta_3(wa_3 < 0) = S6(wa_3 < 0);
    wa_3(wa_3 < 0) = 0;

    VI3 = (wa_3 - 1) .* (soildeep(3) .* soilpar(:, :, 5));
    VI3(wa_3 <= 1) = 0;
    DeepDrain = DeepDrain + VI3;
    theta_3(wa_3 > 1) = S9(wa_3 > 1);
    wa_3(wa_3 > 1) = 1;

    %%

    % Constrains of temperature
    S_tem = StressTemp(Topt, Ta);

    % Calculate the moisture constrains for plant and soil (3 layers)
    [S_plant, S_soil] = StressSoilM(wa, soilpar);

    % % Cal Actrual Transpiration and Soil Evaporation
    Tr = 0;

    for i = 1:3
        Tr = Tr + S_plant(i) .* S_tem .* pTr_ly(i);
    end

    Es = (wet + (1 - wet) .* S_soil(1)) .* pEs; % Only consider about the first layer

    % Output
    wa(:, :, 1) = wa_1;
    wa(:, :, 2) = wa_2;
    wa(:, :, 3) = wa_3;

    theta_new(:, :, 1) = theta_1;
    theta_new(:, :, 2) = theta_2;
    theta_new(:, :, 3) = theta_3;
end
