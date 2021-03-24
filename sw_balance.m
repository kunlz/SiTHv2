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

    % water holding capacity (%)
    whc = wfc - wp;  % field capacity (wfc, %) minus wilting point (wp, %)

    d = [d1 d2 d3]
    d_sum = sum(d);
    zm_s = [soildeep(1), soildeep(1) + soildeep(2), sum(soildeep)];

    % Case 1: groundwater table above layer 1
    if d_sum <= zm_s(1)

        % layer #1 new soil water content in unsaturated zone,
        % and possible exceeded water downward to layer #2
        w1_unsat = wo(1) + IWS ./ (d1 .* whc); % soil water content of unsaturated zone
        
        if w1_unsat > 1
            vw_1 = (w1_unsat - 1) .* (d1 .* whc); % exceeded soil water
            w1_unsat = 1;
        else
            % no additional exceeded soil water
            vw_1 = 0;
        end

        % calculate the adjusted swc#1 with considering the groundwater depth
        w1_adj = ((w1_unsat - w1_sat) * d(1) + w1_sat * soildeep(1)) / soildeep(1);

        % calculate the saturated swc #2 #3
        w2 = w2_sat;
        w3 = w3_sat;
        

        % Soil evaporation and Plant transpiration
        % Calculate the moisture constrains for plant and soil in unsaturated zone
        [S_plant1, S_soil1] = swc_stress(w1_unsat, soilpar);

        % vertical root density for layer #1, PFTpar(1)
        s11 = PFTpar(1) .* (w1_adj / w1_sat).^b; 
        s22 = PFTpar(1) .* (w1_adj / w1_sat).^b + PFTpar(1) + (1 - PFTpar(1) - PFTpar(2)); 
        pTr_1 = s11 ./ s22 .* pEc .* (1 - wet); 

        % for unsaturated zone
        Tr_p_us = pTr_1 .* d1 .* w1_unsat ./ (d1 .* w1_unsat + (zm1 - d1) .* w1_sat);

        % Actrual Transpiration and Soil Evaporation from unsaturated zone
        Tr1_us = S_plant1 .* S_tem .* Tr_p_us;
        Es1_us = (wet + (1 - wet) .* S_soil1) .* pEs; % Only consider about the first layer

        % layer #2

        % layer #3

        % Case 2: groundwater table between layer 1 and layer 2
    elseif d_sum > zm_s(1) && d_sum <= zm_s(2)

        % Case 3: groundwater table between layer 2 and layer 3
    elseif d_sum > zm_s(2) && d_sum < zm_s(3)

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

    % Percolation from the first layer
    Dmin = 0.048; % mm/day
    Dmax = 4.8; % mm/day

    if w1_adj > 0.75

        Perc1 = Dmin * w1_adj + (Dmax - Dmin) * w1_adj^1.5;
    else

        Perc1 = Dmin * w1_adj;
    end

    % update new soil water content
    S6 = soilpar(6);
    S9 = soilpar(9);

    % First layer
    wa_1 = w1_unsat - (Es1 + Tr1 + Perc1) ./ (soildeep(1) .* soilpar(3));

    if wa_1 < 0
        wa_1 = 0;
        d1_1 = soildeep(1);
    end

    % redistributed to unsaterated zone with considering groundwater depth
    wa_11 = 1 + (wa_1 - 1) * soildeep(1) / d1_1;

    % theta_1 = soilpar(6)+soilpar(3).*wa_11;

    theta_1(wa_1 < 0) = S6(wa_1 < 0);
    wa_1(wa_1 < 0) = 0;

    % #2 Second soil layer water content
    w_2 = wo(2) + vw_1 ./ (d(2) .* soilpar(4));

    if w_2 > 1
        vw_2 = (w_2 - 1) .* (d(2) .* soilpar(4)); % downward water 2
        w_2 = 1;
    else
        vw_2 = 0;
    end

    % #3 Third soil layer water content
    w_3 = wo(3) + vw_2 ./ (d(3) .* soilpar(5));

    if w_3 > 1
        vw_3 = (w_3 - 1) .* (d(3) .* soilpar(5)); % downward water 3
        w_3 = 1;
    else
        vw_3 = 0;
    end

    deep_drainage = vw_3;

    wa = [w1_unsat, w_2, w_3];

    % ------------------------------------------------------------
    % The Percolation Process -- ref.to :  Ducharne & Polcher 1998
    % ------------------------------------------------------------
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
