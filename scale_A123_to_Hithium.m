%% scale_A123_to_Hithium_1175Ah_heat_4h.m
% First-order heat-generation estimate for scaling an A123 1.1 Ah LFP ECM
% to a Hithium 1175 Ah LFP prismatic cell under 4-hour constant discharge.
%
% IMPORTANT ASSUMPTION:
% This is an exploratory scaling exercise. The electrical parameters are
% scaled from a small A123 LFP cylindrical cell to a much larger Hithium LFP
% prismatic cell using capacity ratio. This does not replace HPPC testing or
% Hithium-specific parameter identification.
%
% Heat model:
%   Irreversible ECM heat is calculated as electrical dissipation:
%       P_R0  = I^2 * R0
%       P_Ri  = Vrc_i^2 / Ri
%       P_tot = P_R0 + P_R1 + P_R2 + P_R3
%
% ECM convention:
%   Positive current = discharge
%   Vt = OCV(SOC) - I*R0 - Vrc1 - Vrc2 - Vrc3
%
% Outputs:
%   - Figures in hithium_scaled_outputs/
%   - CSV summary table
%   - MAT file with simulation results

clear; clc; close all;

%% ------------------------------------------------------------------------
%  User settings
% -------------------------------------------------------------------------
outDir = "hithium_scaled_outputs";
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% A123 fitted global 3RC parameters at 25 degC
A123.Q_Ah = 1.100;          % nominal capacity used for scaling [Ah]
A123.R0   = 0.158261;       % Ohm
A123.R1   = 0.003202;       % Ohm
A123.tau1 = 3.648;          % s
A123.R2   = 0.020532;       % Ohm
A123.tau2 = 27.800;         % s
A123.R3   = 0.055353;       % Ohm
A123.tau3 = 2000.000;       % s

% Target Hithium cell assumptions
H.Q_Ah      = 1175;         % Hithium nominal capacity [Ah]
H.V_min     = 2.50;         % assumed voltage lower limit [V]
H.V_max     = 3.65;         % assumed voltage upper limit [V]
H.length_m  = 0.580;        % cell dimensions from project note [m]
H.width_m   = 0.075;        % [m]
H.height_m  = 0.216;        % [m]

% Constant discharge case
C_rate_duration_h = 4;      % 4-hour discharge
I_dis_A = H.Q_Ah / C_rate_duration_h;  % positive discharge current [A]
simTime_h = C_rate_duration_h;
dt = 1.0;                  % simulation time step [s]
SOC0 = 1.0;

% OCV model option:
%   "generic_lfp" = built-in smooth LFP-like OCV curve scaled to Hithium voltage range
%   "from_A123_excel" = use A123 low-current OCV file if available in same folder
ocvMode = "from_A123_excel";
ocvFile = "A1-008-OCV-25-20120905.xlsx";
ocvSheet = "Channel_1-005";

% Resistance scaling option:
%   Assumes resistance scales inversely with capacity for cells of same chemistry.
resistanceScale = A123.Q_Ah / H.Q_Ah;

% Time constants:
%   Time constants are retained from the A123 global fit. This implies that
%   the effective capacitances scale inversely with resistance. This is a
%   modeling assumption for exploratory system-level analysis.
keepTauFromA123 = true;

% Optional lumped thermal estimate settings.
% These are uncertain and should be adjusted if actual Hithium mass/cp/h are known.
runLumpedThermalEstimate = true;
T_amb_C = 25;
cellMass_kg = 38;          % placeholder estimate, edit if datasheet mass is known
cp_J_kgK = 1000;           % typical order-of-magnitude for Li-ion cell heat capacity
h_W_m2K = 5;               % natural/weak convection placeholder
useCooling = true;         % true: include h*A*(T-Tamb), false: adiabatic

%% ------------------------------------------------------------------------
%  Scale ECM parameters to Hithium cell
% -------------------------------------------------------------------------
H.R0 = A123.R0 * resistanceScale;
H.R1 = A123.R1 * resistanceScale;
H.R2 = A123.R2 * resistanceScale;
H.R3 = A123.R3 * resistanceScale;

if keepTauFromA123
    H.tau1 = A123.tau1;
    H.tau2 = A123.tau2;
    H.tau3 = A123.tau3;
else
    error('Only keepTauFromA123=true is implemented in this script.');
end

H.C1 = H.tau1 / H.R1;
H.C2 = H.tau2 / H.R2;
H.C3 = H.tau3 / H.R3;

fprintf('--- Scaled Hithium 1175 Ah ECM parameters ---\n');
fprintf('Capacity scaling factor Q_A123/Q_Hithium = %.6g\n', resistanceScale);
fprintf('R0 = %.6g Ohm (%.4f mOhm)\n', H.R0, H.R0*1000);
fprintf('R1 = %.6g Ohm (%.4f mOhm), tau1 = %.3f s\n', H.R1, H.R1*1000, H.tau1);
fprintf('R2 = %.6g Ohm (%.4f mOhm), tau2 = %.3f s\n', H.R2, H.R2*1000, H.tau2);
fprintf('R3 = %.6g Ohm (%.4f mOhm), tau3 = %.3f s\n', H.R3, H.R3*1000, H.tau3);
fprintf('Discharge current for %.1f h discharge: %.2f A\n', C_rate_duration_h, I_dis_A);

%% ------------------------------------------------------------------------
%  Build OCV(SOC)
% -------------------------------------------------------------------------
socGrid = linspace(0, 1, 500);

switch ocvMode
    case "generic_lfp"
        % Smooth LFP-like OCV curve: flat plateau with steeper knees near low/high SOC.
        % Scaled to the Hithium voltage window. This is an approximation only.
        OCV_grid = genericLfpOcv(socGrid, H.V_min, H.V_max);
        ocvFun = @(soc) interp1(socGrid, OCV_grid, min(max(soc,0),1), 'pchip');
        ocvSource = "Generic LFP OCV approximation scaled to Hithium voltage range";

    case "from_A123_excel"
        if ~isfile(ocvFile)
            error('OCV file %s not found. Use ocvMode="generic_lfp" or place file in script folder.', ocvFile);
        end
        [socA, ocvA, Q_ocv] = buildOcvFromA123Excel(ocvFile, ocvSheet);
        % Rescale A123 OCV to Hithium voltage window while preserving OCV shape.
        ocvA_scaled = H.V_min + (ocvA - min(ocvA)) ./ (max(ocvA)-min(ocvA)) .* (H.V_max-H.V_min);
        [socA, ia] = unique(socA, 'stable');
        ocvA_scaled = ocvA_scaled(ia);
        ocvFun = @(soc) interp1(socA, ocvA_scaled, min(max(soc,0),1), 'pchip', 'extrap');
        OCV_grid = ocvFun(socGrid);
        ocvSource = sprintf('A123 OCV file shape rescaled to Hithium voltage window, Q_OCV=%.4f Ah', Q_ocv);

    otherwise
        error('Unknown ocvMode: %s', ocvMode);
end

fprintf('OCV source: %s\n', ocvSource);

%% ------------------------------------------------------------------------
%  Simulate 4-hour constant discharge
% -------------------------------------------------------------------------
t = (0:dt:simTime_h*3600)';
N = numel(t);
I = I_dis_A * ones(N,1);       % positive discharge current

SOC = zeros(N,1);
SOC(1) = SOC0;
Vrc1 = zeros(N,1);
Vrc2 = zeros(N,1);
Vrc3 = zeros(N,1);
Vt = zeros(N,1);
OCV = zeros(N,1);

P_R0 = zeros(N,1);
P_R1 = zeros(N,1);
P_R2 = zeros(N,1);
P_R3 = zeros(N,1);
P_total = zeros(N,1);

for k = 1:N
    OCV(k) = ocvFun(SOC(k));
    Vt(k) = OCV(k) - I(k)*H.R0 - Vrc1(k) - Vrc2(k) - Vrc3(k);

    % Heat generation from ECM dissipative elements.
    P_R0(k) = I(k)^2 * H.R0;
    P_R1(k) = Vrc1(k)^2 / H.R1;
    P_R2(k) = Vrc2(k)^2 / H.R2;
    P_R3(k) = Vrc3(k)^2 / H.R3;
    P_total(k) = P_R0(k) + P_R1(k) + P_R2(k) + P_R3(k);

    if k < N
        dtk = t(k+1) - t(k);

        % Exact discrete update for constant current over each interval.
        a1 = exp(-dtk/H.tau1);
        a2 = exp(-dtk/H.tau2);
        a3 = exp(-dtk/H.tau3);
        Vrc1(k+1) = a1*Vrc1(k) + (1-a1)*H.R1*I(k);
        Vrc2(k+1) = a2*Vrc2(k) + (1-a2)*H.R2*I(k);
        Vrc3(k+1) = a3*Vrc3(k) + (1-a3)*H.R3*I(k);

        SOC(k+1) = SOC(k) - I(k)*dtk/(H.Q_Ah*3600);
        SOC(k+1) = min(max(SOC(k+1),0),1);
    end
end

E_heat_J = trapz(t, P_total);
E_heat_Wh = E_heat_J / 3600;
E_out_Wh = trapz(t, max(Vt,0).*I) / 3600;
P_avg_W = E_heat_J / (t(end)-t(1));
P_end_W = P_total(end);

fprintf('\n--- 4-hour discharge heat estimate ---\n');
fprintf('Average heat generation: %.2f W per cell\n', P_avg_W);
fprintf('End-of-discharge heat generation: %.2f W per cell\n', P_end_W);
fprintf('Total heat over %.1f h: %.2f Wh per cell\n', simTime_h, E_heat_Wh);
fprintf('Electrical output energy estimate: %.2f kWh per cell\n', E_out_Wh/1000);
fprintf('Heat/output-energy ratio: %.2f %%\n', 100*E_heat_Wh/E_out_Wh);
fprintf('Terminal voltage range: %.3f to %.3f V\n', min(Vt), max(Vt));

%% ------------------------------------------------------------------------
%  Optional round-trip-efficiency cross-check
% -------------------------------------------------------------------------
eta_rt = 0.95;
loss_rt_Wh = E_out_Wh*(1/eta_rt - 1);
P_loss_rt_avg_W = loss_rt_Wh / simTime_h;

fprintf('\n--- 95%% round-trip efficiency reference ---\n');
fprintf('If 95%% round-trip efficiency is used as a coarse reference,\n');
fprintf('equivalent total loss over %.1f h discharge energy is %.2f Wh, or %.2f W average.\n', ...
    simTime_h, loss_rt_Wh, P_loss_rt_avg_W);
fprintf('This is not used in the ECM heat calculation; it is only a sanity check.\n');

%% ------------------------------------------------------------------------
%  Optional lumped thermal estimate
% -------------------------------------------------------------------------
if runLumpedThermalEstimate
    A_surface = 2*(H.length_m*H.width_m + H.length_m*H.height_m + H.width_m*H.height_m);
    Tcell = zeros(N,1);
    Tcell(1) = T_amb_C;
    P_cooling = zeros(N,1);

    for k = 1:N-1
        dtk = t(k+1)-t(k);
        if useCooling
            P_cooling(k) = h_W_m2K*A_surface*(Tcell(k)-T_amb_C);
        else
            P_cooling(k) = 0;
        end
        dTdt = (P_total(k) - P_cooling(k))/(cellMass_kg*cp_J_kgK);
        Tcell(k+1) = Tcell(k) + dTdt*dtk;
    end
    if useCooling
        P_cooling(end) = h_W_m2K*A_surface*(Tcell(end)-T_amb_C);
    end

    fprintf('\n--- Lumped thermal estimate ---\n');
    fprintf('Mass = %.1f kg, cp = %.0f J/kgK, h = %.1f W/m2K, area = %.3f m2\n', ...
        cellMass_kg, cp_J_kgK, h_W_m2K, A_surface);
    if useCooling
        fprintf('Estimated cell temperature after %.1f h with simple convection: %.2f degC\n', simTime_h, Tcell(end));
    else
        fprintf('Estimated adiabatic cell temperature after %.1f h: %.2f degC\n', simTime_h, Tcell(end));
    end
else
    Tcell = [];
    P_cooling = [];
    A_surface = NaN;
end

%% ------------------------------------------------------------------------
%  Figures
% -------------------------------------------------------------------------
% 1) Current and SOC
fig = figure('Name','Current and SOC','Color','w');
yyaxis left;
plot(t/3600, I, 'LineWidth', 1.4); ylabel('Discharge current [A]');
yyaxis right;
plot(t/3600, SOC, 'LineWidth', 1.4); ylabel('SOC [-]');
grid on; xlabel('Time [h]'); title('4-hour constant-current discharge: current and SOC');
saveas(fig, fullfile(outDir, '01_current_SOC.png'));

% 2) OCV and terminal voltage
fig = figure('Name','Voltage response','Color','w');
plot(t/3600, OCV, 'LineWidth', 1.3, 'DisplayName','OCV(SOC)'); hold on;
plot(t/3600, Vt, 'LineWidth', 1.3, 'DisplayName','Terminal voltage');
grid on; xlabel('Time [h]'); ylabel('Voltage [V]'); title('Predicted voltage during 4-hour discharge');
legend('Location','best');
saveas(fig, fullfile(outDir, '02_voltage_response.png'));

% 3) OCV curve used
fig = figure('Name','OCV curve used','Color','w');
plot(socGrid, OCV_grid, 'LineWidth', 1.5);
grid on; xlabel('SOC [-]'); ylabel('OCV [V]'); title('OCV(SOC) used for scaled Hithium simulation');
set(gca,'XDir','reverse');
saveas(fig, fullfile(outDir, '03_OCV_curve_used.png'));

% 4) Scaled resistance comparison
fig = figure('Name','Resistance scaling','Color','w');
barVals_mOhm = [A123.R0 A123.R1 A123.R2 A123.R3; H.R0 H.R1 H.R2 H.R3]*1000;
bar(barVals_mOhm');
grid on; ylabel('Resistance [mOhm]'); xticklabels({'R0','R1','R2','R3'});
legend({'A123 original','Hithium scaled'}, 'Location','best');
title('Resistance scaling by capacity ratio');
saveas(fig, fullfile(outDir, '04_resistance_scaling.png'));

% 5) RC branch voltages
fig = figure('Name','RC branch voltages','Color','w');
plot(t/3600, Vrc1, 'LineWidth', 1.2, 'DisplayName','V_{RC1}'); hold on;
plot(t/3600, Vrc2, 'LineWidth', 1.2, 'DisplayName','V_{RC2}');
plot(t/3600, Vrc3, 'LineWidth', 1.2, 'DisplayName','V_{RC3}');
grid on; xlabel('Time [h]'); ylabel('RC polarization voltage [V]');
title('Predicted RC polarization voltages'); legend('Location','best');
saveas(fig, fullfile(outDir, '05_RC_branch_voltages.png'));

% 6) Heat generation components
fig = figure('Name','Heat generation components','Color','w');
plot(t/3600, P_R0, 'LineWidth', 1.2, 'DisplayName','R0 heat'); hold on;
plot(t/3600, P_R1, 'LineWidth', 1.2, 'DisplayName','R1 heat');
plot(t/3600, P_R2, 'LineWidth', 1.2, 'DisplayName','R2 heat');
plot(t/3600, P_R3, 'LineWidth', 1.2, 'DisplayName','R3 heat');
plot(t/3600, P_total, 'k--', 'LineWidth', 1.5, 'DisplayName','Total heat');
grid on; xlabel('Time [h]'); ylabel('Heat generation [W]');
title('ECM irreversible heat generation by component'); legend('Location','best');
saveas(fig, fullfile(outDir, '06_heat_generation_components.png'));

% 7) Cumulative heat energy
fig = figure('Name','Cumulative heat energy','Color','w');
E_heat_cum_Wh = cumtrapz(t, P_total)/3600;
plot(t/3600, E_heat_cum_Wh, 'LineWidth', 1.5);
grid on; xlabel('Time [h]'); ylabel('Cumulative heat [Wh]');
title('Cumulative heat generated during 4-hour discharge');
saveas(fig, fullfile(outDir, '07_cumulative_heat.png'));

% 8) Heat as fraction of output power
fig = figure('Name','Power comparison','Color','w');
P_out = Vt.*I;
plot(t/3600, P_out/1000, 'LineWidth', 1.3, 'DisplayName','Electrical output power [kW]'); hold on;
plot(t/3600, P_total, 'LineWidth', 1.3, 'DisplayName','Heat generation [W]');
grid on; xlabel('Time [h]'); title('Electrical output power and heat generation');
legend('Location','best');
saveas(fig, fullfile(outDir, '08_output_power_and_heat.png'));

% 9) Optional temperature estimate
if runLumpedThermalEstimate
    fig = figure('Name','Lumped thermal estimate','Color','w');
    plot(t/3600, Tcell, 'LineWidth', 1.5); hold on;
    yline(T_amb_C, '--', 'Ambient');
    grid on; xlabel('Time [h]'); ylabel('Cell temperature [degC]');
    title('Simple lumped thermal estimate');
    saveas(fig, fullfile(outDir, '09_lumped_temperature_estimate.png'));
end

%% ------------------------------------------------------------------------
%  Save summary
% -------------------------------------------------------------------------
summary = table;
summary.Parameter = [
    "A123_capacity_Ah";
    "Hithium_capacity_Ah";
    "Capacity_scaling_QA123_over_QH";
    "Discharge_current_A";
    "H_R0_Ohm";
    "H_R1_Ohm";
    "H_R2_Ohm";
    "H_R3_Ohm";
    "tau1_s";
    "tau2_s";
    "tau3_s";
    "Average_heat_W";
    "End_heat_W";
    "Total_heat_Wh";
    "Output_energy_Wh";
    "Heat_output_ratio_percent";
    "RT_eff_reference_loss_Wh";
    "RT_eff_reference_avg_loss_W";
    "Final_temp_C_lumped"
    ];

finalTemp = NaN;
if runLumpedThermalEstimate
    finalTemp = Tcell(end);
end

summary.Value = [
    A123.Q_Ah;
    H.Q_Ah;
    resistanceScale;
    I_dis_A;
    H.R0;
    H.R1;
    H.R2;
    H.R3;
    H.tau1;
    H.tau2;
    H.tau3;
    P_avg_W;
    P_end_W;
    E_heat_Wh;
    E_out_Wh;
    100*E_heat_Wh/E_out_Wh;
    loss_rt_Wh;
    P_loss_rt_avg_W;
    finalTemp
    ];

writetable(summary, fullfile(outDir, 'hithium_scaled_heat_summary.csv'));
save(fullfile(outDir, 'hithium_scaled_heat_results.mat'), ...
    'A123','H','t','I','SOC','OCV','Vt','Vrc1','Vrc2','Vrc3', ...
    'P_R0','P_R1','P_R2','P_R3','P_total','E_heat_Wh','E_out_Wh', ...
    'Tcell','P_cooling','summary','ocvSource');

fprintf('\nOutputs saved in folder: %s\n', outDir);

%% ========================================================================
%  Local functions
% ========================================================================
function OCV = genericLfpOcv(soc, Vmin, Vmax)
    % Smooth generic LFP OCV shape with flat mid-SOC plateau.
    soc = min(max(soc,0),1);
    z = soc;

    % Unscaled LFP-like curve: low/high knees plus mild plateau slope.
    Vraw = 3.22 ...
        + 0.10*z ...
        + 0.20./(1 + exp(-(z-0.92)/0.035)) ...
        - 0.42./(1 + exp((z-0.06)/0.025));

    % Normalize to requested voltage window.
    Vraw0 = min(Vraw);
    Vraw1 = max(Vraw);
    OCV = Vmin + (Vraw - Vraw0)./(Vraw1 - Vraw0).*(Vmax - Vmin);
end

function [SOC_unique, V_unique, Q_Ah] = buildOcvFromA123Excel(filename, sheetName)
    T = readtable(filename, 'Sheet', sheetName, 'VariableNamingRule','preserve');
    timeCol = findColumn(T, {'Test_Time(s)','Step_Time(s)','Time'});
    stepCol = findColumn(T, {'Step_Index','Step'});
    currCol = findColumn(T, {'Current(A)','Current'});
    voltCol = findColumn(T, {'Voltage(V)','Voltage'});

    t = T.(timeCol);
    step = T.(stepCol);
    I = T.(currCol);
    V = T.(voltCol);

    steps = unique(step(isfinite(step)));
    bestStep = NaN;
    bestDur = -Inf;
    for s = steps(:)'
        idx = step == s & isfinite(t) & isfinite(I) & isfinite(V);
        if nnz(idx) < 10, continue; end
        dur = max(t(idx)) - min(t(idx));
        meanI = mean(I(idx),'omitnan');
        if meanI < -0.01 && abs(meanI) < 1.0 && dur > bestDur
            bestStep = s;
            bestDur = dur;
        end
    end
    if isnan(bestStep)
        error('Could not find a low-current discharge step in OCV file.');
    end

    idx = step == bestStep & isfinite(t) & isfinite(I) & isfinite(V);
    t = t(idx); I = I(idx); V = V(idx);
    [t, order] = sort(t);
    I = I(order); V = V(order);

    Q_Ah = trapz(t, abs(I))/3600;
    Ah_cum = cumtrapz(t, abs(I))/3600;
    SOC = 1 - Ah_cum/Q_Ah;
    SOC = min(max(SOC,0),1);

    [SOC_sort, order] = sort(SOC, 'ascend');
    V_sort = V(order);
    [SOC_unique, ia] = unique(SOC_sort, 'stable');
    V_unique = V_sort(ia);
end

function colName = findColumn(T, candidates)
    names = string(T.Properties.VariableNames);
    normNames = lower(regexprep(names, '[^a-zA-Z0-9]', ''));
    for c = 1:numel(candidates)
        cand = lower(regexprep(string(candidates{c}), '[^a-zA-Z0-9]', ''));
        hit = find(contains(normNames, cand), 1, 'first');
        if ~isempty(hit)
            colName = T.Properties.VariableNames{hit};
            return;
        end
    end
    error('Could not find column. Candidates: %s. Available columns: %s', ...
        strjoin(string(candidates), ', '), strjoin(names, ', '));
end
