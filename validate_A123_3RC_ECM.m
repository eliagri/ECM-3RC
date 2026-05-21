%% validate_A123_3RC_ECM.m
% Validate a fixed 3RC ECM against independent A123 dynamic profiles.
%
% This script does NOT refit the ECM parameters. It uses the supplied
% R0/R1/R2/R3/tau values, builds OCV(SOC) from the low-current OCV test,
% simulates the measured dynamic current profiles, and compares simulated
% terminal voltage against measured terminal voltage.
%
% Data files expected in the same folder as this script:
%   - A1-008-OCV-25-20120905.xlsx
%   - A1-008-DST-US06-FUDS-25-20120827.xlsx
%
% Current convention in the A123 files:
%   negative current = discharge
%   positive current = charge
%
% ECM convention used in the simulation:
%   Idis = -Iraw, so positive Idis = discharge
%   Vt = OCV(SOC) - Idis*R0 - Vrc1 - Vrc2 - Vrc3
%
% Author: generated for A123 validation workflow

clear; clc; close all;

%% ------------------------------------------------------------------------
%  User settings
% -------------------------------------------------------------------------
ocvFile = "A1-008-OCV-25-20120905.xlsx";
dynFile = "A1-008-DST-US06-FUDS-25-20120827.xlsx";

% Fixed ECM parameters from previous fitting. Do not change during validation.
R0   = 0.158261;      % Ohm
R1   = 0.003202;      % Ohm
tau1 = 3.648;         % s
R2   = 0.020532;      % Ohm
tau2 = 27.800;        % s
R3   = 0.055353;      % Ohm
tau3 = 2000.000;      % s

% Dynamic Step_Index values to validate.
% In this file, the long dynamic blocks are normally Step_Index 8, 16 and 24.
% If you only want two profiles, set e.g. validationSteps = [16 24];
validationSteps = [8 16 24];

% Initial SOC estimation mode:
%   "ocv"   = estimate SOC from first measured voltage of each dynamic segment
%   numeric = use the same fixed initial SOC for all segments, e.g. 1.0
initialSOCMode = "ocv";

% Optional: ignore the first seconds of each segment when computing metrics.
% This can be useful if the segment begins with a short rest/transition.
metricStartTime_s = 0;

% Output folder
outDir = "validation_outputs";
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

%% ------------------------------------------------------------------------
%  1) Build OCV(SOC) from the low-current OCV file
% -------------------------------------------------------------------------
Tocv = readBatterySheet(ocvFile, "Channel_1-005");

% Select the low-current discharge step automatically:
% negative current, low current magnitude, and longest duration.
ocvStep = findLowCurrentDischargeStep(Tocv);
idxOcv = Tocv.Step_Index == ocvStep;

t_ocv = Tocv.Test_Time_s(idxOcv);
I_ocv_raw = Tocv.Current_A(idxOcv);
V_ocv = Tocv.Voltage_V(idxOcv);

% Capacity from low-current discharge, using absolute current integration.
Q_Ah = trapz(t_ocv, abs(I_ocv_raw)) / 3600;

% SOC during OCV discharge: starts at 1 and ends near 0.
Ah_cum = cumtrapz(t_ocv, abs(I_ocv_raw)) / 3600;
SOC_ocv = 1 - Ah_cum ./ Q_Ah;
SOC_ocv = min(max(SOC_ocv, 0), 1);

% Sort for interpolation and remove duplicate SOC entries.
[SOC_sort, sortIdx] = sort(SOC_ocv, 'ascend');
V_sort = V_ocv(sortIdx);
[SOC_unique, uniqueIdx] = unique(SOC_sort, 'stable');
V_unique = V_sort(uniqueIdx);

% OCV lookup function, clipped to measured SOC range.
ocvFun = @(soc) interp1(SOC_unique, V_unique, min(max(soc, 0), 1), 'linear', 'extrap');

% Inverse OCV lookup used only for initial SOC estimation.
% Voltage is expected to increase with SOC for the discharge OCV curve.
[V_for_inv, invIdx] = sort(V_unique, 'ascend');
SOC_for_inv = SOC_unique(invIdx);
[V_for_inv_unique, vUniqueIdx] = unique(V_for_inv, 'stable');
SOC_for_inv_unique = SOC_for_inv(vUniqueIdx);
socFromOcv = @(v) min(max(interp1(V_for_inv_unique, SOC_for_inv_unique, v, 'linear', 'extrap'), 0), 1);

fprintf('OCV step selected: Step_Index %d\n', ocvStep);
fprintf('Measured capacity from OCV discharge: %.4f Ah\n', Q_Ah);

% Plot OCV curve.
fig = figure('Name', 'OCV curve', 'Color', 'w');
plot(SOC_ocv, V_ocv, '.', 'DisplayName', 'Measured low-current voltage'); hold on;
socGrid = linspace(0, 1, 400);
plot(socGrid, ocvFun(socGrid), 'LineWidth', 1.5, 'DisplayName', 'Interpolated OCV(SOC)');
grid on; xlabel('SOC'); ylabel('Voltage [V]'); title('OCV(SOC) from low-current discharge');
legend('Location', 'best');
set(gca, 'XDir', 'reverse');
saveas(fig, fullfile(outDir, 'OCV_curve.png'));

%% ------------------------------------------------------------------------
%  2) Load dynamic file and validate each selected profile
% -------------------------------------------------------------------------
Tdyn = readBatterySheet(dynFile, "Channel_1-005");

metrics = table();
allResults = struct();

for k = 1:numel(validationSteps)
    stepId = validationSteps(k);
    idx = Tdyn.Step_Index == stepId;

    if ~any(idx)
        warning('Step_Index %d was not found in the dynamic file. Skipping.', stepId);
        continue;
    end

    t = Tdyn.Test_Time_s(idx);
    Iraw = Tdyn.Current_A(idx);
    Vmeas = Tdyn.Voltage_V(idx);

    % Time relative to profile start.
    t = t - t(1);

    % ECM convention: positive Idis means discharge.
    Idis = -Iraw;

    % Initial SOC.
    if isnumeric(initialSOCMode)
        SOC0 = initialSOCMode;
    else
        % Estimate from the first measured voltage in the segment.
        % The first seconds in these A123 dynamic segments are close to rest,
        % making this a reasonable calibration-free SOC initialization.
        SOC0 = socFromOcv(Vmeas(1));
    end

    sim = simulate3RC(t, Idis, SOC0, Q_Ah, ocvFun, R0, R1, tau1, R2, tau2, R3, tau3);

    Vsim = sim.Vt;
    Verr = Vsim - Vmeas;

    metricIdx = t >= metricStartTime_s;
    rmse = sqrt(mean(Verr(metricIdx).^2));
    voltageRange = max(Vmeas(metricIdx)) - min(Vmeas(metricIdx));
    rmse_percent = 100 * rmse / voltageRange;
    mae = mean(abs(Verr(metricIdx)));
    maxAbsErr = max(abs(Verr(metricIdx)));
    bias = mean(Verr(metricIdx));

    profileName = sprintf('Step_%d', stepId);

    newRow = table(stepId, SOC0, rmse, mae, maxAbsErr, bias, min(sim.SOC), max(sim.SOC), ...
        'VariableNames', {'Step_Index','Initial_SOC','RMSE_V','MAE_V','MaxAbsError_V','Bias_V','Min_SOC','Max_SOC'});
    metrics = [metrics; newRow]; %#ok<AGROW>

    allResults(k).profileName = profileName;
    allResults(k).stepId = stepId;
    allResults(k).t = t;
    allResults(k).Iraw = Iraw;
    allResults(k).Idis = Idis;
    allResults(k).Vmeas = Vmeas;
    allResults(k).Vsim = Vsim;
    allResults(k).Verr = Verr;
    allResults(k).SOC = sim.SOC;
    allResults(k).Vrc = sim.Vrc;

    fprintf('\n%s validation:\n', profileName);
    fprintf('  Initial SOC: %.3f\n', SOC0);
    fprintf('  RMSE: %.4f V (%.2f%%) | MAE: %.4f V | Max abs error: %.4f V | Bias: %.4f V\n', ...
        rmse, rmse_percent, mae, maxAbsErr, bias);

    % Voltage fit plot.
    fig = figure('Name', profileName + " voltage", 'Color', 'w');
    plot(t/60, Vmeas, 'LineWidth', 1.1, 'DisplayName', 'Measured voltage'); hold on;
    plot(t/60, Vsim, '--', 'LineWidth', 1.2, 'DisplayName', 'Simulated 3RC voltage');
    grid on; xlabel('Time [min]'); ylabel('Voltage [V]');
    title(sprintf('%s: measured vs simulated voltage', profileName), 'Interpreter', 'none');
    legend('Location', 'best');
    saveas(fig, fullfile(outDir, profileName + "_voltage_fit.png"));

    % Error plot.
    fig = figure('Name', profileName + " error", 'Color', 'w');
    plot(t/60, Verr * 1000, 'LineWidth', 1.1);
    grid on; xlabel('Time [min]'); ylabel('Voltage error [mV]');
    title(sprintf('%s: voltage error, RMSE = %.1f mV', profileName, rmse*1000), 'Interpreter', 'none');
    yline(0, 'k-');
    saveas(fig, fullfile(outDir, profileName + "_voltage_error.png"));

    % Current and SOC plot.
    fig = figure('Name', profileName + " current soc", 'Color', 'w');
    yyaxis left;
    plot(t/60, Iraw, 'LineWidth', 1.0);
    ylabel('Raw current [A]');
    yyaxis right;
    plot(t/60, sim.SOC, 'LineWidth', 1.2);
    ylabel('SOC [-]');
    grid on; xlabel('Time [min]');
    title(sprintf('%s: current and simulated SOC', profileName), 'Interpreter', 'none');
    saveas(fig, fullfile(outDir, profileName + "_current_SOC.png"));
end

%% ------------------------------------------------------------------------
%  3) Save results
% -------------------------------------------------------------------------
writetable(metrics, fullfile(outDir, 'validation_metrics.csv'));
save(fullfile(outDir, 'validation_results.mat'), 'metrics', 'allResults', ...
    'Q_Ah', 'R0', 'R1', 'tau1', 'R2', 'tau2', 'R3', 'tau3', ...
    'SOC_unique', 'V_unique', 'validationSteps');

fprintf('\nValidation complete. Results saved in folder: %s\n', outDir);
disp(metrics);

%% ========================================================================
%  Local functions
% ========================================================================
function T = readBatterySheet(fileName, sheetName)
    % Read Arbin battery sheet and standardize important variable names.
    opts = detectImportOptions(fileName, 'Sheet', sheetName, 'VariableNamingRule', 'preserve');
    Traw = readtable(fileName, opts);

    T = table();
    T.Test_Time_s = getColumn(Traw, {"Test_Time(s)", "Test_Time_s", "Test Time(s)"});
    T.Step_Time_s = getColumn(Traw, {"Step_Time(s)", "Step_Time_s", "Step Time(s)"});
    T.Step_Index = getColumn(Traw, {"Step_Index", "Step Index"});
    T.Current_A = getColumn(Traw, {"Current(A)", "Current_A", "Current"});
    T.Voltage_V = getColumn(Traw, {"Voltage(V)", "Voltage_V", "Voltage"});

    % Remove rows with missing critical values.
    valid = ~isnan(T.Test_Time_s) & ~isnan(T.Step_Index) & ~isnan(T.Current_A) & ~isnan(T.Voltage_V);
    T = T(valid, :);
end

function x = getColumn(T, candidates)
    names = string(T.Properties.VariableNames);
    x = [];
    for i = 1:numel(candidates)
        idx = find(names == candidates{i}, 1);
        if ~isempty(idx)
            x = T{:, idx};
            return;
        end
    end
    error('Could not find any of these columns: %s', strjoin(candidates, ', '));
end

function stepId = findLowCurrentDischargeStep(T)
    steps = unique(T.Step_Index);
    bestScore = -inf;
    stepId = steps(1);

    for i = 1:numel(steps)
        s = steps(i);
        idx = T.Step_Index == s;
        I = T.Current_A(idx);
        t = T.Test_Time_s(idx);

        meanI = mean(I, 'omitnan');
        meanAbsI = mean(abs(I), 'omitnan');
        duration = max(t) - min(t);

        % Low-current discharge should have negative current, long duration,
        % and small current magnitude. For this A123 OCV file this selects Step 5.
        if meanI < 0 && meanAbsI < 0.2 && duration > 1000
            score = duration / max(meanAbsI, 1e-6);
            if score > bestScore
                bestScore = score;
                stepId = s;
            end
        end
    end
end

function sim = simulate3RC(t, Idis, SOC0, Q_Ah, ocvFun, R0, R1, tau1, R2, tau2, R3, tau3)
    n = numel(t);
    SOC = zeros(n,1);
    Vrc = zeros(n,3);
    Vt = zeros(n,1);

    SOC(1) = min(max(SOC0, 0), 1);
    Vt(1) = ocvFun(SOC(1)) - Idis(1)*R0 - sum(Vrc(1,:));

    R = [R1 R2 R3];
    tau = [tau1 tau2 tau3];

    for j = 2:n
        dt = t(j) - t(j-1);
        if dt <= 0 || isnan(dt)
            dt = median(diff(t), 'omitnan');
        end

        % Coulomb counting. Positive Idis discharges the cell.
        SOC(j) = SOC(j-1) - (Idis(j-1) * dt) / (Q_Ah * 3600);
        SOC(j) = min(max(SOC(j), 0), 1);

        % Exact discrete-time update for each RC branch with zero-order hold current.
        for m = 1:3
            alpha = exp(-dt / tau(m));
            Vrc(j,m) = alpha * Vrc(j-1,m) + R(m) * (1 - alpha) * Idis(j-1);
        end

        Vt(j) = ocvFun(SOC(j)) - Idis(j)*R0 - sum(Vrc(j,:));
    end

    sim.SOC = SOC;
    sim.Vrc = Vrc;
    sim.Vt = Vt;
end
