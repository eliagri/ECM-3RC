%% estimate_A123_3RC_ECM.m
% 3RC ECM parameter estimation for A123 LiFePO4 cell at 25 degC.
%
% Dataset expected in same folder as this script:
%   OCV:     A1-007-OCV-25-20120905.xlsx
%   Dynamic: A1-007-DST-US06-FUDS-25-20120827.xlsx
%
% Cell info used here:
%   Capacity rating: 1.1 Ah
%   Chemistry: LiFePO4
%   Format: cylindrical, 25.4 mm x 65 mm
%
% Current sign convention in these files:
%   positive current = charge
%   negative current = discharge
%
% The script estimates a GLOBAL 3RC model at 25 degC from the dynamic profile.
% For true SOC-dependent 3RC parameters, HPPC pulses at several SOC points are
% preferable. This dynamic profile can still identify a useful global ECM.

clear; clc; close all;

%% ---------------- User settings ----------------
Q_Ah_nominal = 1.100;                  % Ah
ocvFile = 'A1-007-OCV-25-20120905.xlsx';
dynFile = 'A1-007-DST-US06-FUDS-25-20120827.xlsx';
ocvSheet = 'Channel_1-006';            % best sheet: full Arbin data with Step_Index
 dynSheet = 'Sheet1';                  % compact dynamic sheet with temperature

SOC_LUT = linspace(0,1,101);           % OCV lookup table resolution
fitMaxPoints = 7000;                   % downsample target for fitting speed
useOnlyDynamicAfterTopRest = true;      % true: skip initial conditioning/charge steps

% Initial guesses for A123 1.1 Ah LFP cell, in Ohm and seconds
p0.R0   = 0.040;
p0.R    = [0.010 0.010 0.010];
p0.tau  = [2 20 200];

%% ---------------- Load OCV data ----------------
Tocv = readtable(ocvFile, 'Sheet', ocvSheet, 'VariableNamingRule','preserve');
colO.time = findColumn(Tocv, {'Test_Time(s)','Step_Time(s)','Time'});
colO.step = findColumn(Tocv, {'Step_Index','Step'});
colO.I    = findColumn(Tocv, {'Current(A)','Current'});
colO.V    = findColumn(Tocv, {'Voltage(V)','Voltage'});

tO = Tocv.(colO.time);
stepO = Tocv.(colO.step);
IO = Tocv.(colO.I);
VO = Tocv.(colO.V);

% Find the long low-current discharge segment, normally Step_Index 5.
ocvSegments = segmentByStep(tO, stepO, IO, VO);
longLowDischarge = ocvSegments.meanI < -0.01 & ocvSegments.duration_s > 10000;
if ~any(longLowDischarge)
    error('Could not find a long low-current discharge segment in the OCV file.');
end
[~,loc] = max(ocvSegments.duration_s .* longLowDischarge);
ocvStep = ocvSegments.step(loc);
idxOCV = stepO == ocvStep;

t_ocv = tO(idxOCV);
I_ocv = IO(idxOCV);
V_ocv = VO(idxOCV);

% Capacity from selected C/20 discharge. This is close to 1.06 Ah in this file.
Ah_dis = cumtrapz(t_ocv, max(-I_ocv,0))/3600;
Q_Ah_measured = Ah_dis(end);
if Q_Ah_measured < 0.5 || Q_Ah_measured > 1.5
    warning('Measured OCV discharge capacity looks unusual: %.3f Ah. Using nominal %.3f Ah.', Q_Ah_measured, Q_Ah_nominal);
    Q_Ah = Q_Ah_nominal;
else
    Q_Ah = Q_Ah_measured;
end
SOC_ocv = 1 - Ah_dis/Q_Ah;

% Clean, sort, and interpolate OCV(SOC). LFP has flat regions, so use pchip.
valid = isfinite(SOC_ocv) & isfinite(V_ocv) & SOC_ocv >= 0 & SOC_ocv <= 1;
SOC_ocv = SOC_ocv(valid); V_ocv = V_ocv(valid);
[SOC_sort, order] = sort(SOC_ocv);
V_sort = V_ocv(order);
[SOC_unique, ia] = unique(SOC_sort, 'stable');
V_unique = V_sort(ia);
OCV_LUT = interp1(SOC_unique, V_unique, SOC_LUT, 'pchip', 'extrap');
OCVfun = @(soc) interp1(SOC_LUT, OCV_LUT, min(max(soc,0),1), 'pchip', 'extrap');

fprintf('Selected OCV step: %g\n', ocvStep);
fprintf('Measured C/20 discharge capacity used: %.4f Ah\n', Q_Ah);

%% ---------------- Load dynamic data ----------------
Tdyn = readtable(dynFile, 'Sheet', dynSheet, 'VariableNamingRule','preserve');
colD.time = findColumn(Tdyn, {'Test_Time(s)','Time'});
colD.step = findColumn(Tdyn, {'Step_Index','Step'});
colD.I    = findColumn(Tdyn, {'Current(A)','Current'});
colD.V    = findColumn(Tdyn, {'Voltage(V)','Voltage'});

t = Tdyn.(colD.time);
I = Tdyn.(colD.I);
V = Tdyn.(colD.V);
stepD = Tdyn.(colD.step);

% Remove NaNs and enforce column vectors.
valid = isfinite(t) & isfinite(I) & isfinite(V);
t = t(valid); I = I(valid); V = V(valid); stepD = stepD(valid);
t = t(:); I = I(:); V = V(:); stepD = stepD(:);

% Use the repeated drive-cycle section after the top-rest step.
% In this file the dynamic cycles start at Step_Index 8 after initial charge/rest.
if useOnlyDynamicAfterTopRest
    firstDyn = find(stepD == 8, 1, 'first');
    if ~isempty(firstDyn)
        t = t(firstDyn:end); I = I(firstDyn:end); V = V(firstDyn:end); stepD = stepD(firstDyn:end);
    else
        warning('Could not find Step_Index 8. Using entire dynamic file.');
    end
end

% Reset time to zero.
t = t - t(1);

% Estimate SOC by coulomb counting. Dataset starts dynamic section near full charge.
SOC0 = 1.0;
SOC = SOC0 + cumtrapz(t, I)/(Q_Ah*3600);
SOC = min(max(SOC,0),1);

% Keep data where OCV map is valid and voltage is not hard cutoff artifact.
valid = isfinite(SOC) & V > 1.9 & V < 3.75;
t = t(valid); I = I(valid); V = V(valid); SOC = SOC(valid);

% Downsample for fitting speed, but keep full data for final simulation.
ds = max(1, floor(numel(t)/fitMaxPoints));
idxFit = 1:ds:numel(t);
tFit = t(idxFit); IFit = I(idxFit); VFit = V(idxFit); SOCFit = SOC(idxFit);

fprintf('Dynamic samples used for fitting: %d of %d\n', numel(tFit), numel(t));
fprintf('Dynamic SOC range: %.3f to %.3f\n', min(SOC), max(SOC));

%% ---------------- Fit global 3RC parameters ----------------
% Parameter vector is log-transformed to keep all values positive:
% x = log([R0 R1 R2 R3 tau1 tau2 tau3])
x0 = log([p0.R0 p0.R p0.tau]);

modelResidual = @(x) simulateResidual3RC(x, tFit, IFit, SOCFit, VFit, OCVfun);

if exist('lsqnonlin','file') == 2
    lb = log([0.001 0.0001 0.0001 0.0001 0.05 1 10]);
    ub = log([0.300 0.2000 0.2000 0.2000 20   200 2000]);
    opts = optimoptions('lsqnonlin', 'Display','iter', ...
        'MaxIterations', 300, 'FunctionTolerance',1e-9, 'StepTolerance',1e-9);
    xhat = lsqnonlin(modelResidual, x0, lb, ub, opts);
else
    warning('Optimization Toolbox not found. Using fminsearch instead of lsqnonlin.');
    obj = @(x) mean(modelResidual(x).^2, 'omitnan');
    opts = optimset('Display','iter', 'MaxIter', 800, 'MaxFunEvals', 4000, 'TolX',1e-8, 'TolFun',1e-10);
    xhat = fminsearch(obj, x0, opts);
end

params = unpackParams(xhat);
[Vhat, states] = simulate3RC(params, t, I, SOC, OCVfun);
err = V - Vhat;
RMSE_mV = sqrt(mean(err.^2, 'omitnan'))*1000;
MAE_mV  = mean(abs(err), 'omitnan')*1000;

fprintf('\nEstimated global 3RC parameters:\n');
fprintf('R0  = %.6f Ohm\n', params.R0);
fprintf('R1  = %.6f Ohm, tau1 = %.3f s\n', params.R(1), params.tau(1));
fprintf('R2  = %.6f Ohm, tau2 = %.3f s\n', params.R(2), params.tau(2));
fprintf('R3  = %.6f Ohm, tau3 = %.3f s\n', params.R(3), params.tau(3));
fprintf('RMSE = %.2f mV, MAE = %.2f mV\n', RMSE_mV, MAE_mV);

%% ---------------- Create Simulink-style LUT variables ----------------
% This global fit repeats the same values over SOC. This allows direct use in
% ECM blocks that expect vectors, while making clear that the dynamic profile
% did not robustly identify SOC-dependent parameters.
R0 = params.R0 * ones(size(SOC_LUT));
R1 = params.R(1) * ones(size(SOC_LUT));
R2 = params.R(2) * ones(size(SOC_LUT));
R3 = params.R(3) * ones(size(SOC_LUT));
tau1 = params.tau(1) * ones(size(SOC_LUT));
tau2 = params.tau(2) * ones(size(SOC_LUT));
tau3 = params.tau(3) * ones(size(SOC_LUT));
Em = OCV_LUT;
Capacity_Ah = Q_Ah;
Temperature_C = 25;

save('A123_3RC_ECM_results_25degC.mat', ...
    'Em','SOC_LUT','R0','R1','R2','R3','tau1','tau2','tau3', ...
    'Capacity_Ah','Temperature_C','params','RMSE_mV','MAE_mV');

%% ---------------- Plots ----------------
figure('Name','OCV curve');
plot(SOC_ocv, V_ocv, '.', 'DisplayName','C/20 discharge data'); hold on;
plot(SOC_LUT, OCV_LUT, 'LineWidth',1.5, 'DisplayName','OCV LUT');
grid on; xlabel('SOC'); ylabel('Voltage / OCV [V]'); title('A123 OCV(SOC) from low-current discharge'); legend('Location','best');

figure('Name','Dynamic fit');
plot(t/3600, V, 'DisplayName','Measured voltage'); hold on;
plot(t/3600, Vhat, 'LineWidth',1.2, 'DisplayName','3RC model voltage');
grid on; xlabel('Time [h]'); ylabel('Voltage [V]'); title(sprintf('3RC ECM fit, RMSE %.2f mV', RMSE_mV)); legend('Location','best');

figure('Name','Voltage error');
plot(t/3600, err*1000);
grid on; xlabel('Time [h]'); ylabel('Error [mV]'); title('Measured - model voltage error');

figure('Name','Current and SOC');
yyaxis left; plot(t/3600, I); ylabel('Current [A]');
yyaxis right; plot(t/3600, SOC); ylabel('SOC');
grid on; xlabel('Time [h]'); title('Dynamic profile and estimated SOC');

fprintf('\nSaved: A123_3RC_ECM_results_25degC.mat\n');

%% ---------------- Local functions ----------------
function name = findColumn(T, candidates)
    names = T.Properties.VariableNames;
    lowerNames = lower(names);
    for c = 1:numel(candidates)
        cand = lower(candidates{c});
        match = strcmpi(names, candidates{c}) | contains(lowerNames, cand);
        if any(match)
            name = names{find(match,1,'first')};
            return;
        end
    end
    error('Could not find column. Tried: %s', strjoin(candidates, ', '));
end

function S = segmentByStep(t, step, I, V)
    u = unique(step, 'stable');
    S.step = zeros(numel(u),1);
    S.duration_s = zeros(numel(u),1);
    S.meanI = zeros(numel(u),1);
    S.vStart = zeros(numel(u),1);
    S.vEnd = zeros(numel(u),1);
    for k = 1:numel(u)
        idx = step == u(k);
        S.step(k) = u(k);
        tk = t(idx); Ik = I(idx); Vk = V(idx);
        S.duration_s(k) = max(tk) - min(tk);
        S.meanI(k) = mean(Ik, 'omitnan');
        S.vStart(k) = Vk(1);
        S.vEnd(k) = Vk(end);
    end
    S = struct2table(S);
end

function params = unpackParams(x)
    p = exp(x(:)).';
    params.R0 = p(1);
    Rraw = p(2:4);
    tauraw = p(5:7);
    % Sort by time constant for interpretability: fast, medium, slow.
    [tauSorted, ord] = sort(tauraw);
    params.tau = tauSorted;
    params.R = Rraw(ord);
end

function r = simulateResidual3RC(x, t, I, SOC, V, OCVfun)
    params = unpackParams(x);
    Vhat = simulate3RC(params, t, I, SOC, OCVfun);
    r = Vhat - V;
    r = r(isfinite(r));
end

function [Vhat, xrc] = simulate3RC(params, t, I, SOC, OCVfun)
    n = numel(t);
    xrc = zeros(n,3);
    Vhat = zeros(n,1);
    Vhat(1) = OCVfun(SOC(1)) + I(1)*params.R0 + sum(xrc(1,:));
    for k = 2:n
        dt = max(t(k)-t(k-1), 0);
        for j = 1:3
            a = exp(-dt/params.tau(j));
            xrc(k,j) = a*xrc(k-1,j) + params.R(j)*(1-a)*I(k-1);
        end
        Vhat(k) = OCVfun(SOC(k)) + I(k)*params.R0 + sum(xrc(k,:));
    end
end
