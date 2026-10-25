clc;
clear;
close all;

%% =========================================================================
% 1. READ RAW DATA FROM EXCEL
% =========================================================================
file_excel = '..\..\File_path\TiO2.txt';

try
    fprintf('Reading XRD data file: %s ...\n', file_excel);
    data_excel = readtable(file_excel);
catch ME
    error('Failed to read Excel file. Ensure the file path is correct and the file exists.\nError Message: %s', ME.message);
end

% Extract Column 1 (2-Theta) and Column 2 (Total Intensity)
twoTheta  = data_excel{:, 1};
Intensity = data_excel{:, 2};

% Validate Data Integrity
if isempty(twoTheta) || isempty(Intensity) || length(twoTheta) ~= length(Intensity)
    error('2-Theta and Intensity data are invalid or have mismatched row lengths.');
end

%% =========================================================================
% 2. LIGHT SMOOTHING OF RAW DATA
% =========================================================================
span     = round(0.01 * length(Intensity)); 
I_smooth = smoothdata(Intensity, 'gaussian', span);

%% =========================================================================
% 3. AMORPHOUS BASELINE (ASYMMETRIC WEIGHTED LEAST SQUARES)
% =========================================================================
N      = length(I_smooth);
p      = 0.002;   % Balanced asymmetry penalty
lambda = 1e5;     % Smoothness parameter to follow the broad amorphous halo

e = ones(N, 1);
D = spdiags([e -2*e e], 0:2, N-2, N);
w = ones(N, 1);

for iter = 1:15
    W = spdiags(w, 0, N, N);
    C = chol(W + lambda * (D' * D));
    I_amorphous = C \ (C' \ (w .* I_smooth));
    w = p * (I_smooth > I_amorphous) + (1 - p) * (I_smooth <= I_amorphous);
end

I_amorphous   = min(I_amorphous, I_smooth);
I_amorphous   = max(0, I_amorphous);
I_crystalline = max(0, I_smooth - I_amorphous);

%% =========================================================================
% 4. CALCULATE PERCENT CRYSTALLINITY INDEX (CI %)
% =========================================================================
area_cryst = trapz(twoTheta, I_crystalline);
area_amorf = trapz(twoTheta, I_amorphous);
area_total = trapz(twoTheta, I_smooth);

CI = (area_cryst / area_total) * 100;

%% =========================================================================
% 5. PEAK DETECTION & CRYSTALLITE SIZE ESTIMATION (SCHERRER EQUATION)
% =========================================================================
% Scherrer Constants
K          = 0.9;     % Shape factor (Cubic/Spherical approximation)
lambda_xrd = 0.15406; % X-ray wavelength Cu-Kalpha in nanometers (nm)

min_height = max(I_crystalline) * 0.05; % Peak detection threshold (5% of main peak)

% Custom Local Peak Detection
locs_idx = find(I_crystalline(2:end-1) > I_crystalline(1:end-2) & ...
                I_crystalline(2:end-1) > I_crystalline(3:end) & ...
                I_crystalline(2:end-1) >= min_height) + 1;

pks      = I_crystalline(locs_idx);
numPeaks = length(pks);

fitResults = struct();
for i = 1:numPeaks
    idx = locs_idx(i);
    x0  = twoTheta(idx); % 2-Theta peak value (degrees)
    
    % Calculation of FWHM (Full Width at Half Maximum)
    half_height = pks(i) / 2;
    left_idx    = find(I_crystalline(1:idx) <= half_height, 1, 'last');
    right_idx   = idx + find(I_crystalline(idx:end) <= half_height, 1, 'first') - 1;
    
    if ~isempty(left_idx) && ~isempty(right_idx)
        fwhm_deg = twoTheta(right_idx) - twoTheta(left_idx);
    else
        fwhm_deg = 0.5;
    end
    
    % Degree to Radian Conversion for Scherrer Formula
    theta_rad = deg2rad(x0 / 2);   % Theta in radians
    beta_rad  = deg2rad(fwhm_deg); % FWHM (beta) in radians
    
    % Calculation of Crystallite Size D (Scherrer Equation)
    D_size_nm = (K * lambda_xrd) / (beta_rad * cos(theta_rad));
    
    fitResults.peak(i).peak_number         = i;
    fitResults.peak(i).center_2theta       = x0;
    fitResults.peak(i).fwhm_deg            = fwhm_deg;
    fitResults.peak(i).height              = pks(i);
    fitResults.peak(i).crystallite_size_nm = D_size_nm;
end

%% =========================================================================
% 6. DISPLAY RESULTS IN COMMAND WINDOW
% =========================================================================
fprintf('\n=========================================================================\n');
fprintf('       CRYSTALLINITY & CRYSTALLITE SIZE RESULTS (SCHERRER EQUATION)       \n');
fprintf('=========================================================================\n');
fprintf('Crystallinity Index (CI) : %.2f %%\n', CI);
fprintf('X-Ray Wavelength         : %.5f nm (Cu-K\\alpha)\n', lambda_xrd);
fprintf('-------------------------------------------------------------------------\n');
fprintf('Peak # | 2-Theta (deg) | FWHM (deg) | Height (a.u.) | Size D (nm)\n');
fprintf('-------------------------------------------------------------------------\n');
for i = 1:numPeaks
    fprintf('%6d | %13.3f | %10.3f | %13.1f | %11.2f\n', ...
        i, fitResults.peak(i).center_2theta, fitResults.peak(i).fwhm_deg, ...
        fitResults.peak(i).height, fitResults.peak(i).crystallite_size_nm);
end
fprintf('=========================================================================\n');

%% =========================================================================
% 7. MAIN VISUALIZATION (ANALYSIS & PEAK DETECTION)
% =========================================================================
figure('Name', 'XRD Analysis - Crystallinity & Peak Detection', 'Color', 'w', 'Position', [100, 100, 950, 650]);

% Subplot 1: Raw Data vs Amorphous Baseline
subplot(2, 1, 1);
plot(twoTheta, Intensity, 'Color', [0.7 0.7 0.7], 'LineWidth', 0.8, 'DisplayName', 'Raw Data'); hold on;
plot(twoTheta, I_smooth, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Smoothed Data');
plot(twoTheta, I_amorphous, 'r-', 'LineWidth', 1.8, 'DisplayName', 'Amorphous Baseline');
fill_x = [twoTheta; flipud(twoTheta)];
fill_y = [I_smooth; flipud(I_amorphous)];
fill(fill_x, fill_y, 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', 'Crystalline Area');
title(sprintf('Crystallinity Analysis | CI = %.2f%%', CI), 'FontSize', 11);
xlabel('2\theta (degrees)'); ylabel('Intensity (a.u.)');
legend('Location', 'northeast'); grid on;

% Subplot 2: Peak Detection & Crystallite Size Labels
subplot(2, 1, 2);
plot(twoTheta, I_crystalline, 'g-', 'LineWidth', 1.2, 'DisplayName', 'Crystalline Profile (I_{cryst})'); hold on;

peakCenters = [fitResults.peak.center_2theta];
peakHeights = [fitResults.peak.height];
plot(peakCenters, peakHeights, 'r*', 'MarkerSize', 8, 'LineWidth', 1.5, 'DisplayName', 'Detected Peaks');

% Text labels for position & crystallite size above peaks
for i = 1:numPeaks
    x0   = fitResults.peak(i).center_2theta;
    h    = fitResults.peak(i).height;
    d_nm = fitResults.peak(i).crystallite_size_nm;
    text(x0, h + max(I_crystalline)*0.03, sprintf('%.2f^\\circ\n(%.1f nm)', x0, d_nm), ...
        'FontSize', 8, 'HorizontalAlignment', 'center', 'Color', 'b', 'FontWeight', 'bold');
end

title('Peak Detection & Crystallite Size D (nm)', 'FontSize', 11);
xlabel('2\theta (degrees)'); ylabel('Intensity (a.u.)');
legend('Location', 'northeast'); grid on;
ylim([0, max(I_crystalline)*1.2]); % Add space for text labels above peaks

%% =========================================================================
% 8. EXPORT ANALYSIS RESULTS TO EXCEL FILE
% =========================================================================
output_excel = 'XRD_Analysis_Results_Complete.xlsx';
try
    T_peaks = struct2table(fitResults.peak);
    writetable(T_peaks, output_excel, 'Sheet', 'Peak_Analysis');
    
    T_summary = table(CI, area_cryst, area_amorf, area_total, ...
        'VariableNames', {'CI_Percent', 'Area_Crystalline', 'Area_Amorphous', 'Area_Total'});
    writetable(T_summary, output_excel, 'Sheet', 'Summary');
    fprintf('\nComplete analysis results successfully exported to: %s\n', output_excel);
catch ME
    fprintf('\nNote: Excel export skipped/failed (%s).\n', ME.message);
end
