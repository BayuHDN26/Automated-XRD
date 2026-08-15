clc;
clear;
close all;

%% =========================================================================
% 1. MEMBACA & PRE-PROCESSING DATA (ZERO-LEVEL SHIFT)
% =========================================================================
file_excel = '/MATLAB Drive/New Folder/Kak Nunu data/GO.xlsx';

if exist(file_excel, 'file')
    fprintf('Membaca file data XRD: %s ...\n', file_excel);
    data_excel = readtable(file_excel);
    twoTheta  = data_excel{:, 1};
    Intensity = data_excel{:, 2};
else
    fprintf('File tidak ditemukan. Menggunakan data simulasi dengan drift negatif...\n');
    twoTheta  = (15:0.05:70)';
    I_amorf   = 80 * exp(-((twoTheta - 22)/4.5).^2);
    I_cryst   = 45*exp(-((twoTheta-17.2)/0.4).^2) + 20*exp(-((twoTheta-20.2)/0.3).^2) + ...
                22*exp(-((twoTheta-22.5)/0.3).^2) + 18*exp(-((twoTheta-23.7)/0.3).^2) + ...
                15*exp(-((twoTheta-26.3)/0.4).^2) + 80*exp(-((twoTheta-32.1)/0.4).^2) + ...
                75*exp(-((twoTheta-34.7)/0.45).^2) + 160*exp(-((twoTheta-36.6)/0.5).^2) + ...
                12*exp(-((twoTheta-39.1)/0.4).^2) + 30*exp(-((twoTheta-57.0)/0.4).^2);
    drift     = -100 * (1 ./ (1 + exp(-(twoTheta - 42)/5))); % Drift ke nilai negatif
    noise     = 3 * randn(size(twoTheta));
    Intensity = I_amorf + I_cryst + drift + noise;
end

% A. ZERO-LEVEL OFFSET RESTORATION (Mengangkat spektrum ke I >= 0 secara fisis)
min_raw_val = min(Intensity);
if min_raw_val < 0
    I_raw_shifted = Intensity - min_raw_val;
else
    I_raw_shifted = Intensity;
end

% B. GAUSSIAN SMOOTHING
span     = round(0.01 * length(I_raw_shifted));
I_smooth = smoothdata(I_raw_shifted, 'gaussian', span);

%% =========================================================================
% 2. BASELINE AMORF (ASYMMETRIC LEAST SQUARES)
% =========================================================================
N      = length(I_smooth);
p      = 0.002;
lambda = 1e5;

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
% 3. KALKULASI PERSEN KRISTALINITAS (CI %)
% =========================================================================
area_cryst = trapz(twoTheta, I_crystalline);
area_amorf = trapz(twoTheta, I_amorphous);
CI = (area_cryst / (area_cryst + area_amorf)) * 100;

%% =========================================================================
% 4. DETEKSI PUNCAK & PERSAMAAN SCHERRER
% =========================================================================
K          = 0.9;
lambda_xrd = 0.15406; % nm (Cu-Kalpha)
min_height = max(I_crystalline) * 0.04;

locs_idx = find(I_crystalline(2:end-1) > I_crystalline(1:end-2) & ...
                I_crystalline(2:end-1) > I_crystalline(3:end) & ...
                I_crystalline(2:end-1) >= min_height) + 1;

pks      = I_crystalline(locs_idx);
numPeaks = length(pks);
fitResults = struct();

for i = 1:numPeaks
    idx = locs_idx(i);
    x0  = twoTheta(idx);
    
    % FWHM Calculation
    half_h = pks(i) / 2;
    l_idx  = find(I_crystalline(1:idx) <= half_h, 1, 'last');
    r_idx  = idx + find(I_crystalline(idx:end) <= half_h, 1, 'first') - 1;
    
    if ~isempty(l_idx) && ~isempty(r_idx)
        fwhm_deg = twoTheta(r_idx) - twoTheta(l_idx);
    else
        fwhm_deg = 0.4;
    end
    
    theta_rad = deg2rad(x0 / 2);
    beta_rad  = deg2rad(fwhm_deg);
    D_size_nm = (K * lambda_xrd) / (beta_rad * cos(theta_rad));
    
    fitResults.peak(i).center_2theta       = x0;
    fitResults.peak(i).fwhm_deg            = fwhm_deg;
    fitResults.peak(i).height              = pks(i);
    fitResults.peak(i).crystallite_size_nm = D_size_nm;
end

%% =========================================================================
% 5. VISUALISASI SIAP PUBLIKASI DENGAN ANOTASI HIBRIDA
% =========================================================================
figure('Name', 'XRD Analysis - Smart Hybrid Annotation', 'Color', 'w', 'Position', [100, 50, 1000, 750]);

% --- Subplot 1: Profil Total, Amorf, & Area Kristalin ---
subplot(2, 1, 1); hold on;
plot(twoTheta, I_raw_shifted, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8, 'DisplayName', 'Data Mentah (Shifted)');
plot(twoTheta, I_smooth, 'k-', 'LineWidth', 1.2, 'DisplayName', 'Smoothed Data');
plot(twoTheta, I_amorphous, 'r-', 'LineWidth', 2.0, 'DisplayName', 'Amorphous Baseline');

fill_x = [twoTheta; flipud(twoTheta)];
fill_y = [I_smooth; flipud(I_amorphous)];
fill(fill_x, fill_y, [0.2 0.7 0.3], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'DisplayName', 'Crystalline Area');

title(sprintf('Analisis Kristalinitas | CI = %.2f%% (Ac = %.1f, Aa = %.1f)', CI, area_cryst, area_amorf), ...
    'FontSize', 12, 'FontWeight', 'bold');
xlabel('2\theta (derajat)', 'FontWeight', 'bold');
ylabel('Intensitas (a.u.)', 'FontWeight', 'bold');
xlim([min(twoTheta) max(twoTheta)]);
ylim([0, max(I_smooth)*1.18]);
grid on; grid minor;
legend('Location', 'northeast', 'Box', 'off');
set(gca, 'LineWidth', 1.2, 'TickDir', 'in');

% --- Subplot 2: Profil Kristalin dengan Anotasi Adaptif ---
subplot(2, 1, 2); hold on;
plot(twoTheta, I_crystalline, 'g-', 'LineWidth', 1.3, 'DisplayName', 'Profil Kristalin (I_{cryst})');

peak_x = [fitResults.peak.center_2theta];
peak_y = [fitResults.peak.height];
plot(peak_x, peak_y, 'r.', 'MarkerSize', 8, 'DisplayName', 'Puncak Terdeteksi');

max_cryst = max(I_crystalline);

% Threshold: Puncak di atas 35% tinggi maksimum dikategorikan sebagai puncak tinggi
high_peak_threshold = 0.35 * max_cryst;

% Level vertikal 4 tingkat khusus untuk puncak-puncak rendah
low_tier_levels = [0.18, 0.36, 0.54, 0.72];
low_peak_count = 0;

for i = 1:numPeaks
    x0   = fitResults.peak(i).center_2theta;
    h0   = fitResults.peak(i).height;
    d_nm = fitResults.peak(i).crystallite_size_nm;
    label_str = sprintf('%.2f^\\circ\n(%.1f nm)', x0, d_nm);
    
    if h0 >= high_peak_threshold
        % =================================================================
        % A. PUNCAK TINGGI: Teks langsung di atas puncak (Tanpa Garis)
        % =================================================================
        text(x0, h0 + (max_cryst * 0.04), label_str, ...
            'FontSize', 8.5, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Color', [0 0.15 0.75], 'FontWeight', 'bold');
    else
        % =================================================================
        % B. PUNCAK RENDAH: Garis Penunjuk Bertingkat (4-Tier Staggering)
        % =================================================================
        low_peak_count = low_peak_count + 1;
        tier_idx = mod(low_peak_count - 1, length(low_tier_levels)) + 1;
        
        target_y = h0 + (max_cryst * low_tier_levels(tier_idx));
        
        % 1. Garis penunjuk putus-putus biru halus
        plot([x0, x0], [h0, target_y - max_cryst*0.015], 'b:', 'LineWidth', 0.9, 'HandleVisibility', 'off');
        plot(x0, target_y - max_cryst*0.015, 'b.', 'MarkerSize', 3.5, 'HandleVisibility', 'off');
        
        % 2. Teks informasi di ujung garis penunjuk
        text(x0, target_y, label_str, ...
            'FontSize', 7.5, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'Color', [0.1 0.3 0.65], 'FontWeight', 'bold');
    end
end

title('Deteksi Puncak & Ukuran Kristalit D (nm) [Smart Hybrid Annotation]', 'FontSize', 12, 'FontWeight', 'bold');
xlabel('2\theta (derajat)', 'FontWeight', 'bold');
ylabel('Intensitas (a.u.)', 'FontWeight', 'bold');
xlim([min(twoTheta) max(twoTheta)]);
ylim([0, max_cryst * 1.35]); % Ruang atas grafik diperlebar agar tier 4 tidak terpotong
grid on; grid minor;
legend('Location', 'northeast', 'Box', 'off');
set(gca, 'LineWidth', 1.2, 'TickDir', 'in');

%% =========================================================================
% 6. CETAK HASIL KE COMMAND WINDOW
% =========================================================================
fprintf('\n=======================================================================================\n');
fprintf('                HASIL KRISTALINITAS & UKURAN KRISTAL (PERSAMAAN SCHERRER)              \n');
fprintf('=======================================================================================\n');
fprintf('Crystallinity Index (CI) : %.2f %%\n', CI);
fprintf('Area Kristalin (Ac)      : %.2f\n', area_cryst);
fprintf('Area Amorf (Aa)          : %.2f\n', area_amorf);
fprintf('Panjang Gelombang X-Ray  : %.5f nm (Cu-K\\alpha)\n', lambda_xrd);
fprintf('--------------------------------------------------------------------\n');
fprintf('Peak # | 2-Theta (deg) | FWHM (deg) | Height (a.u.) | Size D (nm) |\n');
fprintf('--------------------------------------------------------------------\n');
for i = 1:numPeaks
    if fitResults.peak(i).height >= high_peak_threshold
        tipe_str = ' ';
    else
        tipe_str = ' ';
    end
    fprintf('%6d | %13.2f | %10.3f | %13.1f | %11.1f | %s\n', ...
        i, fitResults.peak(i).center_2theta, fitResults.peak(i).fwhm_deg, ...
        fitResults.peak(i).height, fitResults.peak(i).crystallite_size_nm, tipe_str);
end
fprintf('=======================================================================================\n');