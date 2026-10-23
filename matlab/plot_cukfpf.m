% plot_capstone_cukfpf.m

clear all; close all; clc; format long; 

%% inputs
load("colors.mat");
LU  = 389703; TU = 382981; U = [LU; LU; LU; (LU / TU) * 1000; (LU / TU) * 1000; (LU / TU) * 1000];
mu  = 1.215058446919780E-2; 
d = 6; 
a0 = 0.003; 
orbit = "NRHO";
len_str = "30d"; 
n_str = "5000"; 
HZ_str = "1-047"; 
didx = 50; 

%% b, HZ thresholds
nu = d * (d + 1) * (d + 2) / 6; 
N = 2 * d * d + 2 ^ d + 1;  
z1_t = chi2inv(1 - a0, nu); 
b1_t = z1_t * (6 * ((N + 1) * (d + 1) - 6)) / ((d + 1) * (N + 1) * (N + 3)); 
z2_t = norminv(1 - (a0 / 2));
b2_t = z2_t * sqrt(8 * d * (d + 2) / N) + (d * (d + 2));  

n = str2double(n_str); 
b = 1/(sqrt(2))*((2*d + 1)/4)^(1/(d + 4))*(n^(1/(d + 4))); % smoothing parameter
a = 1 + 2*b^2;
wb = (1 + b^2)*(1 + 3*b^2); 
mu_h = 1 - a^(- d/2)*(1 + d*b^2/a + (d*(d + 2)*(b^4))/(2*a^2)); % HZ mean
s2_h  = 2*(1 + 4*b^2)^(- d/2) + 2*a^(-d)*(1 + (2*d*b^4)/a^2 + (3*d*...
       (d + 2)*b^8)/(4*a^4)) - 4*wb^(-d/2)*(1 + (3*d*b^4)/(2*wb) + (d*...
       (d + 2)*b^8)/(2*wb^2)); % HZ variance
HZ_t = mu_h * (1 + (s2_h / mu_h^2))^(-1/2) * exp(norminv(1 - a0) * sqrt(log(1 + (s2_h / mu_h^2))));

%% reading in trajectories

% traj
FILE = "../results/" + orbit + "/traj/" + len_str + "/" + orbit + "_traj_vars.mat"; load(FILE);  
t_traj = t_traj .* TU / 86400;
x_traj = x_traj .* U'; 
x_traj_M = x_traj; x_traj_M(:, 1) = x_traj_M(:, 1) - LU * (1 - mu);
dist_center = vecnorm(x_traj_M(:, 1:3)'); 
tCA = find(islocalmin(dist_center));
tFA = find(islocalmax(dist_center)); 
max_dist = max(dist_center); 

% cukfpf
FILE = "../results/" + orbit + "/cukfpf/" + len_str + "_" + n_str + "_HZ_" + HZ_str + "/" + orbit + "_cukfpf_vars.mat"; load(FILE); 
for i = 1:numel(t_cukfpf_cukf)
    t_cukfpf_cukf{i} = t_cukfpf_cukf{i} .* TU / 86400;
    x_cukfpf_cukf{i} = x_cukfpf_cukf{i} .* permute(repmat(U, 1, size(x_cukfpf_cukf{i}, 1), size(x_cukfpf_cukf{i}, 2)), [2 3 1]); 
    t_0 = t_cukfpf_cukf{i}(1, 1); t_f = t_cukfpf_cukf{i}(end, 1); 
    idx_0 = find(t_traj == t_0); idx_f = find(t_traj == t_f); 
    x_traj_i = x_traj(idx_0:idx_f, :); 
    x_traj_i = reshape(x_traj_i, 1, size(x_traj_i, 1), size(x_traj_i, 2));
    dx_cukfpf_cukf{i} = x_cukfpf_cukf{i} - x_traj_i; 
end

for i = 1:numel(t_cukfpf_pf)
    t_cukfpf_pf{i} = t_cukfpf_pf{i} .* TU / 86400;
    x_cukfpf_pf{i} = x_cukfpf_pf{i} .* permute(repmat(U, 1, size(x_cukfpf_pf{i}, 1), size(x_cukfpf_pf{i}, 2)), [2 3 1]); 
    t_0 = t_cukfpf_pf{i}(1); t_f = t_cukfpf_pf{i}(end); 
    idx_0 = find(t_traj == t_0); idx_f = find(t_traj == t_f); 
    x_traj_i = x_traj(idx_0:idx_f, :); 
    x_traj_i = reshape(x_traj_i, 1, size(x_traj_i, 1), size(x_traj_i, 2));
    dx_cukfpf_pf{i} = x_cukfpf_pf{i} - x_traj_i; 
end

idxs = {}; 
for i = 1:numel(t_b)
    t_b{i} = t_b{i} .* TU / 86400; 
    % b{i}(:, 1) = b{i}(:, 1) ./ b1_t;
    % b{i}(:, 2) = (b{i}(:, 2) ./ b2_t) - 1;

    idx_i = 1:didx:numel(t_b{i})'; 
    for j = 1:numel(t_cukfpf_pf)

        t_j = t_cukfpf_pf{j}(1); 
        cross_idx = find(t_b{i} == t_j); 
        if ~isempty(cross_idx)
            idx_i(end + 1) = cross_idx(:); 
        end
        t_j = t_cukfpf_pf{j}(end); 
        cross_idx = find(t_b{i} == t_j); 
        if ~isempty(cross_idx)
            idx_i(end + 1) = cross_idx(:); 
        end

    end
    idxs{i} = sort(unique(idx_i)); 
end
for i = 1:numel(t_HZ), t_HZ{i} = t_HZ{i} .* TU / 86400; end

% pf
FILE = "../results/" + orbit + "/pf/" + len_str + "_" + n_str + "/" + orbit + "_pf_vars.mat"; load(FILE); 
for i = 1:numel(x_pf)
    t_pf{i} = t_pf{i} .* TU / 86400;
    x_pf{i} = x_pf{i} .* permute(repmat(U, 1, size(x_pf{i}, 1), size(x_pf{i}, 2)), [2 3 1]); 
    t_0 = t_pf{i}(1); t_f = t_pf{i}(end); 
    idx_0 = find(t_traj == t_0); idx_f = find(t_traj == t_f); 
    x_traj_i = x_traj(idx_0:idx_f, :); 
    x_traj_i = reshape(x_traj_i, 1, size(x_traj_i, 1), size(x_traj_i, 2));
    dx_pf{i} = x_pf{i} - x_traj_i; 
end

% cukf
FILE = "../results/" + orbit + "/cukf/" + len_str + "/" + orbit + "_cukf_vars.mat"; load(FILE); 
for i = 1:numel(x_cukf)
    t_cukf{i} = t_cukf{i} .* TU / 86400;
    x_cukf{i} = x_cukf{i} .* permute(repmat(U, 1, size(x_cukf{i}, 1), size(x_cukf{i}, 2)), [2 3 1]); 
    t_0 = t_cukf{i}(1); t_f = t_cukf{i}(end); 
    idx_0 = find(t_traj == t_0); idx_f = find(t_traj == t_f); 
    x_traj_i = x_traj(idx_0:idx_f, :); 
    x_traj_i = reshape(x_traj_i, 1, size(x_traj_i, 1), size(x_traj_i, 2));
    dx_cukf{i} = x_cukf{i} - x_traj_i; 
end

%% initialize figures
% 
% f = figure; clf; f.Position = [323,249,1600,1000];  
% t = tiledlayout(3, 1, "TileSpacing", "compact", "Padding", "compact");
% 
% nexttile(1); hold on; box on; grid on; ax = gca;
% set(ax, 'FontName', 'times', 'FontSize', 28, 'LineWidth', 2);
% ylabel(ax, "$\delta x$ (km)", 'FontSize', 36, "Interpreter", "latex");
% xlim([0, t_traj(end)]); 
% % ylim([-130, 130]);
% 
% nexttile(2); hold on; box on; grid on; ax = gca;
% set(ax, 'FontName', 'times', 'FontSize', 28, 'LineWidth', 2);
% if ~isempty(HZ)
%     yyaxis left; hold on; 
%     ax.YAxis(1).Color = hanred;
% end
% % yline(1, "--", "$\tilde{b}^*_{i}$", ...
% yline(b1_t, "--", ...
%     "FontSize", 24, "LabelHorizontalAlignment", "left", ...
%     "LabelVerticalAlignment", "top", "LineWidth", 3, "LineStyle", "--", ...
%     "Color", hanred, "HandleVisibility", "off", "Interpreter", "latex"); 
% yline(b2_t, "--", ...
%     "FontSize", 24, "LabelHorizontalAlignment", "left", ...
%     "LabelVerticalAlignment", "top", "LineWidth", 3, "LineStyle", "--", ...
%     "Color", hanred, "HandleVisibility", "off", "Interpreter", "latex"); 
% ylabel("Normalized $\tilde{b}_{i}$", 'FontSize', 32, "Interpreter", "latex");
% xlim([0, t_traj(end)]); 
% % b2_min = inf; b2_max = -inf; 
% % for i = 1:numel(b)
% %     b2_min = min(b2_min, min(b{i}(:, 2))); b2_max = max(b2_max, max(b{i}(:, 2)));
% % end
% % db2 = b2_max - b2_min;
% % ylim([b2_min - 0.05 * db2, b2_max + 0.1 * db2]);
% % ylim([-0.5, 3.25]); 
% 
% if ~isempty(HZ)
%     yyaxis right; hold on; 
%     ax.YAxis(2).Color = hanblue;
%     % yline(HZ_t, "--", "$\mathrm{HZ}^*$", ...
%     yline(HZ_t, "--", ...
%           "FontSize", 28, "LabelHorizontalAlignment", "right", ...
%           "LabelVerticalAlignment", "bottom", "LineWidth", 3, "LineStyle", "--", ...
%           "Color", hanblue, "HandleVisibility", "off", "Interpreter", "latex", "HandleVisibility", "off"); 
%     ylabel("HZ", 'FontSize', 36, "Interpreter", "latex");
%     HZ_min = inf; HZ_max = -inf; 
%     for i = 1:numel(HZ)
%         HZ_min = min(HZ_min, min(HZ{i})); HZ_max = max(HZ_max, max(HZ{i}));
%     end
%     dHZ = HZ_max - HZ_min;
%     ylim([-10, HZ_max]);
% end
% 
% nexttile(3); hold on; box on; grid on; ax = gca;
% set(ax, 'FontName', 'times', 'FontSize', 28, 'LineWidth', 2);
% xlabel(ax, "$t$ (days)", 'FontSize', 36, "Interpreter", "latex");
% ylabel("$r_{\mathrm{M}}$ ($10^4$ km)", 'FontSize', 36, "Interpreter", "latex");
% xlim([0, t_traj(end)]); 
% ylim([0, 1.1 * max_dist / 1e4]);
% 
% % close approaches
% for j = 1:numel(tCA) - 1
%     for k = 1:3
%         nexttile(k); 
%         xline(t_traj(tCA(j)), ":", "LineWidth", 3, "Color", "k", "HandleVisibility", "off", "Interpreter", "latex"); 
%         % xline(t_traj(j), ":", "$t_{\mathrm{CA}}$", "LineWidth", 3, "Color", "k", ...
%         %   "LabelHorizontalAlignment", "right", "LabelOrientation", "aligned", ...
%         %   "FontSize", 22, "FontName", "times", "HandleVisibility", "off", "Interpreter", "latex");
%     end
% end
% 
% %% trajectories
% for i = 1:numel(x_cukfpf_pf)
%     nexttile(1); 
%     plot(t_cukfpf_pf{i}, squeeze(dx_cukfpf_pf{i}(:, :, 1))', "LineWidth", 3, "Color", (1 - 0.6)* hanblue + 0.6*[1 1 1]); 
% 
%     for j = 1:3
%         nexttile(j); 
%         if (i ~= (numel(x_cukfpf_pf) - 1))&&(i ~= (numel(x_cukfpf_pf)))
%             xline(t_cukfpf_pf{i}(end), "Color", "k", "LineStyle", "--", "LineWidth", 3, "HandleVisibility", "off");
%         end
%         if i ~= numel(x_cukfpf_pf)
%             xline(t_cukfpf_pf{i}(1), "Color", "k", "LineStyle", "--", "LineWidth", 3, "HandleVisibility", "off");
%         end
%     end
% end
% 
% for i = 1:numel(x_cukfpf_cukf)
%     nexttile(1); 
%     plot(t_cukfpf_cukf{i}, squeeze(dx_cukfpf_cukf{i}(:, :, 1))', "LineWidth", 2, "Color", hanred);
% end
% 
% for i = 1:3
%     nexttile(i); 
%     xline(t_cukf{1}(end), "Color", [0 0 0], "LineStyle", "-.", "LineWidth", 3, "HandleVisibility", "off"); 
% end
% 
% %% HZ-b
% nexttile(2); legend("Location", "northwest", "Interpreter", "latex"); 
% for i = 1:numel(t_b)
%     if ~isempty(HZ), yyaxis left; end
%     plot(t_b{i}(idxs{i}), b{i}(idxs{i}, 1), "Color", hanred, "LineWidth", 3, "LineStyle", "-", "Marker", "x", "MarkerSize", 15, "DisplayName", "$i=1$ "); 
%     plot(t_b{i}(idxs{i}), b{i}(idxs{i}, 2), "Color", hanred, "LineWidth", 3, "LineStyle", "-", "Marker", "o", "MarkerSize", 15, "DisplayName", "$i=2$ "); 
% end
% 
% for i = 1:numel(t_HZ)
%     if ~isempty(HZ), yyaxis right; end
%     plot(t_HZ{i}, HZ{i}, "Color", hanblue, "LineWidth", 4, "LineStyle", "-", "HandleVisibility", "off"); 
% end
% 
% nexttile(3); 
% plot(t_traj, dist_center ./ 1e4, "Color", "k", "LineWidth", 3, "LineStyle", "-", "DisplayName", "Nominal");
% 
% clear L; clear LH;
% LH(1) = plot(nan, nan,  "Color", "k", "Marker", "none", "LineStyle", "-", "LineWidth", 3);
% L{1} = "Nominal";
% LH(2) = plot(nan, nan,  "Color", hanred, "Marker", "none", "LineStyle", "-", "LineWidth", 3);
% L{2} = "CUKF6";
% LH(3) = plot(nan, nan,  "Color", hanblue, "Marker", "none", "LineStyle", "-", "LineWidth", 3);
% L{3} = "PF";
% LH(4) = plot(nan, nan,  "Color", "k", "LineWidth", 2, "LineStyle", "--", "Marker", "none");
% L{4} = "GVB transition";
% LH(5) = plot(nan, nan,  "Color", "k", "LineWidth", 2, "LineStyle", "-.", "Marker", "none");
% L{5} = "Measurement update";
% LH(6) = plot(nan, nan,  "Color", "k", "LineWidth", 2, "LineStyle", ":", "Marker", "none");
% L{6} = "$t_{\mathrm{CA}}$";
% leg = legend(LH, L, "FontSize", 28, "FontName", "times", "Interpreter", "latex",...
%              "LineWidth", 1, "Orientation", "horizontal");
% leg.Layout.Tile = 'north'; 
% drawnow; 

%% distributions
var_lbls = {"$\delta x$ (km)", "$\delta y$ (km)", "$\delta z$ (km)", ...
            "$\delta v_x$ (m/s)", "$\delta v_y$ (m/s)", "$\delta v_z$ (m/s)"};

% t meas -
clear p; clear L; clear LH;
f = figure; clf; f.Position = [276,31,800,800]; hold on;
pause(1.0); 
p.type = "scatter"; p.color = [0.7 0.7 0.7]; p.hist = 1; p.m_alpha = 0.25; 
plot_corner_pdf(squeeze(dx_pf{1}(:, end, :)), 'lbls', var_lbls, 'p', p); 
p.type = "cukf"; p.color = hanred; p.ls = "-."; p.hold_lims = 0; 
plot_corner_pdf(squeeze(dx_cukf{1}(:, end, :)), 'p', p, 'P', W_cukf{1});
p.type = "scatter_grid"; p.color = hanblue; p.ls = "-"; 
plot_corner_pdf(squeeze(dx_cukfpf_pf{end - 1}(:, end, :)), 'P', W_cukfpf_pf{end - 1}, 'isovalue', 1, 'p', p); 
LH(1) = scatter(nan, nan, 100, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
L{1} = "PF";
LH(2) = plot(nan, nan,  "Color", hanred, "LineWidth", 2, "LineStyle", "-.");
L{2} = "CUKF6";
LH(3) = plot(nan, nan,  "Color", hanblue, "LineWidth", 2);
L{3} = "CUKF6/\textbf{PF}";
leg = legend(LH, L, "FontSize", 18, "FontName", "times", "Interpreter", "latex",...
             "LineWidth", 1, "Orientation", "horizontal", "EdgeColor", "w", "TextColor", "w");
leg.Layout.Tile = 'south'; 
sgtitle(sprintf('$t^- = %.4f$ days', t_pf{1}(end)), "Color", "w", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
drawnow; pause(0.1); 

% t meas +
clear p; clear L; clear LH;
f = figure; clf; f.Position = [276,31,800,800]; hold on;
p.type = "scatter"; p.color = [0.7 0.7 0.7]; p.hist = 1; p.m_alpha = 0.25; 
plot_corner_pdf(squeeze(dx_pf{2}(:, 1, :)), 'lbls', var_lbls, 'p', p); 
p.type = "cukf"; p.color = hanred; p.ls = "-."; p.hold_lims = 0; 
plot_corner_pdf(squeeze(dx_cukf{2}(:, 1, :)), 'p', p, 'P', W_cukf{2});
p.type = "scatter_grid"; p.color = hanblue; p.ls = "-"; p.Nx_bins = 10; p.Ny_bins = 10; 
plot_corner_pdf(squeeze(dx_cukfpf_pf{end}(:, 1, :)), 'P', W_cukfpf_pf{end}, 'isovalue', 1, 'p', p); 
LH(1) = scatter(nan, nan, 100, 'o', 'MarkerFaceColor', [0.7 0.7 0.7], 'MarkerEdgeColor', 'none');
L{1} = "PF";
LH(2) = plot(nan, nan,  "Color", hanred, "LineWidth", 2, "LineStyle", "-.");
L{2} = "CUKF6";
LH(3) = plot(nan, nan,  "Color", hanblue, "LineWidth", 2);
L{3} = "CUKF6/\textbf{PF}";
leg = legend(LH, L, "FontSize", 18, "FontName", "times", "Interpreter", "latex",...
             "LineWidth", 1, "Orientation", "horizontal", "EdgeColor", "w", "TextColor", "w");
leg.Layout.Tile = 'south'; 
sgtitle(sprintf('$t^+ = %.4f$ days', t_pf{2}(1)), "Color", "k", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
drawnow; pause(0.1); 