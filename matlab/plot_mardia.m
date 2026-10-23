% plot_mardia.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

%% inputs
model = "PCR3BP"; 
a0 = 0.001;

%% parameters
if strcmp(model, "PCR3BP")
    d = 4; 
    m_epochs = [0, 0.931792004789438, 0.931792004789438 * 2, 0.931792004789438 * 3]; 
else
    m_epochs = [1.133];
    d = 6; 
end
nu = d * (d + 1) * (d + 2) / 6; 
N = 2 * d * d + 2 ^ d + 1;  
z1_s = chi2inv(1 - a0, nu); 
b1_s = z1_s * (6 * ((N + 1) * (d + 1) - 6)) / ((d + 1) * (N + 1) * (N + 3)); 
z2_s = norminv(1 - (a0 / 2));
b2_s = z2_s * sqrt(8 * d * (d + 2) / N) + (d * (d + 2));  
TU = 18913;
C = [238, 102, 119;  % Red
     68,  119, 170;  % Blue
     255, 140, 0;    % Orange
     34,  139, 34;]; % Green
C  = C/255;

%% figures
FILE = "../results/" + model + "/cukf_gbees/" + model + "_cukf_gbees_vars.mat"; load(FILE); 
f = figure; clf; f.Position = [250, 250, 1500, 800]; 
tiledlayout(2, 1, "Padding", "compact"); 
nexttile(1); hold on; grid on; box on; 
set(gca, 'FontName', 'Times','FontSize', 24, "LineWidth", 2); 
leg = legend(gca, 'Orientation', 'Horizontal', 'FontSize', 24, ...
             'FontName', 'times', 'Interpreter', 'latex', "LineWidth", 2);
leg.Layout.Tile = 'north'; 
ylabel("Estimated $\beta_{1, " + num2str(d) + "}$", "Interpreter", "latex", "FontSize", 32);
yline(b1_s, "-", "$b^*_{1, " + num2str(d) + "}(\alpha_0 =" + num2str(a0) + ")$", ...
              "FontSize", 28,"LabelHorizontalAlignment", "left", ...
              "LabelVerticalAlignment", "top", "LineWidth", 3, ...
              "Color", [0.5 0.5 0.5], "HandleVisibility", "off", "Interpreter", "latex"); 
xline(nan, ":", "LineWidth", 3, "Color", "k", "DisplayName", "Measurement update");
xline(nan, "--", "LineWidth", 3, "Color", "k", "DisplayName", "GVB Transition");

nexttile(2); hold on; grid on; box on; 
set(gca, 'FontName', 'Times','FontSize', 24, "LineWidth", 2); 
xlabel("$t$ (hours)", "Interpreter", "latex", "FontSize", 36);
ylabel("Estimated $\beta_{2, " + num2str(d) + "}$", "Interpreter", "latex", "FontSize", 32);
yline(b2_s, "-", "$b^*_{2, " + num2str(d) + "}(\alpha_0 =" + num2str(a0) + ")$", ...
              "FontSize", 28,"LabelHorizontalAlignment", "left", ...
              "LabelVerticalAlignment", "top", "LineWidth", 3, ...
              "Color", [0.5 0.5 0.5], "HandleVisibility", "off", "Interpreter", "latex"); 

for i = 1:numel(b_cukf_gbees)
    nexttile(1); 
    plot((t_cukf_gbees{i} + m_epochs(i)) .* TU / 3600, b_cukf_gbees{i}(:, 1), "LineWidth", 3, "Color", C(i, :), "DisplayName", "$p(\mathbf{x},t^{(" + num2str(i - 1) + ")})\,\,\,$");
    xline((t_cukf_gbees{i}(end) + m_epochs(i)) .* TU / 3600, "--", "LineWidth", 3, "Color", C(i, :), "HandleVisibility", "off");
    if i ~= 1
        xline(m_epochs(i) * TU / 3600, ":", "LineWidth", 3, "Color", "k", "HandleVisibility", "off");
    end
    nexttile(2); 
    plot((t_cukf_gbees{i} + m_epochs(i)) .* TU / 3600, b_cukf_gbees{i}(:, 2), "LineWidth", 3, "Color", C(i, :));
    xline((t_cukf_gbees{i}(end) + m_epochs(i)) .* TU / 3600, "--", "LineWidth", 3, "Color", C(i, :), "HandleVisibility", "off");
    if i ~= 1
        xline(m_epochs(i) * TU / 3600, ":", "LineWidth", 3, "Color", "k", "HandleVisibility", "off");
    end
end