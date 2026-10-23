% GVT_diff_NEA.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

close all; clc; 

%% parameters
x_var = "periods"; % "days", "periods"

%% metaparameters
DIR  = "../c/outputs"; 
load hancolors.mat;

%% initializing figure
f1 = figure(1); clf; f1.Position = [100,100,1300,650];
tiledlayout(1, 2, "Padding", "compact"); 
nexttile(1); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; axis square; 
xlabel("MC ($n=5000$) GVT -- HZ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
ylabel("UT ($N=13$) GVT -- NED (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
set(gca, 'ColorScale', 'log');

% nexttile(2); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; axis square; 
% xlabel("MC ($n=5000$) GVT -- HZ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
% ylabel("CUT6 ($N=137$) GVT -- $\tilde{b}_1$, $\tilde{b}_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 

nexttile(2); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; axis square; 
xlabel("MC ($n=5000$) GVT -- $b_1$, $b_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
ylabel("CUT6 ($N=137$) GVT -- $\tilde{b}_1$, $\tilde{b}_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 

set(gca, 'ColorScale', 'log');
cmap = colormap(slanCM("bubblegum")); 
cb = colorbar('eastoutside');  
cb.Label.String = '$||\Sigma^{(0)}||_F$';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 32;  

%% Get GVT values
gvt_list = get_GVT(DIR + "/gauss_detect_NEA_output_500.log"); 
N = size(gvt_list, 1); 

%% Plot GVT values
T = gvt_list(:, 1); 
A = gvt_list(:, 2); 
HZ_GVT = gvt_list(:, 3); 
Mardia_GVT = gvt_list(:, 4); 
NED_GVT = gvt_list(:, 5); 
CUT_GVT = gvt_list(:, 6); 
if strcmp(x_var, "periods")
    HZ_GVT = HZ_GVT ./ T; 
    Mardia_GVT = Mardia_GVT ./ T; 
    NED_GVT = NED_GVT ./ T; 
    CUT_GVT = CUT_GVT ./ T; 
end
   
nexttile(1); 
scatter(HZ_GVT, NED_GVT, 70, A, "filled", "o", "HandleVisibility", "off");
clim([min(A), max(A)]);

% nexttile(2); 
% scatter(HZ_GVT_i, CUT_GVT_i, 70, A, "filled", "o", "HandleVisibility", "off");

nexttile(2); 
scatter(Mardia_GVT, CUT_GVT, 70, A, "filled", "o", "HandleVisibility", "off");
clim([min(A), max(A)]);

% Add 1:1 reference line and compute R^2 relative to y = x
for k = 1:2
    nexttile(k); 
    
    % Get data points for tile k
    if k == 1
        x = HZ_GVT;
        y = NED_GVT;
    % elseif k == 2
    %     x = HZ_GVT;
    %     y = CUT_GVT;
    elseif k == 2
        x = Mardia_GVT;
        y = CUT_GVT;
    end
    
    % Calculate R^2 with respect to 1:1 line (y = x)
    SS_res = sum((y - x).^2);               % Sum of squared residuals relative to y = x
    SS_tot = sum((y - mean(y)).^2);         % Total sum of squares
    R2_1to1 = 1 - (SS_res / SS_tot);
    
    % Draw 1:1 reference line
    xL = xlim; yL = ylim; 
    p1 = plot([-100, 100], [-100, 100], "Color", [0.6 0.6 0.6], "Linestyle", "--", "LineWidth", 3, "HandleVisibility", "off");
    uistack(p1, "bottom"); 
    xlim(xL); ylim(yL); 

    % Place text in top-left using normalized axes coordinates (5% from left, 90% from bottom)
    text(0.05, 0.90, sprintf('$R^2 = %.4f$', R2_1to1), ...
         'Units', 'normalized', ...
         'Interpreter', 'latex', ...
         'FontSize', 22, ...
         'FontName', 'Times', ...
         'BackgroundColor', 'w', ...      % White background box so points behind it don't overlap
         'EdgeColor', 'none');
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gvt_list = get_GVT(filename)
    % Open log file for reading
    fileID = fopen(filename, 'r');
    
    % List to hold aggregated data per family
    gvt_list = []; 
    
    % Read line-by-line
    while ~feof(fileID)
        line = strtrim(fgetl(fileID));
        
        % Skip empty lines or header info without data
        if isempty(line) || startsWith(line, 'GVT') || ~contains(line, '-')
            continue;
        end
        
        % Split the line into [Header, Data] around the '-' dash
        parts = strsplit(line, '=');
        header_part = parts{1};
        data_part   = parts{2};
        
        % Extract Family and Period from header using regex
        NEA = regexp(header_part, 'NEA:\s*([\d\.]+)', 'tokens');
        NEA = str2double(NEA{1}{1});
        period = regexp(header_part, 'Period \(d\):\s*([\d\.]+)', 'tokens');
        period = str2double(period{1}{1});
        token = regexp(header_part, 'Frob.\s*norm:\s*(\S+)', 'tokens', 'once');
        frob_norm = str2double(token{1});
            
        % Extract numerical values from the right side of the dash
        vals = sscanf(data_part, '%f %f %f %f')';
        
        gvt_list = [gvt_list; period frob_norm vals NEA];
    end
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%