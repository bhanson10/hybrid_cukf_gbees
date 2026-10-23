% GVT_diff_SatEnc.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

close all; clc; 

%% parameters
x_var = "periods"; % "TU", "hrs", "periods"
fams  = {"DPO", "DRO", "L1_HALO_NORTHERN", "L1_HALO_SOUTHERN",...
         "L2_HALO_NORTHERN", "L2_HALO_SOUTHERN", "L2_BUTTERFLY_NORTHERN",...
         "L2_BUTTERFLY_SOUTHERN", "DRAGONFLY_NORTHERN", "DRAGONFLY_SOUTHERN"};
N = 50; 

%% metaparameters
TU = 18913;
DIR  = "../c/outputs"; 
load hancolors.mat;

%% initializing figure
f1 = figure(1); clf; f1.Position = [100,100,1300,650]; 
tiledlayout(1, 2, "Padding", "compact"); 
nexttile(1); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; 
xlabel("MC ($n=5000$) GVT -- HZ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
ylabel("UT ($N=13$) GVT -- NED (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 

c = {};
names = {};
for i = 1:numel(fams)
   
    if strcmp(fams{i}, "DPO")
        name = "Distant Prograde"; 
    elseif strcmp(fams{i}, "DPO_2")||strcmp(fams{i}, "DPO_3")
        name = "Distant Prograde ($+0.1T$)"; 
    elseif strcmp(fams{i}, "DRAGONFLY_NORTHERN")
        name = "Northern Dragonfly"; 
    elseif strcmp(fams{i}, "DRAGONFLY_SOUTHERN")
        name = "Southern Dragonfly"; 
    elseif strcmp(fams{i}, "DRO")
        name = "Distant Retrograde"; 
    elseif strcmp(fams{i}, "L1_HALO_NORTHERN")
        name = "$L_1$ Northern Halo"; 
    elseif strcmp(fams{i}, "L1_HALO_SOUTHERN")
        name = "$L_1$ Southern Halo"; 
    elseif strcmp(fams{i}, "L2_BUTTERFLY_NORTHERN")
        name = "Northern Butterfly"; 
    elseif strcmp(fams{i}, "L2_BUTTERFLY_SOUTHERN")
        name = "Southern Butterfly"; 
    elseif strcmp(fams{i}, "L2_HALO_NORTHERN")
        name = "$L_2$ Northern Halo"; 
    elseif strcmp(fams{i}, "L2_HALO_SOUTHERN")
        name = "$L_2$ Southern Halo"; 
    end
    
    c{i} = C(i, :);
    scatter(nan, nan, 50, c{i}, "filled", "o", "DisplayName", name);
end
leg = legend("FontSize", 22, "FontName", "times", "Interpreter", "latex",...
         "LineWidth", 1, "Orientation", "horizontal", "NumColumns", 5);
leg.Layout.Tile = 'north';

% nexttile(2); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; axis square; 
% xlabel("MC ($n=5000$) GVT -- HZ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
% ylabel("CUT6 ($N=137$) GVT -- $\tilde{b}_1$, $\tilde{b}_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 

nexttile(2); hold on; set(gca, "LineWidth", 2, 'FontName' , 'Times','FontSize', 24); grid on; box on; 
xlabel("MC ($n=5000$) GVT -- $b_1$, $b_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 
ylabel("CUT6 ($N=137$) GVT -- $\tilde{b}_1$, $\tilde{b}_2$ (" + x_var + ")", "FontSize", 24, "Interpreter", "latex"); 

%% Get GVT values
gvt_map = get_GVT(DIR + "/gauss_detect_SatEnc_output.log"); 

%% Plot GVT values
HZ_GVT = [];
Mardia_GVT = [];
NED_GVT = [];
CUT_GVT = [];
for i = 1:numel(fams)

    fam_data = gvt_map(fams{i}); 

    HZ_GVT_i = fam_data.values(:, 1); 
    Mardia_GVT_i = fam_data.values(:, 2); 
    NED_GVT_i = fam_data.values(:, 3); 
    CUT_GVT_i = fam_data.values(:, 4); 

    if strcmp(x_var, "hrs")
        HZ_GVT_i = HZ_GVT_i .* TU / 3600; 
        Mardia_GVT_i = Mardia_GVT_i .* TU / 3600; 
        NED_GVT_i = NED_GVT_i .* TU / 3600; 
        CUT_GVT_i = CUT_GVT_i .* TU / 3600; 
    elseif strcmp(x_var, "periods")
        HZ_GVT_i = HZ_GVT_i ./ fam_data.periods; 
        Mardia_GVT_i = Mardia_GVT_i ./ fam_data.periods; 
        NED_GVT_i = NED_GVT_i ./ fam_data.periods;
        CUT_GVT_i = CUT_GVT_i ./ fam_data.periods; 
    end

    nexttile(1); 
    scatter(HZ_GVT_i, NED_GVT_i, 70, "filled", "o", "MarkerFaceColor", c{i}, "HandleVisibility", "off");

    % nexttile(2); 
    % scatter(HZ_GVT_i, CUT_GVT_i, 70, "filled", "o", "MarkerFaceColor", c{i}, "HandleVisibility", "off");

    nexttile(2); 
    scatter(Mardia_GVT_i, CUT_GVT_i, 70, "filled", "o", "MarkerFaceColor", c{i}, "HandleVisibility", "off");

    HZ_GVT = [HZ_GVT; HZ_GVT_i(:)];
    Mardia_GVT = [Mardia_GVT; Mardia_GVT_i(:)];
    NED_GVT = [NED_GVT; NED_GVT_i(:)];
    CUT_GVT = [CUT_GVT; CUT_GVT_i(:)];
end

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

if strcmp(x_var, "periods")
    nexttile(1); xlim([0, 0.7]); ylim([0, 0.7]);
    nexttile(2); xlim([0, 0.7]); ylim([0, 0.7]);
end

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gvt_map = get_GVT(filename)
    % Open log file for reading
    fileID = fopen(filename, 'r');
    
    % Map to hold aggregated data per family
    % Structure per key: struct('Period', [], 'Values', [])
    gvt_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    
    % Read line-by-line
    while ~feof(fileID)
        line = strtrim(fgetl(fileID));
        
        % Skip empty lines or header info without data
        if isempty(line) || startsWith(line, 'GVT') || ~contains(line, '-')
            continue;
        end
        
        % Split the line into [Header, Data] around the '-' dash
        parts = strsplit(line, '-');
        header_part = parts{1};
        data_part   = parts{2};
        
        % Extract Family and Period from header using regex
        fam_match = regexp(header_part, 'Family:\s*([^,]+)', 'tokens');
        per_match = regexp(header_part, 'Period \(TU\):\s*([\d\.]+)', 'tokens');
        
        % Extract numerical values from the right side of the dash
        vals = sscanf(data_part, '%f %f %f %f')';
        
        if ~isempty(fam_match) && ~isempty(per_match) && length(vals) == 4
            fam_name = strtrim(fam_match{1}{1});
            period  = str2double(per_match{1}{1});
            
            % Retrieve or initialize struct for this family
            if isKey(gvt_map, fam_name)
                entry = gvt_map(fam_name);
            else
                entry = struct('periods', [], 'values', []);
            end
            
            % Append period and values
            entry.periods(end+1, 1) = period;
            entry.values(end+1, :)  = vals;
            
            gvt_map(fam_name) = entry;
        end
    end
    fclose(fileID);
    
    % Create dynamic variables in the workspace for each family
    families = keys(gvt_map);
    for i = 1:length(families)
        fam_name = families{i};
        
        % Sanitize family name to create a valid MATLAB variable name
        var_name = genvarname(fam_name); 
        
        % Store structured data into workspace dynamic variable
        assignin('caller', var_name, gvt_map(fam_name));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%