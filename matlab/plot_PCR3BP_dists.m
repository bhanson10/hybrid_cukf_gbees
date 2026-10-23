% plot_PCR3BP_dists.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

%% parameters
model = "PCR3BP"; 
filters = {"ukf", "cukf"};
mu = 1.901109735892602E-7;
LU = 238529;
TU = 18913;
U = [LU; LU; (LU / TU); (LU / TU)];
sec_r = 252.1;
C = [238, 102, 119;  % Red
     68,  119, 170;  % Blue
     255, 140, 0;    % Orange
     34,  139, 34;]; % Green
C  = C/255;

%% nominal trajectory
FILE = "../results/PCR3BP_rv.txt"; 
[t_traj, x_traj] = parse_traj(FILE); 
t_traj = t_traj .* TU / 86400;
x_traj = x_traj .* U'; 

axs = initialize_figures(filters, mu, LU, sec_r); 
for nm = 1:numel(filters)
    nexttile(nm); 
    plot(x_traj(:,1), x_traj(:,2),"k-","LineWidth", 2,"HandleVisibility","off");
    drawnow;
    nexttile(numel(filters) + nm); 
    plot(x_traj(:,3), x_traj(:,4),"k-","LineWidth", 2,"HandleVisibility","off");
    drawnow;
end

for nm = 1:numel(filters)
    clear p; 
    switch filters{nm}
        case "pf"
            %% PF
            alphavalue = 0.1; 
            FILE = "../results/" + model + "/pf/" + model + "_pf_vars.mat"; load(FILE); 
            for i = 1:numel(x_pf)
                t_pf{i} = t_pf{i} .* TU / 86400;
                x_pf{i} = x_pf{i} .* permute(repmat(U, 1, size(x_pf{i}, 1), size(x_pf{i}, 2)), [2 3 1]); 
            end

            count = 1;
            for i=1:numel(x_pf)

                p.color = C(i, :);

                x_pf_i = x_pf{i}; 

                nexttile(nm); 
                scatter(squeeze(x_pf_i(:, 1, 1)), squeeze(x_pf_i(:, 1, 2)), 10, p.color, "filled", ...
                    "MarkerEdgeAlpha", alphavalue, "MarkerFaceAlpha", alphavalue);
                scatter(squeeze(x_pf_i(:, end, 1)), squeeze(x_pf_i(:, end, 2)), 10, p.color, "filled", ...
                    "MarkerEdgeAlpha", alphavalue, "MarkerFaceAlpha", alphavalue);
                nexttile(numel(filters) + nm); 
                scatter(squeeze(x_pf_i(:, 1, 3)), squeeze(x_pf_i(:, 1, 4)), 10, p.color, "filled", ...
                    "MarkerEdgeAlpha", alphavalue, "MarkerFaceAlpha", alphavalue);
                scatter(squeeze(x_pf_i(:, end, 3)), squeeze(x_pf_i(:, end, 4)), 10, p.color, "filled", ...
                    "MarkerEdgeAlpha", alphavalue, "MarkerFaceAlpha", alphavalue);

                drawnow; 
            end
        case "ukf"
            %% UKF
            FILE = "../results/" + model + "/ukf/" + model + "_ukf_vars.mat"; load(FILE); 
            for i = 1:numel(x_ukf)
                t_ukf{i} = t_ukf{i} .* TU / 86400;
                x_ukf{i} = x_ukf{i} .* permute(repmat(U, 1, size(x_ukf{i}, 1), size(x_ukf{i}, 2)), [2 3 1]); 
            end

            count = 1;
            for i=1:numel(x_ukf)

                p.color = C(i, :);

                x_ukf_i = x_ukf{i}; W_ukf_i = W_ukf{i};
                x_ukf_i0 = squeeze(x_ukf_i(:, 1, :)); 
                x_ukf_if = squeeze(x_ukf_i(:, end, :)); 
                [mu0, S0] = calc_moments_ukf(x_ukf_i0, W_ukf_i);
                [muf, Sf] = calc_moments_ukf(x_ukf_if, W_ukf_i);
                
                nexttile(nm); 
                scatter(x_ukf_i0(:, 1), x_ukf_i0(:, 2), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(1:2), S0(1:2, 1:2), "p", p);
                scatter(x_ukf_if(:, 1), x_ukf_if(:, 2), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(muf(1:2), Sf(1:2, 1:2), "p", p);

                nexttile(numel(filters) + nm); 
                scatter(x_ukf_i0(:, 3), x_ukf_i0(:, 4), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(3:4), S0(3:4, 3:4), "p", p);
                scatter(x_ukf_if(:, 3), x_ukf_if(:, 4), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(muf(3:4), Sf(3:4, 3:4), "p", p);

                drawnow; 
            end
        case "cukf"
            %% CUKF
            FILE = "../results/" + model + "/cukf/" + model + "_cukf_vars.mat"; load(FILE); 
            for i = 1:numel(x_cukf)
                t_cukf{i} = t_cukf{i} .* TU / 86400;
                x_cukf{i} = x_cukf{i} .* permute(repmat(U, 1, size(x_cukf{i}, 1), size(x_cukf{i}, 2)), [2 3 1]);  
            end

            count = 1;
            for i = 1:numel(x_cukf)

                p.color = C(i, :);

                x_cukf_i = x_cukf{i}; W_cukf_i = W_cukf{i};
                x_cukf_i0 = squeeze(x_cukf_i(:, 1, :)); 
                mu0 = reconstruct_cut(x_cukf_i0, W_cukf_i, 1, 6); 
                S0  = reconstruct_cut(x_cukf_i0, W_cukf_i, 2, 6);
                x_cukf_if = squeeze(x_cukf_i(:, end, :)); 
                muf = reconstruct_cut(x_cukf_if, W_cukf_i, 1, 6); 
                Sf  = reconstruct_cut(x_cukf_if, W_cukf_i, 2, 6);
                
                p.type = "line"; 

                nexttile(nm);
                scatter(x_cukf_i0(:, 1), x_cukf_i0(:, 2), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(1:2), S0(1:2, 1:2), "p", p);
                scatter(x_cukf_if(:, 1), x_cukf_if(:, 2), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(muf(1:2), Sf(1:2, 1:2), "p", p);   

                nexttile(numel(filters) + nm);
                scatter(x_cukf_i0(:, 3), x_cukf_i0(:, 4), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(3:4), S0(3:4, 3:4), "p", p);
                scatter(x_cukf_if(:, 3), x_cukf_if(:, 4), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(muf(3:4), Sf(3:4, 3:4), "p", p);       

                drawnow; 
            end  

        case "cukf_pme"
            %% CUKF_PME
            FILE = "../results/" + model + "/cukf_pme/" + model + "_cukf_pme_vars.mat"; load(FILE); 
            for i = 1:numel(x_cukf_pme)
                t_cukf_pme{i} = t_cukf_pme{i} .* TU / 86400;
                x_cukf_pme{i} = x_cukf_pme{i} .* permute(repmat(U, 1, size(x_cukf_pme{i}, 1), size(x_cukf_pme{i}, 2)), [2 3 1]);  
                X_cukf_pme{i} = X_cukf_pme{i} .* U'; 
            end

            count = 1;
            for i = 1:numel(x_cukf_pme)

                p.color = C(i, :);

                x_cukf_i = x_cukf_pme{i}; W_cukf_i = W_cukf_pme{i};
                X_cukf_i = X_cukf_pme{i}; P_cukf_i = P_cukf_pme{i};
                x_cukf_i0 = squeeze(x_cukf_i(:, 1, :)); 
                mu0 = reconstruct_cut(x_cukf_i0, W_cukf_i, 1, 6); 
                S0  = reconstruct_cut(x_cukf_i0, W_cukf_i, 2, 6);
                x_cukf_if = squeeze(x_cukf_i(:, end, :)); 
                muf = reconstruct_cut(x_cukf_if, W_cukf_i, 1, 6); 
                Sf  = reconstruct_cut(x_cukf_if, W_cukf_i, 2, 6);
                
                nexttile(nm);
                scatter(x_cukf_i0(:, 1), x_cukf_i0(:, 2), 10, p.color, "filled", "d");
                p.type = "line"; 
                plot_gaussian_ellipsoid(mu0(1:2), S0(1:2, 1:2), "p", p);
                scatter(x_cukf_if(:, 1), x_cukf_if(:, 2), 10, p.color, "filled", "d");
                p.type = "grid";
                plot_nongaussian_surface(X_cukf_i(:, 1:2), "P", P_cukf_i, "p", p);      

                nexttile(numel(filters) + nm);
                scatter(x_cukf_i0(:, 3), x_cukf_i0(:, 4), 10, p.color, "filled", "d");
                p.type = "line"; 
                plot_gaussian_ellipsoid(mu0(3:4), S0(3:4, 3:4), "p", p);
                scatter(x_cukf_if(:, 3), x_cukf_if(:, 4), 10, p.color, "filled", "d");
                p.type = "grid";
                plot_nongaussian_surface(X_cukf_i(:, 3:4), "P", P_cukf_i, "p", p);         

                drawnow; 
            end  
        case "gbees"
            %% GBEES
            FILE = "../results/" + model + "/gbees/" + model + "_gbees_vars.mat"; load(FILE); 
            for i = 1:numel(X_gbees)
                X_gbees{i} = X_gbees{i} .* U';
            end

            count = 1;
            for i = 1:2:numel(X_gbees)

                p.color = C(ceil(i / 2), :);
                
                nexttile(nm); 
                p.type = "grid";
                plot_nongaussian_surface(X_gbees{i}(:, 1:2), "P", P_gbees{i}, "p", p);
                plot_nongaussian_surface(X_gbees{i + 1}(:, 1:2), "P", P_gbees{i + 1}, "p", p);

                nexttile(numel(filters) + nm); 
                plot_nongaussian_surface(X_gbees{i}(:, 3:4), "P", P_gbees{i}, "p", p);
                plot_nongaussian_surface(X_gbees{i + 1}(:, 3:4), "P", P_gbees{i + 1}, "p", p);

                drawnow; 
            end              

        case "cukf_gbees"
            %% CUKF GBEES
            FILE = "../results/" + model + "/cukf_gbees/" + model + "_cukf_gbees_vars.mat"; load(FILE); 
            for i = 1:numel(x_cukf_gbees)
                t_cukf_gbees{i} = t_cukf_gbees{i} .* TU / 86400;
                x_cukf_gbees{i} = x_cukf_gbees{i} .* permute(repmat(U, 1, size(x_cukf_gbees{i}, 1), size(x_cukf_gbees{i}, 2)), [2 3 1]);  
            end

            for i = 1:numel(X_cukf_gbees) 
                X_cukf_gbees{i} = X_cukf_gbees{i} .* U';
            end

            count = 1;
            for i = 1:2:numel(X_cukf_gbees)

                p.color = C(ceil(i / 2), :);
                p.type = "grid"; 

                nexttile(nm); 
                plot_nongaussian_surface(X_cukf_gbees{i + 1}(:, 1:2), "P", P_cukf_gbees{i + 1}, "p", p);

                nexttile(numel(filters) + nm); 
                plot_nongaussian_surface(X_cukf_gbees{i + 1}(:, 3:4), "P", P_cukf_gbees{i + 1}, "p", p);

                drawnow; 
            end           

            count = 1;
            for i = 1:numel(x_cukf_gbees)

                p.color = C(i, :);

                x_cukf_i = x_cukf_gbees{i}; W_cukf_i = W_cukf_gbees{i};
                x_cukf_i0 = squeeze(x_cukf_i(:, 1, :)); 
                mu0 = reconstruct_cut(x_cukf_i0, W_cukf_i, 1, 6); 
                S0  = reconstruct_cut(x_cukf_i0, W_cukf_i, 2, 6);
                x_cukf_if = squeeze(x_cukf_i(:, end, :)); 
                muf = reconstruct_cut(x_cukf_if, W_cukf_i, 1, 6); 
                Sf  = reconstruct_cut(x_cukf_if, W_cukf_i, 2, 6);
                
                p.type = "line"; 

                nexttile(nm);
                scatter(x_cukf_i0(:, 1), x_cukf_i0(:, 2), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(1:2), S0(1:2, 1:2), "p", p);
                scatter(x_cukf_if(:, 1), x_cukf_if(:, 2), 10, p.color, "filled", "d");   
                plot_gaussian_ellipsoid(muf(1:2), Sf(1:2, 1:2), "p", p);

                nexttile(numel(filters) + nm);
                scatter(x_cukf_i0(:, 3), x_cukf_i0(:, 4), 10, p.color, "filled", "d");
                plot_gaussian_ellipsoid(mu0(3:4), S0(3:4, 3:4), "p", p);
                scatter(x_cukf_if(:, 3), x_cukf_if(:, 4), 10, p.color, "filled", "d"); 
                plot_gaussian_ellipsoid(muf(3:4), Sf(3:4, 3:4), "p", p);

                drawnow; 
            end  
    end
end

clear L; clear LH; 
LH(1) = plot(NaN, NaN, "o", "MarkerSize", 10, "MarkerFaceColor","c","MarkerEdgeColor","k");
L{1} = "Enceladus\,\,\,";
LH(2) = plot(NaN, NaN, "k-", "LineWidth",1);
L{2} = "Nominal\,\,\,";
LH(3) = scatter(NaN, NaN, 10, "k", "filled", "d");
L{3} = "CUT Sigma Point\,\,\,";
LH(4) = fill(nan, nan, nan, "FaceAlpha", 0.7, "FaceColor", C(1,:), "EdgeColor", "none");
L{4} = "$p(\mathbf{x},t^{(0)})\,\,\,$";
LH(5) = fill(nan, nan, nan, "FaceAlpha", 0.7, "FaceColor", C(2,:), "EdgeColor", "none");
L{5} = "$p(\mathbf{x},t^{(1)})\,\,\,$";
LH(6) = fill(nan, nan, nan, "FaceAlpha", 0.7, "FaceColor", C(3,:), "EdgeColor", "none");
L{6} = "$p(\mathbf{x},t^{(2)})\,\,\,$";
LH(7) = fill(nan, nan, nan, "FaceAlpha", 0.7, "FaceColor", C(4,:), "EdgeColor", "none");
L{7} = "$p(\mathbf{x},t^{(3)})\,\,\,$";
leg = legend(axs{1}, LH, L, "Orientation", "Horizontal", "FontSize", 24, "FontName", "times", "Interpreter", "latex");
leg.Layout.Tile = "south";
drawnow; 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, x] = parse_traj(filename)
    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
    data = data(2:end); 
    A = cellfun(@(x) str2num(x), data, 'UniformOutput', false); A = cell2mat(A); 
    t = A(:,1); x = A(:,2:end); 
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function axs = initialize_figures(filters, mu, LU, sec_r)
        
    f = figure(1); clf; hold on; f.Position = [250 100 numel(filters)*400 1200];
    tiledlayout(2, numel(filters), "TileSpacing", "compact", "Padding", "compact");

    for i = 1:numel(filters)
        if i == 1
            axs{1} = nexttile(i); hold on; box on; 
        else
            nexttile(i); hold on; box on; 
        end
        axis("equal")
        set(gca, "FontName", "Times", "FontSize", 22);
        xlabel("$x$ (km)", "Interpreter","latex", "FontSize", 24);
        ylabel("$y$ (km)", "Interpreter","latex", "FontSize", 24);
        xlim([237200, 240000]);
        ylim([-1750, 1750]);
        if strcmp(filters{i}, "cukf")
            title("CUKF6", "FontWeight", "bold", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
        elseif strcmp(filters{i}, "cukf_pme")
            title("CUKF6+PME ($M = 3$)", "FontWeight", "bold", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
        elseif strcmp(filters{i}, "cukf_gbees")
            title("CUKF6/GBEES", "FontWeight", "bold", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
        else
            title(upper(filters{i}), "FontWeight", "bold", "FontSize", 28, "FontName", "times", "Interpreter", "latex");
        end
        [C, ~, ~] = imread("./enceladus.jpeg");
        radius = size(C,1)/2; 
        center = [size(C,1)/2, size(C,2)/2]; 
        for j = 1:size(C,1)
            for k = 1:size(C,2)
                x = [j,k]; 
                if norm(x - center, 2) > radius
                    C(j,k,:) = uint8([255, 255, 255]); 
                end
            end
        end
        center = [(1 - mu) * LU, 0]; 
        x = [center(1) - (sec_r), center(1) + (sec_r)];
        y = [center(2) - (sec_r), center(2) + (sec_r)];
        image(x, y, C)
    end

    for i = 1:numel(filters)
        nexttile(numel(filters) + i); hold on; box on; 
        axis("equal")
        set(gca, "FontName", "Times", "FontSize", 22);
        xlabel("$\dot{x}$ (km/s)", "Interpreter","latex", "FontSize", 24);
        ylabel("$\dot{y}$ (km/s)", "Interpreter","latex", "FontSize", 24);
        xlim([-0.23, 0.23]);
        ylim([-0.29, 0.29]);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
