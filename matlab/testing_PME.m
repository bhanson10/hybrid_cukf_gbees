% testing_PME.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

mu = 1.901109735892602E-7;
LU = 238529;
TU = 18913;
U = [LU; LU; (LU / TU); (LU / TU)];

[X1, P1] = parse_pdf("../results/PCR3BP/testing/prior.txt");
[X2, P2] = parse_pdf("../results/PCR3BP/testing/posterior.txt");

X1_D = X1 .* U';
X2_D = X2 .* U';

p1.color = "r";
p2.color = "g";

figure; hold on; axis equal; 

plot_nongaussian_surface(X1_D(:, 1:2), "P", P1, "p", p1); 
plot_nongaussian_surface(X2_D(:, 1:2), "P", P2, "p", p2); 
xlim([2.37e5, 2.4e5]);
ylim([-1750, 1750]);

figure; hold on; axis equal; 
plot_nongaussian_surface(X1_D(:, 3:4), "P", P1, "p", p1); 
plot_nongaussian_surface(X2_D(:, 3:4), "P", P2, "p", p2); 
xlim([-0.2, 0.2]);
ylim([-0.4, 0.4]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [X, P] = parse_pdf(filename)

    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
    data = data(2:end); 
    pdf = cellfun(@(x) str2num(x), data, 'UniformOutput', false); pdf = cell2mat(pdf); 
    P = pdf(:, 1); X = pdf(:, 2:end); 
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%