clear all; close all; clc; format long;

% parameters
Nmc = 10000; 
rng(0);

%% 2D
d = 2; 
mu = randn(d,1); A = randn(d); S = A*A';  
% mu = zeros(d, 1); S = eye(d); 
X = mvnrnd(mu, S, Nmc); 
[Z, W] = cut6(mu, S); 

mu_mc = mean(X)'; 
S_mc = cov(X, 1); 
Xc = X - mu_mc';             
Sinv_XcT = S_mc \ Xc';       
Q = Xc * Sinv_XcT;    
b1_mc = mean(Q(:).^3);
b2_mc = mean(diag(Q).^2);

mu_cut = reconstruct_cut(Z, W, 1, 6); 
S_cut  = reconstruct_cut(Z, W, 2, 6); 
Zc = Z - mu_cut';         
Sinv_ZcT = S_cut \ Zc';    
Q = Zc * Sinv_ZcT;         
q_diag = diag(Q);
W_outer = W * W';          
b1_cut = sum(sum(W_outer .* (Q.^3)));
b2_cut = sum(W .* (q_diag.^2));

fprintf("2D:\n")
fprintf("   Analytical B1: 0\n");
fprintf("           MC B1: %0.5f\n", b1_mc);
fprintf("          CUT B1: %0.5f\n\n", b1_cut);
fprintf("   Analytical B2: %d\n", d * (d + 2));
fprintf("           MC B2: %0.5f\n", b2_mc);
fprintf("          CUT B2: %0.5f\n\n", b2_cut);

%% 3D
d = 3; 
mu = randn(d,1); A = randn(d); S = A*A';  
% mu = zeros(d, 1); S = eye(d); 
X = mvnrnd(mu, S, Nmc); 
[Z, W] = cut6(mu, S); 

mu_mc = mean(X)'; 
S_mc = cov(X, 1); 
Xc = X - mu_mc';             
Sinv_XcT = S_mc \ Xc';       
Q = Xc * Sinv_XcT;     
b1_mc = mean(Q(:).^3);  
b2_mc = mean(diag(Q).^2);    

mu_cut = reconstruct_cut(Z, W, 1, 6); 
S_cut  = reconstruct_cut(Z, W, 2, 6); 
Zc = Z - mu_cut';          
Sinv_ZcT = S_cut \ Zc';    
Q = Zc * Sinv_ZcT;         
q_diag = diag(Q);
W_outer = W * W';         
b1_cut = sum(sum(W_outer .* (Q.^3)));
b2_cut = sum(W .* (q_diag.^2));

fprintf("3D:\n")
fprintf("   Analytical B1: 0\n");
fprintf("           MC B1: %0.5f\n", b1_mc);
fprintf("          CUT B1: %0.5f\n\n", b1_cut);
fprintf("   Analytical B2: %d\n", d * (d + 2));
fprintf("           MC B2: %0.5f\n", b2_mc);
fprintf("          CUT B2: %0.5f\n\n", b2_cut);


%% 4D
d = 4; 
mu = randn(d,1); A = randn(d); S = A*A';  
% mu = zeros(d, 1); S = eye(d); 
X = mvnrnd(mu, S, Nmc); 
[Z, W] = cut6(mu, S); 

mu_mc = mean(X)'; 
S_mc = cov(X, 1); 
Xc = X - mu_mc';             
Sinv_XcT = S_mc \ Xc';       
Q = Xc * Sinv_XcT;       
b1_mc = mean(Q(:).^3);   
b2_mc = mean(diag(Q).^2);    

mu_cut = reconstruct_cut(Z, W, 1, 6); 
S_cut  = reconstruct_cut(Z, W, 2, 6); 
Zc = Z - mu_cut';          
Sinv_ZcT = S_cut \ Zc';    
Q = Zc * Sinv_ZcT;        
q_diag = diag(Q);
W_outer = W * W';         
b1_cut = sum(sum(W_outer .* (Q.^3)));
b2_cut = sum(W .* (q_diag.^2));

fprintf("4D:\n")
fprintf("   Analytical B1: 0\n");
fprintf("           MC B1: %0.5f\n", b1_mc);
fprintf("          CUT B1: %0.5f\n\n", b1_cut);
fprintf("   Analytical B2: %d\n", d * (d + 2));
fprintf("           MC B2: %0.5f\n", b2_mc);
fprintf("          CUT B2: %0.5f\n\n", b2_cut);


%% 5D
d = 5; 
mu = randn(d,1); A = randn(d); S = A*A';  
% mu = zeros(d, 1); S = eye(d); 
X = mvnrnd(mu, S, Nmc); 
[Z, W] = cut6(mu, S); 

mu_mc = mean(X)'; 
S_mc = cov(X, 1); 
Xc = X - mu_mc';             
Sinv_XcT = S_mc \ Xc';       
Q = Xc * Sinv_XcT;       
b1_mc = mean(Q(:).^3);   
b2_mc = mean(diag(Q).^2);    

mu_cut = reconstruct_cut(Z, W, 1, 6); 
S_cut  = reconstruct_cut(Z, W, 2, 6); 
Zc = Z - mu_cut';          
Sinv_ZcT = S_cut \ Zc';    
Q = Zc * Sinv_ZcT;        
q_diag = diag(Q);
W_outer = W * W';         
b1_cut = sum(sum(W_outer .* (Q.^3)));
b2_cut = sum(W .* (q_diag.^2));

fprintf("5D:\n")
fprintf("   Analytical B1: 0\n");
fprintf("           MC B1: %0.5f\n", b1_mc);
fprintf("          CUT B1: %0.5f\n\n", b1_cut);
fprintf("   Analytical B2: %d\n", d * (d + 2));
fprintf("           MC B2: %0.5f\n", b2_mc);
fprintf("          CUT B2: %0.5f\n\n", b2_cut);

%% 6D
d = 6; 
mu = randn(d,1); A = randn(d); S = A*A';  
% mu = zeros(d, 1); S = eye(d); 
X = mvnrnd(mu, S, Nmc); 
[Z, W] = cut6(mu, S); 

mu_mc = mean(X)'; 
S_mc = cov(X, 1); 
Xc = X - mu_mc';             
Sinv_XcT = S_mc \ Xc';       
Q = Xc * Sinv_XcT;       
b1_mc = mean(Q(:).^3);   
b2_mc = mean(diag(Q).^2);    

mu_cut = reconstruct_cut(Z, W, 1, 6); 
S_cut  = reconstruct_cut(Z, W, 2, 6); 
Zc = Z - mu_cut';          
Sinv_ZcT = S_cut \ Zc';    
Q = Zc * Sinv_ZcT;        
q_diag = diag(Q);
W_outer = W * W';         
b1_cut = sum(sum(W_outer .* (Q.^3)));
b2_cut = sum(W .* (q_diag.^2));

fprintf("6D:\n")
fprintf("   Analytical B1: 0\n");
fprintf("           MC B1: %0.5f\n", b1_mc);
fprintf("          CUT B1: %0.5f\n\n", b1_cut);
fprintf("   Analytical B2: %d\n", d * (d + 2));
fprintf("           MC B2: %0.5f\n", b2_mc);
fprintf("          CUT B2: %0.5f\n\n", b2_cut);


