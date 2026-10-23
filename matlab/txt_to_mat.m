% txt_to_mat.m
% Copyright 2026 by Benjamin Hanson, published under BSD 3-Clause License.

clear all; close all; clc; 

filter = {"cukfpf"}; 
orbit = "NRHO";
len = "30d"; 
n_str = "5000"; 
HZ_str = "1-047"; 
dim = 6; 

for nm = 1:numel(filter)
    PATH = "../results/" + orbit + "/" + filter{nm} + "/" + len;

    if strcmp(filter{nm}, "traj")
        
        FILE = PATH + "/" + orbit + "_rv.txt",
        [t_traj, x_traj] = parse_traj(FILE); 

        save(PATH + "/" + orbit + "_traj_vars", "t_traj", "x_traj", "-v7.3");

    elseif strcmp(filter{nm}, "pf")

        PATH = PATH + "_" + n_str; 
        
        % particles
        FILE_LIST = dir(fullfile(PATH + "/particles", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/particles/" + orbit + "_rv_pf_" + num2str(i) + ".txt",
            [t_pf{i + 1}, W_pf{i + 1}, x_pf{i + 1}] = parse_pf(FILE, dim);
        end

        save(PATH + "/" + orbit + "_pf_vars", "t_pf", "W_pf", "x_pf", "-v7.3");

    elseif strcmp(filter{nm}, "ukf")
        FILE_LIST = dir(fullfile(PATH + "/sigmas", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/sigmas/" + orbit + "_rv_ukf_" + num2str(i) + ".txt",
            [t_ukf{i + 1}, W_ukf{i + 1}, x_ukf{i + 1}] = parse_ukf(FILE);
        end

        save(PATH + "/" + orbit + "_ukf_vars", "t_ukf", "W_ukf", "x_ukf", "-v7.3");

    elseif strcmp(filter{nm}, "cukf")

        FILE_LIST = dir(fullfile(PATH + "/sigmas", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/sigmas/" + orbit + "_rv_cukf_" + num2str(i) + ".txt",
            [t_cukf{i + 1}, W_cukf{i + 1}, x_cukf{i + 1}, ~] = parse_cukf(FILE, dim, 0);
        end

        save(PATH + "/" + orbit + "_cukf_vars", "t_cukf", "W_cukf", "x_cukf", "-v7.3");
    
     elseif strcmp(filter{nm}, "cukf_pme")
        FILE_LIST = dir(fullfile(PATH + "/sigmas", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/sigmas/" + orbit + "_rv_cukf_pme_" + num2str(i) + ".txt",
            [t_cukf_pme{i + 1}, W_cukf_pme{i + 1}, x_cukf_pme{i + 1}, ~] = parse_cukf(FILE, dim, 0);
        end
        
        FILE_LIST = dir(fullfile(PATH + "/pdfs", '*.bin'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/pdfs/" + orbit + "_rv_pdf_" + num2str(i) + ".bin",
            [P_cukf_pme{i + 1}, X_cukf_pme{i + 1}] = parse_pdf(FILE, 4, 1, "double"); 
            % [X_cukf_pme{i + 1}, P_cukf_pme{i + 1}] = parse_pdf(FILE);
        end

        save(PATH + "/" + orbit + "_cukf_pme_vars", "t_cukf_pme", "W_cukf_pme", "x_cukf_pme", "X_cukf_pme", "P_cukf_pme", "-v7.3");

    elseif strcmp(filter{nm}, "cukfpf")
        
        PATH = PATH + "_" + n_str + "_HZ_" + HZ_str; 
        
        % sigmas
        FILE_LIST = dir(fullfile(PATH + "/sigmas", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/sigmas/" + orbit + "_rv_cukf_" + num2str(i) + ".txt",
            [t_cukfpf_cukf{i + 1}, W_cukfpf_cukf{i + 1}, x_cukfpf_cukf{i + 1}, ~] = parse_cukf(FILE, dim, 0);
        end

        % particles
        FILE_LIST = dir(fullfile(PATH + "/particles", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/particles/" + orbit + "_rv_pf_" + num2str(i) + ".txt",
            [t_cukfpf_pf{i + 1}, W_cukfpf_pf{i + 1}, x_cukfpf_pf{i + 1}] = parse_pf(FILE, dim);
        end
        
        [t_b, b, t_HZ, HZ] = parse_b_HZ("../c/examples/" + orbit + "/outputs/" + orbit + "_" + filter{nm} + "_output.log"); 

        save(PATH + "/" + orbit + "_cukfpf_vars", "t_cukfpf_cukf", "W_cukfpf_cukf", "x_cukfpf_cukf", ...
                                          "t_cukfpf_pf", "W_cukfpf_pf", "x_cukfpf_pf", ...
                                          "t_b", "b", "t_HZ", "HZ", "-v7.3");

    elseif strcmp(filter{nm}, "gbees")
        DIR_LIST = dir(fullfile(PATH + "/pdfs"));  % List only .txt files
        N_dirs = sum([DIR_LIST.isdir] & ~ismember({DIR_LIST.name}, {'.', '..'}));
        count = 1; 
        for i = 0:N_dirs - 1
            PATH_i = PATH + "/pdfs/P" + num2str(i); 
            FILE_LIST = dir(fullfile(PATH_i, '*.bin'));  % List only .bin files
            for j = 0:numel(FILE_LIST) - 1
                FILE = PATH_i + "/P" + num2str(i) + "_pdf_" + num2str(j) + ".bin",
                [P_gbees{count}, X_gbees{count}] = parse_pdf(FILE, 4, 1, "float"); 
                count = count + 1; 
            end 
        end
        save(PATH + "/" + orbit + "_gbees_vars", "P_gbees", "X_gbees", "-v7.3");

    elseif strcmp(filter{nm}, "cukf_gbees")

        FILE_LIST = dir(fullfile(PATH + "/sigmas", '*.txt'));  % List only .txt files
        N_sets = numel(FILE_LIST);
        for i = 0:N_sets - 1 
            FILE = PATH + "/sigmas/" + orbit + "_rv_cukf_" + num2str(i) + ".txt",
            [t_cukf_gbees{i + 1}, W_cukf_gbees{i + 1}, x_cukf_gbees{i + 1}, b_cukf_gbees{i + 1}] = parse_cukf(FILE, dim, 1);
        end

        DIR_LIST = dir(fullfile(PATH + "/pdfs"));  % List only .txt files
        N_dirs = sum([DIR_LIST.isdir] & ~ismember({DIR_LIST.name}, {'.', '..'}));
        count = 1; 
        for i = 0:N_dirs - 1
            PATH_i = PATH + "/pdfs/P" + num2str(i); 
            FILE_LIST = dir(fullfile(PATH_i, '*.bin'));  % List only .bin files
            for j = 0:numel(FILE_LIST) - 1
                FILE = PATH_i + "/P" + num2str(i) + "_pdf_" + num2str(j) + ".bin",
                [P_cukf_gbees{count}, X_cukf_gbees{count}] = parse_pdf(FILE, 4, 1, "float"); 
                count = count + 1; 
            end 
        end
        save(PATH + "/" + orbit + "_cukf_gbees_vars", "t_cukf_gbees", "W_cukf_gbees", ...
            "x_cukf_gbees", "b_cukf_gbees", "P_cukf_gbees", "X_cukf_gbees", "-v7.3");
    end
end

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
function [t, W, x] = parse_pf(filename, dim)

    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
    data = data(2:end); 
    A = cellfun(@(x) str2num(x), data, 'UniformOutput', false); A = cell2mat(A); 
    t = A(:,1); 
    A = reshape(A(:,2:end).', dim+1,[],size(A,1));
    W = squeeze(A(1,:,1)).';
    x = permute(A(2:dim+1,:,:),[2 3 1]);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, W, x] = parse_ukf(filename, dim)

    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
    data = data(2:end); 
    A = cellfun(@(x) str2num(x), data, 'UniformOutput', false); A = cell2mat(A); 
    t = A(:,1); 
    A = reshape(A(:,2:end).', dim+2,[],size(A,1));
    W = squeeze(A(1:2,:,1)).';
    x = permute(A(3:dim+2,:,:),[2 3 1]);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, W, x, b] = parse_cukf(filename, dim, flag)

    fileID = fopen(filename, 'r');
    data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
    data = data(2:end); 
    A = cellfun(@(x) str2num(x), data, 'UniformOutput', false); A = cell2mat(A); 
    t = A(:,1); 
    if flag
        b = A(:, end - 1:end);
        A = reshape(A(:,2:end - 2).', dim + 1,[],size(A,1));
    else
        b = NaN; 
        A = reshape(A(:,2:end).', dim + 1,[],size(A,1));
    end
    W = squeeze(A(1,:,1)).';
    x = permute(A(2:dim + 1,:,:),[2 3 1]);
    fclose(fileID);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P, x] = parse_pdf(filename, dim, binary, precision)
    if nargin < 3
        precision = "single"; 
    end
    if nargin < 3
        binary = 0; 
    end
    if(binary)
        fid = fopen(filename, 'rb');
        assert(fid ~= -1, 'Could not open file.');
        t = fread(fid, 1, "double");
        n = fread(fid, 1, 'uint32');
        raw = fread(fid, [dim+1, n], precision);
        P = raw(1,:)';
        x = raw(2:end,:)';
        fclose(fid);
    else
        fileID = fopen(filename, 'r');
        data = textscan(fileID, '%s', 'Delimiter', '\n'); data = data{1}; 
        data = data(2:end); 
        pdf = cellfun(@(x) str2num(x), data, 'UniformOutput', false); pdf = cell2mat(pdf); 
        P = pdf(:,1); x = pdf(:,2:2+dim-1);
        fclose(fileID);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t_b, b, t_HZ, HZ] = parse_b_HZ(filename)
% PARSE_LOG_FILE Parses log files containing time, b1, b2, PF mode, and HZ data.
%
% Outputs:
%   t_b  : N1x1 cell array where each cell contains M1x1 time values for b1 & b2
%   b    : N1x1 cell array where each cell contains M1x2 [b1, b2] values
%   t_HZ : N2x1 cell array where each cell contains M2x1 time values for valid HZ
%   HZ   : N2x1 cell array where each cell contains M2x1 valid HZ values (excluding NaNs)

    % Open the file for reading
    fileID = fopen(filename, 'r');
    if fileID == -1
        error('Cannot open file: %s', filename);
    end
    
    % Initialize temporary holding variables
    curr_t_b = [];
    curr_b   = [];
    
    curr_t_HZ = [];
    curr_HZ   = [];
    
    t_b  = {};
    b    = {};
    t_HZ = {};
    HZ   = {};
    
    % Read line by line
    while ~feof(fileID)
        line = strtrim(fgetl(fileID));
        
        % Skip empty lines
        if isempty(line)
            continue;
        end
        
        % 1. Check for b1 and b2 lines using Regular Expression
        % Matching pattern: t: <float>, b1: <float>, b2: <float>
        b_tokens = regexp(line, 't:\s*([\d\.]+),\s*b1:\s*([\d\.\-\+eE]+),\s*b2:\s*([\d\.\-\+eE]+)', 'tokens');
        
        if ~isempty(b_tokens)
            vals = str2double(b_tokens{1});
            curr_t_b = [curr_t_b; vals(1)];
            curr_b   = [curr_b; vals(2), vals(3)];
            continue;
        end
        
        % 2. Check for HZ lines
        % Matching pattern: t: <float>, HZ: <float_or_NaN>
        hz_tokens = regexp(line, 't:\s*([\d\.]+),\s*HZ:\s*([\d\.\-\+eE]+|NaN|nan)', 'tokens');
        
        if ~isempty(hz_tokens)
            t_val  = str2double(hz_tokens{1}{1});
            hz_val = str2double(hz_tokens{1}{2});
            
            % Save valid (non-NaN) HZ values
            if ~isnan(hz_val)
                curr_t_HZ = [curr_t_HZ; t_val];
                curr_HZ   = [curr_HZ; hz_val];
            end
            continue;
        end
        
        % 3. Check for section/run separators if multiple runs exist in one file
        % (e.g., lines starting with '---', 'Run', 'Header', etc.)
        if startsWith(line, '---') || startsWith(line, 'Run')
            if ~isempty(curr_t_b)
                t_b{end+1, 1} = curr_t_b;
                b{end+1, 1}   = curr_b;
                curr_t_b = []; curr_b = [];
            end
            if ~isempty(curr_t_HZ)
                t_HZ{end+1, 1} = curr_t_HZ;
                HZ{end+1, 1}   = curr_HZ;
                curr_t_HZ = []; curr_HZ = [];
            end
        end
    end
    
    fclose(fileID);
    
    % Push final accumulated data into cell arrays
    if ~isempty(curr_t_b)
        t_b{end+1, 1} = curr_t_b;
        b{end+1, 1}   = curr_b;
    end
    
    if ~isempty(curr_t_HZ)
        t_HZ{end+1, 1} = curr_t_HZ;
        HZ{end+1, 1}   = curr_HZ;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [t_b, b, t_HZ, HZ] = parse_b_HZ(filename)
% 
%     fileID = fopen(filename, 'r');
%     t_b = {}; b = {}; t_HZ = {}; HZ = {};
%     current_type = '';  
%     temp_t = []; temp_val = [];
% 
%     % Read file line by line
%     while ~feof(fileID)
%         line = fgetl(fileID);
% 
%         % Skip empty lines
%         if isempty(line) || all(isspace(line))
%             continue;
%         end
% 
%         % Parse the line
%         tokens = regexp(line, 't:\s*([\d.]+),\s*(b|HZ):\s*([\d.]+)', 'tokens');
% 
%         if isempty(tokens)
%             continue; % skip lines that don't match
%         end
% 
%         tokens = tokens{1};
%         t = str2double(tokens{1});
%         val_type = tokens{2};
%         val = str2double(tokens{3});
% 
%         % Check if type has changed
%         if isempty(current_type)
%             current_type = val_type;
%         elseif ~strcmp(val_type, current_type)
%             % Save the previous block
%             if strcmp(current_type, 'b')
%                 t_b{end+1} = temp_t;
%                 b{end+1} = temp_val;
%             else
%                 t_HZ{end+1} = temp_t;
%                 HZ{end+1} = temp_val;
%             end
%             % Reset temporary storage
%             temp_t = [];
%             temp_val = [];
%             current_type = val_type;
%         end
% 
%         % Append to current block
%         temp_t(end+1) = t;
%         temp_val(end+1) = val;
%     end
% 
%     % Save the last block
%     if ~isempty(temp_t)
%         if strcmp(current_type, 'b')
%             t_b{end+1} = temp_t;
%             b{end+1} = temp_val;
%         else
%             t_HZ{end+1} = temp_t;
%             HZ{end+1} = temp_val;
%         end
%     end
% 
%     fclose(fileID);
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%