%% Run PopPK (RIF)
% 
% DESCRIPTION:
% This script uses the equations and parameters found in Ramachandran &
% Gadgil, 2023, in order to compare an oral dose of rifampin to an
% inhaled/lung dose of rifampin. 
% 
% The dosing regimen of each may be altered with oral/lung_dose_RIF and
% oral/lung_dose_freq_RIF. The compartments which will be plotted and for
% which PK metrics (AUC_24, C_avg, and C_max) will be calculated may be
% altered with relevant_compts_RIF.
%
% PLOTS GENERATED:
% - Concentration-time courses of all patients in relevant compts.
% - 10th, 50th, and 90th percentiles of patient timecourses in relevant
%   compts.
% - MIC distribution and AUC/MIC probability of target attainment (PTA) 
%   comparison
% - MIC distribution and C_max/MIC PTA comparison
%
% FILES GENERATED:
% - popPK_analysis.xlsx: compares PK metrics across patient populations for
%   oral vs. lung dosing
% - popPK_CFRs.xlsx: contains estimated cumulative fractions of response 
%   (CFRs) for AUC/MIC and C_max/MIC PTA target values
%
% AUTHORS:
% Noah Strawhacker (nstrawha@purdue.edu)
% Alexis Hoerter

clc;
clearvars;

addpath("Oral_Dose_ODEs/");
addpath("Lung_Dose_ODEs/");
addpath("Methods/");


%% Set parameters

% modeling parameters
n_pts_RIF = 50;
n_days_RIF = 4;
days_to_plot_RIF = 1;
relevant_compts_RIF = {"Lung", "Liver"};

oral_dose_RIF = 600;    % mg
lung_dose_RIF = 600;    % mg
oral_dose_freq_RIF = 1; % doses/day
lung_dose_freq_RIF = 1; % doses/day

tstep_RIF = 0.01;

% RIF-specific parameters
ka_oral_RIF = 1.08;     % absorption rate [1/h]
kdiss_lung_RIF = 50;    % dissolution rate [1/h] from Himstedt et al.
kF_RIF = 0.252;         % gut transit rate
kr_RIF = 0.17;          % gut reabsorption rate [1/h]

CL_RIF = 7.86;          % systemic clearance [L/h]
fR_RIF = 0.1830;        % fractional renal clearance

br_frac_RIF = 9/49;     % proportion of RIF absorbed by bronchi from Himstedt et al.
effRB_RIF = 2.56;       % bronchi efflux ratio from Himstedt et al.
effRA_RIF = 3.67;       % alveolar efflux ratio from Himstedt et al.

% model compts      name                index in concentration array
compt_list_RIF =   ["Plasma";           % 1
                    "Arterial Blood";   % 2
                    "Lung";             % 3
                    "Pleura";           % 4
                    "Brain";            % 5
                    "Adipose Tissue";   % 6
                    "Heart";            % 7
                    "Muscle";           % 8
                    "Skin";             % 9
                    "Other Tissue";     % 10
                    "Bone";             % 11
                    "Spleen";           % 12
                    "Kidney";           % 13
                    "Gut";              % 14
                    "Liver";            % 15
                    "Lymph Node";       % 16
                    "Gut Lumen";        % 17
                    "ELFb";             % 18
                    "ELFa";             % 19
                    "Absorption"];      % 20;

ncompts_total_RIF = length(compt_list_RIF);
toxic_compts_RIF = ["Liver", "Kidney"];

% initialize cell arrays to store params and output
param_store_RIF = cell(1, n_pts_RIF);
Cs_oral_store_RIF = cell(1, n_pts_RIF);
Cs_lung_store_RIF = cell(1, n_pts_RIF);


%% Iterate through patient equations

% sample physiological parameters
[vol_PDs_RIF, vol_frac_PDs_RIF, ...
    flow_PDs_RIF, flow_frac_PDs_RIF] = getParamPDs(0.2, 0.3); % CVs from Lyons et al.


ts_RIF = 0:tstep_RIF:(24 * n_days_RIF - tstep_RIF);

parfor pt_idx = 1:n_pts_RIF
    
    [phys_RIF, pt_RIF]  = loadPhysParams("RIF", vol_PDs_RIF, vol_frac_PDs_RIF, ...
                                            flow_PDs_RIF, flow_frac_PDs_RIF, ...
                                            0); % unfixed weight

    % package params                % index
    params_RIF =   {oral_dose_RIF;  % 1
                    lung_dose_RIF;  % 2
                    ka_oral_RIF;    % 3
                    CL_RIF;         % 4
                    fR_RIF;         % 5
                    kr_RIF;         % 6
                    pt_RIF;         % 7
                    phys_RIF;       % 8
                    kF_RIF;         % 9
                    kdiss_lung_RIF; % 10
                    effRB_RIF;      % 11
                    effRA_RIF;      % 12
                    br_frac_RIF;    % 13
                    tstep_RIF};     % 14

    % solve ODEs
    [C_oraldose_RIF_pt, ...
        C_lungdose_RIF_pt] = solveODEs("RIF", params_RIF, ncompts_total_RIF, n_days_RIF, ...
                                        oral_dose_freq_RIF, lung_dose_freq_RIF);
    % store results
    param_store_RIF{pt_idx} = params_RIF;
    Cs_oral_store_RIF{pt_idx} = C_oraldose_RIF_pt;
    Cs_lung_store_RIF{pt_idx} = C_lungdose_RIF_pt;

    % track progress
    disp(append("pt. #", num2str(pt_idx), " calculated..."))

end


%% Plot resulting timecourses

% iterate through compartments
for compt_idx = 1:length(relevant_compts_RIF)
    compt = relevant_compts_RIF{compt_idx};
    [idx_to_plot, ~] = find(string(compt_list_RIF) == compt);

    % pull all concentration timecourses for the current compartment
    current_Cs_orals = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_oral_store_RIF, "UniformOutput", false));
    current_Cs_lungs = cell2mat(cellfun(@(x) x(:, idx_to_plot), Cs_lung_store_RIF, "UniformOutput", false));

    plotTimecourses("RIF", oral_dose_RIF, lung_dose_RIF, ...
                    compt, n_days_RIF, days_to_plot_RIF, ...
                    oral_dose_freq_RIF, lung_dose_freq_RIF, ...
                    ts_RIF, current_Cs_orals, current_Cs_lungs);

end


%% Analysis

% compare PK metrics
[AUCs_oral_RIF, ...
    AUCs_lung_RIF, ...
    Cmaxs_oral_RIF, ...
    Cmaxs_lung_RIF] = trackPKMetrics("RIF", compt_list_RIF, relevant_compts_RIF, toxic_compts_RIF, ...
                                        oral_dose_RIF, oral_dose_freq_RIF, ...
                                        lung_dose_RIF, lung_dose_freq_RIF, ...
                                        n_days_RIF, tstep_RIF, ts_RIF, n_pts_RIF, ...
                                        Cs_oral_store_RIF, Cs_lung_store_RIF);


%% Calculate PTAs for oral and lung dosing

% MIC dist formatted as [conc., # isolates]
MICs_TB_RIF =  {[0.031, 28];
                [0.062, 47];
                [0.120, 51];
                [0.250, 76];
                [0.500, 81];
                [1.000, 68];
                [2.000, 0];
                [4.000, 5];
                [8.000, 0];
                [16.00, 0]};

RIF_AUC_target = 271; % AUC/MIC, from Jayaram et al.
RIF_Cmax_target = 175; % Cmax/MIC, from Gumbo et al.
% RIF_AUC_target = 1360; % AUC/MIC, from Gumbo (more conservative, UNUSED)

% record nontoxic compartments to calculate PTA for
nontoxic_compts_RIF = [];
for compt_idx = 1:length(relevant_compts_RIF)
    current_compt = relevant_compts_RIF{compt_idx};

    if ~ismember(current_compt, toxic_compts_RIF)
        nontoxic_compts_RIF = [nontoxic_compts_RIF, current_compt];
    end

end

% calculate and compare PTAs
plotPTAs("RIF", nontoxic_compts_RIF, ...
            AUCs_oral_RIF, AUCs_lung_RIF, ...
            Cmaxs_oral_RIF, Cmaxs_lung_RIF, ...
            MICs_TB_RIF, ...
            RIF_AUC_target, RIF_Cmax_target, ...
            oral_dose_RIF, oral_dose_freq_RIF, ...
            lung_dose_RIF, lung_dose_freq_RIF, ...
            n_pts_RIF, n_days_RIF);