%% Optimize Regimen (RIF)
% 
% DESCRIPTION:
%   Adaptation of run_popPK_RIF intended to compare dosing regimens. 
%   100-1000 mg amounts of inhaled rifampin are compared in 50 mg 
%   increments, with frequencies of 1-4 x/day, to the standard oral dose of
%   600 mg 1x/day. Measures resulting AUC_24 and C_max for each 
%   regimen.
%
% PLOTS GENERATED:
% - Level curves of all PK metrics
% 
% FILES GENERATED:
% - [PKMetric]_comparisons.xlsx: records data about all average PK metrics
%   for each dosing method and regimen
%
% AUTHORS:
% Noah Strawhacker (nstrawha@purdue.edu)

clc;
clearvars;

addpath("Oral_Dose_ODEs/");
addpath("Lung_Dose_ODEs/");
addpath("Methods/");


%% Set parameters

% Modeling parameters
n_days_RIF = 1;
days_to_plot_RIF = 1;
relevant_compts_RIF = {"Lung", "Liver"};
tstep_RIF = 0.01;

% Dosing regimen info
oral_dose_RIF = 600;        % mg
oral_dose_freq_RIF = 1;     % doses/day

lung_dose_min_RIF = 500;    % mg
lung_dose_max_RIF = 700;    % mg
lung_dose_inc_RIF = 50;     % mg
lung_dose_freq_min_RIF = 1; % doses/day
lung_dose_freq_max_RIF = 3; % doses/day

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

% sample physiological parameters
[vol_PDs_RIF, vol_frac_PDs_RIF, ...
    flow_PDs_RIF, flow_frac_PDs_RIF] = getParamPDs(0.2, 0.3); % CVs from Lyons et al.

[phys_RIF, pt_RIF]  = loadPhysParams("RIF", vol_PDs_RIF, vol_frac_PDs_RIF, ...
                                        flow_PDs_RIF, flow_frac_PDs_RIF, ...
                                        0); % unfixed weight


%% Solve oral equations

lung_dose_RIF = 1; % dummy
lung_dose_freq_RIF = 1; % dummy

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

ts_RIF = 0:tstep_RIF:(24 * n_days_RIF - tstep_RIF);

% solve ODEs
[C_oraldose_RIF, ~] = solveODEs("RIF", params_RIF, ncompts_total_RIF, n_days_RIF, ...
                                oral_dose_freq_RIF, lung_dose_freq_RIF);

disp("oral pt. calculated")

%% Iterate through lung regimens

dose_freqs_RIF = lung_dose_freq_min_RIF:lung_dose_freq_max_RIF;
lung_doses_RIF = lung_dose_min_RIF:lung_dose_inc_RIF:lung_dose_max_RIF;
regimen_store_RIF = cell(1, length(lung_doses_RIF) * length(dose_freqs_RIF));
storage_idx = 1;

for lung_dose_freq_RIF = dose_freqs_RIF
    for lung_dose_RIF = lung_doses_RIF
    
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
        [~, C_lungdose_RIF_pt] = solveODEs("RIF", params_RIF, ncompts_total_RIF, n_days_RIF, ...
                                            oral_dose_freq_RIF, lung_dose_freq_RIF);
        % store results
        regimen_store_RIF{storage_idx} = C_lungdose_RIF_pt;
        storage_idx = storage_idx + 1;
    
        % track progress
        disp(append("freq: ", num2str(lung_dose_freq_RIF), " x/day; dose: ", num2str(lung_dose_RIF), " mg"));
    
    end
end


%% Calculate and store metrics

% set up storage
AUCs_store_RIF  = cell(1, length(relevant_compts_RIF));
Cmaxs_store_RIF = cell(1, length(relevant_compts_RIF));

last_day_start = length(ts_RIF) - 24 / tstep_RIF + 1;

for compt_idx = 1:length(relevant_compts_RIF)
    current_compt = relevant_compts_RIF{compt_idx};
    [idx_to_calc, ~] = find(string(compt_list_RIF) == current_compt);

    % pull timecourses for the current compt
    current_Cs_oral = C_oraldose_RIF(:, idx_to_calc);
    current_Cs_lung = cell2mat(cellfun(@(x) x(:, idx_to_calc), regimen_store_RIF, "UniformOutput", false));

    % calculate metrics
    oral_AUC  = trapz(ts_RIF(last_day_start:end), current_Cs_oral(last_day_start:end, :));
    lung_AUCs = trapz(ts_RIF(last_day_start:end), current_Cs_lung(last_day_start:end, :));

    oral_Cmax  = max(current_Cs_oral(last_day_start:end, :));
    lung_Cmaxs = max(current_Cs_lung(last_day_start:end, :));

    % calculate percent differences
    AUC_percs_higher  = (lung_AUCs - oral_AUC) ./ oral_AUC .* 100;
    Cmax_percs_higher = (lung_Cmaxs - oral_Cmax) ./ oral_Cmax .* 100;

    % store results
    AUCs_store_RIF{compt_idx}  = AUC_percs_higher;
    Cmaxs_store_RIF{compt_idx} = Cmax_percs_higher;

end


%% Plot level curves

for compt_idx = 1:length(relevant_compts_RIF)
    current_AUC_percs = AUCs_store_RIF{compt_idx};
    current_AUC_percs = reshape(current_AUC_percs, length(lung_doses_RIF), length(dose_freqs_RIF));

    current_Cmax_percs = Cmaxs_store_RIF{compt_idx};
    current_Cmax_percs = reshape(current_Cmax_percs, length(lung_doses_RIF), length(dose_freqs_RIF));

    % create plot
    fig = figure();
    tlayout = tiledlayout(1, 2);
    [freqplot, doseplot] = meshgrid(dose_freqs_RIF, lung_doses_RIF);
    contour_levels_RIF = -1000:50:1000;

    % AUC
    nexttile
    [C_temp, h_temp] = contour(freqplot, doseplot, current_AUC_percs, contour_levels_RIF);
    clabel(C_temp, h_temp)
    title("AUC")
    set(gca, "xtick", dose_freqs_RIF)
    set(gca, "ytick", lung_doses_RIF)

    % Cmax
    nexttile
    [C_temp, h_temp] = contour(freqplot, doseplot, current_Cmax_percs, contour_levels_RIF);
    clabel(C_temp, h_temp)
    title("C_{max}")
    set(gca, "xtick", dose_freqs_RIF)
    set(gca, "ytick", lung_doses_RIF)

    xlabel(tlayout, "Dose Frequency (doses/day)");
    ylabel(tlayout, "Dose Amount (mg)");
    sgtitle(append(relevant_compts_RIF{compt_idx}, ...
                    " Compartment; % Increase from Oral Dose PK Metrics to Lung Dose PK Metrics"), ...
            "FontSize", 15);

    saveas(fig, append("Outputs/RIF/Figures/", relevant_compts_RIF{compt_idx}, "_regimen_comparison_RIF.png"));

end


%% TODO: Write output file

metrics_RIF = ["AUC", "Cmax"];
labels_RIF =   ["AUC_24; % Higher than Oral Dose", ...
                "Cmax; % Higher than Oral Dose"];

% reshape all comparison arrays
for compt_idx = 1:length(relevant_compts_RIF)
    AUCs_store_RIF{compt_idx}  = reshape(AUCs_store_RIF{compt_idx}, length(lung_doses_RIF), length(dose_freqs_RIF));
    Cmaxs_store_RIF{compt_idx} = reshape(Cmaxs_store_RIF{compt_idx}, length(lung_doses_RIF), length(dose_freqs_RIF));
end

all_metrics_store_RIF = {AUCs_store_RIF, Cmaxs_store_RIF};

% iterate through metric types
for metric_idx = 1:length(all_metrics_store_RIF)
    current_metrics = all_metrics_store_RIF{metric_idx};

    % iterate through compartments
    for compt_idx = 1:length(relevant_compts_RIF)
        metrics_for_current_compt = current_metrics{compt_idx};

        current_cell = cell(length(lung_doses_RIF) + 2, length(dose_freqs_RIF) + 2);
        current_cell{1, 1} = labels_RIF(metric_idx);
        current_cell{3, 1} = "Dose amt. (mg)";
        current_cell{1, 3} = "Dose freq. (doses/day)";
    
        % iterate through frequencies
        for freq_idx = 1:length(dose_freqs_RIF)
    
            % iterate through dose amts
            for dose_idx = 1:length(lung_doses_RIF)
                current_cell{dose_idx + 2, freq_idx + 2} = round(metrics_for_current_compt(dose_idx, freq_idx), 2);
    
            end
        end
    
        % add labels
        current_cell(3:end, 2) = num2cell(lung_doses_RIF(:));
        current_cell(2, 3:end) = num2cell(dose_freqs_RIF(:));

        % write output
        writecell(current_cell, ...
            append("Outputs/RIF/Tables/", metrics_RIF(metric_idx), "_comparisons.xlsx"), ...
            "Sheet", relevant_compts_RIF{compt_idx});

    end
end
