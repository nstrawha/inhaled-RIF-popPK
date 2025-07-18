function [vol_PDs, vol_frac_PDs, flow_PDs, flow_frac_PDs] = getParamPDs(v_cv, q_cv)
% GETPARAMPDS - Returns the probability distributions of physiological
% parameters given two CVs, one for volume parameters and one for flow
% parameters.
%
% DESCRIPTION:
%   This function sets up a normal probability distribution for each
%   physiological parameter based on its mean (taken from Ramachandran &
%   Gadgil, 2023) and CV (taken from Lyons et al., 2013). All distributions
%   are truncated at +/- 3 SDs to ensure biological feasibility.
%
% INPUTS:
% - v_cv (double): The coefficient of variation for all physiological 
%   volume parameters
% - q_cv (double): The coefficient of variation for all physiological 
%   flow parameters
%
% OUTPUTS:
% - vol_PDs (cell array): Contains all probability distributions for raw
%   volume parameters not based on body weight
% - vol_frac_PDs (cell array): Contains all probability distributions for
%   volume parameters that are calculated as fractions of body weight 
% - flow_PDs (cell array): Contains all probability distributions for raw
%   flow parameters not based on cardiac output
% - flow_PDs_fracs (cell array): Contains all probability distributions for
%   flow parameters that are calculated as fractions of cardiac output


%% Volumes

% raw volumes
v_v_mean  = 3.6; % venous
v_a_mean  = 1.8; % arterial
v_ln_mean = 0.274;

% vol fracs
v_lu_mean_frac      = 0.0076;
v_brain_mean_frac   = 0.02;
v_heart_mean_frac   = 0.0047;
v_adipose_mean_frac = 0.2142 / 0.916;
v_muscle_mean_frac  = 0.4;
v_skin_mean_frac    = 0.0371;
v_kidney_mean_frac  = 0.0044;
v_bone_mean_frac    = 0.1429 / 1.92;
v_spleen_mean_frac  = 0.0026;
v_gut_mean_frac     = 0.0171;
v_liver_mean_frac   = 0.0257;
v_others_mean_frac  = 0.04264;
v_pl_mean_frac      = 0.3 / 1000; % mL/kg to L

v_means =  [v_v_mean, v_a_mean, v_ln_mean];

v_means_fracs = [v_lu_mean_frac, v_brain_mean_frac, v_heart_mean_frac, v_adipose_mean_frac, ...
                v_muscle_mean_frac, v_skin_mean_frac, v_kidney_mean_frac, ...
                v_bone_mean_frac, v_spleen_mean_frac, v_gut_mean_frac, v_liver_mean_frac, ...
                v_others_mean_frac, v_pl_mean_frac];

v_sds = v_means .* v_cv;
v_sds_fracs = v_means_fracs .* v_cv;

% create cell array of truncated probability dists for each vol param
vol_PDs = cell(1, length(v_means));
for param_idx = 1:length(vol_PDs)
    current_mean = v_means(param_idx);
    current_sd = v_sds(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    vol_PDs{param_idx} = current_pd;
end

vol_frac_PDs = cell(1, length(v_means_fracs));
for param_idx = 1:length(vol_frac_PDs)
    current_mean = v_means_fracs(param_idx);
    current_sd = v_sds_fracs(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    vol_frac_PDs{param_idx} = current_pd;
end


%% Flows

% raw flows
q_pl_mean   = 0.15/1000; % Pleural fluid flow [L/h]
q_bELF_mean = (5.25 * 100) * (5.75 * 10^-5) * (1/1000) * 3600; % from Himstedt et al.
q_aELF_mean = (171 * 100) * (5.75 * 10^-5) * (1/1000) * 3600; % from Himstedt et al.

% flow fracs
q_la_mean_frac  = 0.06;
q_sp_mean_frac  = 77/5200; % spleen flow in L/h
q_gu_mean_frac  = 1100/5200; % gut flow in L/h
q_br_mean_frac  = 0.12;
q_hr_mean_frac  = 0.04;
q_ad_mean_frac  = 0.05;
q_mu_mean_frac  = 0.17;
q_bo_mean_frac  = 0.05;
q_sk_mean_frac  = 0.05;
q_oth_mean_frac = 0.04365;
q_kd_mean_frac  = 0.19;

% handle fracs and actual means separately for later processing
q_means = [q_pl_mean, q_bELF_mean, q_aELF_mean];

q_means_fracs = [q_la_mean_frac, q_sp_mean_frac, q_gu_mean_frac, ...
                q_br_mean_frac, q_hr_mean_frac, q_ad_mean_frac, ...
                q_mu_mean_frac, q_bo_mean_frac, q_sk_mean_frac, ...
                q_oth_mean_frac, q_kd_mean_frac];

q_sds = q_means .* q_cv;
q_sds_fracs = q_means_fracs .* q_cv;

% create cell arrays of truncated probability dists for each flow param
flow_PDs = cell(1, length(q_means));
for param_idx = 1:length(flow_PDs)
    current_mean = q_means(param_idx);
    current_sd = q_sds(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    flow_PDs{param_idx} = current_pd;
end

flow_frac_PDs = cell(1, length(q_means_fracs));
for param_idx = 1:length(flow_frac_PDs)
    current_mean = q_means_fracs(param_idx);
    current_sd = q_sds_fracs(param_idx);

    current_pd = makedist("Normal", "mu", current_mean, "sigma", current_sd);
    current_pd = truncate(current_pd, current_mean - 3 * current_sd, current_mean + 3 * current_sd);

    flow_frac_PDs{param_idx} = current_pd;
end


end