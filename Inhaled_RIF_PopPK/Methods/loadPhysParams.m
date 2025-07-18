function [phys, pt] = loadPhysParams(drug, vol_PDs, vol_frac_PDs, flow_PDs, flow_frac_PDs, fixed_weight)
% LOADPHYSPARAMS: Calculates random physiological parameter values based on
% given probability distributions
%
% DESCRIPTION:
%   This function repeatedly samples from a given array of probability
%   distributions in order to get a set of random parameters that describes 
%   a certain patient. It handles the distributions for raw values and
%   distributions that are for fractions of body weight or cardiac output
%   separately, making sure to renormalize the sampled fractions so their 
%   sum is equal to 1.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - vol_PDs (cell array): Contains all probability distributions for raw
%   volume parameters not based on body weight
% - vol_frac_PDs (cell array): Contains all probability distributions for
%   volume parameters that are calculated as fractions of body weight 
% - flow_PDs (cell array): Contains all probability distributions for raw
%   flow parameters not based on cardiac output
% - flow_frac_PDs (cell array): Contains all probability distributions for
%   flow parameters that are calculated as fractions of cardiac output
% - fixed_weight (logical): Decides whether or not to sample weight as well
%   (set to 1 when using population means in optimize_regimen_[drug], 0 
%   otherwise)
%
% OUTPUTS:
% - phys (struct): Contains all relevant patient physiological volume and 
%   flow parameters
% - pt (struct): Contains all relevant partition coefficients (not randomly
%   sampled)

if ~fixed_weight
    body_weight_PD = makedist("Normal", "mu", 70, "sigma", 1.10); % SD from Dugas et al.
    body_weight_PD = truncate(body_weight_PD, 70 - 3 * 1.10, 70 + 3 * 1.10);
    BW = random(body_weight_PD);
else
    BW = 70;
end


%% Volumes

% sample all params from their probability dists
rand_params = zeros(1, length(vol_PDs));
for param_idx = 1:length(rand_params)
    rand_params(param_idx) = random(vol_PDs{param_idx});
end

rand_params_fracs = zeros(1, length(vol_frac_PDs));
for param_idx = 1:length(rand_params_fracs)
    rand_params_fracs(param_idx) = random(vol_frac_PDs{param_idx});
end

% renormalize fractions
rand_params_fracs = rand_params_fracs / sum(rand_params_fracs);
rand_params = [rand_params, rand_params_fracs .* BW];

% unpack and place into struct
phys.V.V = rand_params(1);
phys.V.A = rand_params(2);
phys.V.LN = rand_params(3);
phys.V.Lu = rand_params(4);
phys.V.Brain = rand_params(5);
phys.V.Heart = rand_params(6);
phys.V.Adipose = rand_params(7);
phys.V.Muscle = rand_params(8);
phys.V.Skin = rand_params(9);
phys.V.Kidney = rand_params(10);
phys.V.Bone = rand_params(11);
phys.V.Spleen = rand_params(12);
phys.V.Gut = rand_params(13);
phys.V.Liver = rand_params(14);
phys.V.Others = rand_params(15);

phys.V.Pl = rand_params(16);

%% Flows

QC = 5200 / 1000 * 60;

% sample all params from their probability dists
rand_params = zeros(1, length(flow_PDs));
for param_idx = 1:length(rand_params)
    rand_params(param_idx) = random(flow_PDs{param_idx});
end

rand_params_fracs = zeros(1, length(flow_frac_PDs));
for param_idx = 1:length(rand_params_fracs)
    rand_params_fracs(param_idx) = random(flow_frac_PDs{param_idx});
end

% renormalize fractions
rand_params_fracs = rand_params_fracs / sum(rand_params_fracs);
rand_params = [rand_params, rand_params_fracs .* QC];

% unpack and place intro struct
phys.Q.total = QC;

phys.Q.Pl   = rand_params(1);
phys.Q.bELF = rand_params(2);
phys.Q.aELF = rand_params(3);
phys.Q.LA   = rand_params(4);
phys.Q.Sp   = rand_params(5);
phys.Q.Gu   = rand_params(6);
phys.Q.Br   = rand_params(7);
phys.Q.Hr   = rand_params(8);
phys.Q.Ad   = rand_params(9);
phys.Q.Mu   = rand_params(10);
phys.Q.Bo   = rand_params(11);
phys.Q.Sk   = rand_params(12);
phys.Q.Oth  = rand_params(13);
phys.Q.Kd   = rand_params(14);
phys.Q.Li   = phys.Q.Sp + phys.Q.Gu + phys.Q.LA;


%% Lymph flow in L/h
lymph_total = 8 / 24; % L/day to L/h
phys.L.Lu   = lymph_total * 0.03;
phys.L.Br   = lymph_total * 0.0105;
phys.L.Hr   = lymph_total * 0.01;
phys.L.Ad   = lymph_total * 0.128;
phys.L.Mu   = lymph_total * 0.16;
phys.L.Bo   = 0;
phys.L.Sk   = lymph_total * 0.0703;
phys.L.Oth  = lymph_total * 0.0562;
phys.L.Kd   = lymph_total * 0.085;
phys.L.Li   = lymph_total * 0.33;
phys.L.Sp   = 0;
phys.L.Gu   = lymph_total * 0.12;
phys.L.LN   = lymph_total;

%% Partition coefficients
if strcmpi(drug, "RIF")
    pt.Lu             = 1.7115;
    pt.LN             = 1.2081;
    pt.Tissue.Brain   = 0.2285;
    pt.Tissue.Adipose = 0.1885;
    pt.Tissue.Heart   = 1.0158;
    pt.Tissue.Muscle  = 0.6949;
    pt.Tissue.Skin    = 0.6265;
    pt.Tissue.Others  = 1.047;
    pt.Tissue.Bone    = 0.3157;
    pt.Tissue.Spleen  = 1.3950;
    pt.Tissue.Kidney  = 2.1725;
    pt.Tissue.Gut     = 1.0781;
    pt.Tissue.Liver   = 1.9646;

else
    error('Partition coefficients for %s not defined.', drug);

end


end