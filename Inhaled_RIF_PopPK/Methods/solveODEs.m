function [Cs_oral, Cs_lung] = solveODEs(drug, params, ncompts_total, n_days, oral_dose_freq, lung_dose_freq)
% SOLVEODES - Sets up and solves the corresponding system of ODEs for a 
% certain drug.
%
% DESCRIPTION:
%   This function takes as input relevant drug and model parameters and
%   calls a drug-specific system of differential equations in order to
%   solve them and return concentration vs. time arrays for both an oral
%   and lung dosing method of the drug.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - params (cell array): A list of drug-specific parameters (dosing,
%   clearance rates, etc.)
% - ncompts_total (int): The number of compartments for the drug model
% - n_days (int): The number of days for which to calculate
%   concentration-time courses
% - oral_dose_freq (int): The number of times per day an oral dose is to be
%   administered
% - lung_dose_freq (int): The number of times per day an inhaled/lung dose
%   is to be administered
% 
% OUTPUTS:
% - Cs_oral (matrix): Contains full concentration courses for each
%   compartment as its columns for an oral dose of medication
% - Cs_lung (matrix): Contains full concentration courses for each
%   compartment as its columns for a lung dose of medication

options = odeset("RelTol", 1e-6, "AbsTol", 1e-8);

%% Unpack parameters
if drug == "RIF"
    oral_dose = params{1};
    lung_dose = params{2};
    ka_oral = params{3};
    CL = params{4};
    fR = params{5};
    kr = params{6};
    pt = params{7};
    phys = params{8};
    kF = params{9};
    ka_lung = params{10};
    effRB = params{11};
    effRA = params{12};
    br_frac = params{13};
    tstep = params{14};

else
    disp("Invalid drug specified");
    return

end

%% Set up matrices to track solutions
C0_oraldose = zeros(1, ncompts_total); C0_oraldose(end) = oral_dose;
C0_lungdose = zeros(1, ncompts_total); C0_lungdose(end) = lung_dose;

timepts_oral = 0:tstep:(24 / oral_dose_freq);
timepts_lung = 0:tstep:(24 / lung_dose_freq);

Cs_oral = zeros((length(timepts_oral) - 1) * oral_dose_freq * n_days, ncompts_total);
Cs_lung = zeros((length(timepts_lung) - 1) * lung_dose_freq * n_days, ncompts_total);

%% Solve ODEs
if drug == "RIF"
    % oral eqs
    for dose_idx = 1:(n_days * oral_dose_freq)
        
        [~, Cs_oral_temp] = ode23s(@(t, C) RIF_oral_ODEs(t, C, ka_oral, kr, kF, ...
                            CL, fR, phys, pt), timepts_oral, C0_oraldose, options);
        C0_oraldose = [Cs_oral_temp(end, 1:(ncompts_total - 1))'; 
                                    Cs_oral_temp(end, ncompts_total) + oral_dose];
        Cs_oral_temp(end, :) = []; % remove initial condition

        % store result
        temp_row_start = (dose_idx - 1) * (length(timepts_oral) - 1) + 1;
        temp_row_end = dose_idx * (length(timepts_oral) - 1);
        Cs_oral(temp_row_start:temp_row_end, :) = Cs_oral_temp;
    end

    % lung eqs
    for dose_idx = 1:(n_days * lung_dose_freq)
        
        [~, Cs_lung_temp] = ode23s(@(t, C) RIF_lung_ODEs(t, C, ka_lung, kr, kF, effRB, ...
                            effRA, br_frac, CL, fR, phys, pt), timepts_lung, C0_lungdose, options);
        C0_lungdose = [Cs_lung_temp(end, 1:(ncompts_total - 1))'; 
                                    Cs_lung_temp(end, ncompts_total) + lung_dose];
        Cs_lung_temp(end, :) = []; % remove initial condition

        % store result
        temp_row_start = (dose_idx - 1) * (length(timepts_lung) - 1) + 1;
        temp_row_end = dose_idx * (length(timepts_lung) - 1);
        Cs_lung(temp_row_start:temp_row_end, :) = Cs_lung_temp;
    end

else
    disp("Invalid drug specified");
    return

end


end