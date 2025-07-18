function plotTimecourses(drug, oral_dose, lung_dose, compt_to_plot, n_days, days_to_plot, oral_dose_freq, lung_dose_freq, ts, Cs_oraldose, Cs_lungdose)
% PLOTTIMECOURSES - Plots the raw concentration-time courses of each
% patient, as well as 10th, 50th, and 90th percentiles for the patient
% population for both an oral and lung dosing method.
%
% DESCRIPTION:
%   This function takes the concentration-time courses from solving the
%   model ODEs and plots them. Resulting figures show a side-by-side
%   comparison of dosing methods for the same day in the same compartment.
%   Both a plot of raw patient data as well as the population percentiles
%   are generated for each compartment, with single patient being
%   represented by a single colored line in the respective plot.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - oral_dose (int): The oral dose amt. of the drug (in mg)
% - lung_dose (int): The lung dose amt. of the drug (in mg)
% - compt_to_plot (str): The name of the compartment to be plotted
% - n_days (int): The number of days contained within the entire 
%   concentration-time course
% - days_to_plot (int): The number of days (backwards from the end of the
%   timecourse) that are to be included in the plot
% - oral_dose_freq (int): The number of times per day an oral dose is to be
%   administered
% - lung_dose_freq (int): The number of times per day an inhaled/lung dose
%   is to be administered
% - ts (double array): Contains the timepoints at which each concentration
%   is calculated and plotted
% - C_oraldose (matrix): Contains full concentration courses for each
%   compartment as its columns for an oral dose of medication
% - C_lungdose (matrix): Contains full concentration courses for each
%   compartment as its columns for a lung dose of medication
%
% OUTPUTS:
% - none
%
% PLOTS GENERATED:
% - [compt]_raw_pt_timecourses.png: A plot of all patient timecourses
% - [compt]_timecourses_percentiles.png: A plot of population-level
%   timecourses for the 10th, 50th, and 90th percentiles

% check for input error
if days_to_plot > n_days
    disp("Error: cannot plot more days than days calculated for")
    return
end

%% Raw timecourses

% set up figure as: oral conc | lung conc
fig = figure();
tlayout = tiledlayout(1, 2);
set(0, "DefaultFigureWindowStyle", "docked");

% oral dose
nexttile
plot(ts, Cs_oraldose, ...
    "LineWidth", 0.5);

% area(ts, Cs_oraldose, ...
%    "LineWidth", 0.5, ...
%    "FaceColor", "#AECCE4");

xlim([(n_days - 1 * days_to_plot) * 24, n_days * 24]);
title(append("Oral Dose (", num2str(oral_dose), " mg, ", num2str(oral_dose_freq), "x/day)"));
set(gca,'FontSize', 12);
grid on;

% lung dose
nexttile
plot(ts, Cs_lungdose, ...
    "LineWidth", 0.5);

% area(ts, Cs_lungdose, ...
%    "LineWidth", 0.5, ...
%    "FaceColor", "#E78587");

xlim([(n_days - 1 * days_to_plot) * 24, n_days * 24]);
title(append("Lung Dose (", num2str(lung_dose), " mg, ", num2str(lung_dose_freq), "x/day)"));
set(gca,'FontSize', 12);
grid on;

xlabel(tlayout, "Time Post-Initiation (h)");
ylabel(tlayout, "Concentration (\mug/mL)");

if n_days < 4
    sgtitle(append(compt_to_plot, " Concentration, Day ", num2str(n_days)), "FontSize", 20);

else
    sgtitle(append(compt_to_plot, " Concentration, Steady State"), "FontSize", 20);

end

saveas(fig, append("Outputs/", drug, "/Figures/", compt_to_plot, "_raw_pt_timecourses.png"));


%% Percentiles

% set up figure as: oral conc | lung conc
fig = figure();
tlayout = tiledlayout(1, 2);
set(0, "DefaultFigureWindowStyle", "docked");

% oral dose
nexttile
oral_p10 = prctile(Cs_oraldose, 10, 2);
oral_p50 = prctile(Cs_oraldose, 50, 2);
oral_p90 = prctile(Cs_oraldose, 90, 2);

hold on;
plot(ts, oral_p10, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "Black", ...
    "DisplayName", "10th  percentile");

plot(ts, oral_p50, ...
    "LineWidth", 2, ...
    "LineStyle", "-", ...
    "Color", "Black", ...
    "DisplayName", "50th  percentile");

plot(ts, oral_p90, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "Black", ...
    "DisplayName", "90th  percentile");

xlim([(n_days - 1 * days_to_plot) * 24, n_days * 24]);
title(append("Oral Dose (", num2str(oral_dose), " mg, ", num2str(oral_dose_freq), "x/day)"));
set(gca,'FontSize', 12);
legend show;
grid on;
hold off;

% lung dose
nexttile
lung_p10 = prctile(Cs_lungdose, 10, 2);
lung_p50 = prctile(Cs_lungdose, 50, 2);
lung_p90 = prctile(Cs_lungdose, 90, 2);

hold on;
plot(ts, lung_p10, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "Black", ...
    "DisplayName", "10th  percentile");

plot(ts, lung_p50, ...
    "LineWidth", 2, ...
    "LineStyle", "-", ...
    "Color", "Black", ...
    "DisplayName", "50th  percentile");

plot(ts, lung_p90, ...
    "LineWidth", 2, ...
    "LineStyle", "--", ...
    "Color", "Black", ...
    "DisplayName", "90th  percentile");
xlim([(n_days - 1 * days_to_plot) * 24, n_days * 24]);
title(append("Lung Dose (", num2str(lung_dose), " mg, ", num2str(lung_dose_freq), "x/day)"));
set(gca,'FontSize', 12);
legend show;
grid on;
hold off;

xlabel(tlayout, "Time Post-Initiation (h)");
ylabel(tlayout, "Concentration (\mug/mL)");

if n_days < 4
    sgtitle(append(compt_to_plot, " Concentration, Day ", num2str(n_days)), "FontSize", 20);

else
    sgtitle(append(compt_to_plot, " Concentration, Steady State"), "FontSize", 20);

end

saveas(fig, append("Outputs/", drug, "/Figures/", compt_to_plot, "_timecourses_percentiles.png"));

end