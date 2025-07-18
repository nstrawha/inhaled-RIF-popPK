function [AUCs_oral_store, AUCs_lung_store, Cmaxs_oral_store, Cmaxs_lung_store] = trackPKMetrics(drug, compt_list, relevant_compts, toxic_compts, odose, odose_freq, ldose, ldose_freq, n_days, tstep, ts, n_pts, Cs_oral, Cs_lung)
% TRACKPKMETRICS - Records, analyzes, and stores AUC and Cmax for a
% comparison between oral and lung dosing methods.
%
% DESCRIPTION:
%   This function takes the concentration-time courses for each patient and
%   computes important PK metrics (AUC and Cmax) from them. It then
%   performs a series of pairwise 2-sample t-tests between oral and lung
%   dosing methods for each metric in order to determine if a statistically
%   significant difference exists between the two, recording the resulting
%   p-value and effect size.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - compt_list (str array): A list of all compartments for the model
% - relevant_compts (str array): A list of all compartments for which
%   metrics are to be calculated and analyzed
% - odose (int): The oral dose amt. of the drug (in mg)
% - odose_freq (int): The number of times per day an oral dose is to be
%   administered
% - ldose (int): The lung dose amt. of the drug (in mg)
% - ldose_freq (int): The number of times per day an inhaled/lung dose
%   is to be administered
% - n_days (int): The number of days for which concentration-time courses
%   have been calculated
% - tstep (double): The time interval at which each concentration-time
%   course is calculated for
% - ts (double array): Contains the timepoints at which each concentration
%   is calculated and plotted
% - n_pts (int): The number of patients for which concentration-time
%   courses have been calculated
% - Cs_oral (cell array): Contains full concentration courses for each
%   patient as matrices whose columns represent different compartments for 
%   an oral dose of medication
% - Cs_lung (cell array): Contains full concentration courses for each
%   patient as matrices whose columns represent different compartments for 
%   a lung dose of medication
%
% OUTPUTS:
% - AUCs_oral_store (cell array): Contains AUCs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - AUCs_lung_store (cell array): Contains AUCs for all patients given a 
%   lung dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_oral_store (cell array): Contains Cmaxs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_lung_store (cell array): Contains Cmaxs for all patients given a
%   lung dose in each cell, with each cell representing a different 
%   compartment
%
% FILES GENERATED:
% - popPK_metric_analysis.xlsx: Contains the means and standard deviations
%   for each PK metric in each compartment Also contains the p-value and
%   effect size from each t-test between oral and lung dosing methods.


%% Set up cell array formatting
cells_store = cell(1, length(relevant_compts));
tbl_format = cell(8, 4);

tbl_format{1, 1} = append("Day ", num2str(n_days), ", n = ", num2str(n_pts));

tbl_format{2, 1} = "Oral Dose";
tbl_format{3, 1} = append(num2str(odose), " mg, ", num2str(odose_freq), " x/day");
tbl_format{4, 1} = "Lung Dose";
tbl_format{5, 1} = append(num2str(ldose), " mg, ", num2str(ldose_freq), " x/day");
tbl_format{7, 1} = "Better Dose";

tbl_format{1, 2} = "Metric Type";
tbl_format{2, 2} = "AUC_24";
tbl_format{3, 2} = "C_max (ug/mL)";
tbl_format{4, 2} = "AUC_24";
tbl_format{5, 2} = "C_max (ug/mL)";
tbl_format{7, 2} = "AUC_24";
tbl_format{8, 2} = "C_max (ug/mL)";

tbl_format{1, 3} = "Mean";
tbl_format{1, 4} = "SD";
tbl_format{6, 3} = "p-value";
tbl_format{6, 4} = "Effect Size";

last_day_start = length(ts) - 24 / tstep + 1;

% set up arrays to store metrics for each compartment to return
AUCs_oral_store = cell(1, length(relevant_compts));
AUCs_lung_store = cell(1, length(relevant_compts));
Cmaxs_oral_store = cell(1, length(relevant_compts));
Cmaxs_lung_store = cell(1, length(relevant_compts));


%% Iterate through compartments
for compt_idx = 1:length(relevant_compts)

    current_compt = relevant_compts{compt_idx};
    current_cell = tbl_format;
    [idx_to_calc, ~] = find(string(compt_list) == current_compt);

    % pull compt specific timecourses
    current_Cs_oral = cell2mat(cellfun(@(x) x(:, idx_to_calc), Cs_oral, "UniformOutput", false));
    current_Cs_lung = cell2mat(cellfun(@(x) x(:, idx_to_calc), Cs_lung, "UniformOutput", false));

    % calculate metrics
    AUCs_oral = trapz(ts(last_day_start:end), current_Cs_oral(last_day_start:end, :));
    AUCs_lung = trapz(ts(last_day_start:end), current_Cs_lung(last_day_start:end, :));

    Cmaxs_oral = max(current_Cs_oral(last_day_start:end, :));
    Cmaxs_lung = max(current_Cs_lung(last_day_start:end, :));

    % calculate metric info for comparison
    AUCs_oral_mean  = mean(AUCs_oral);
    AUCs_oral_sd    = std(AUCs_oral);
    AUCs_lung_mean  = mean(AUCs_lung);
    AUCs_lung_sd    = std(AUCs_lung);

    Cmaxs_oral_mean  = mean(Cmaxs_oral);
    Cmaxs_oral_sd    = std(Cmaxs_oral);
    Cmaxs_lung_mean  = mean(Cmaxs_lung);
    Cmaxs_lung_sd    = std(Cmaxs_lung);

    % calculate better method for each metric
    if ~ismember(current_compt, toxic_compts)

        if AUCs_oral_mean >= AUCs_lung_mean
            better_AUC = "Oral";
        else
            better_AUC = "Lung";
        end

        if Cmaxs_oral_mean >= Cmaxs_lung_mean
            better_Cmax = "Oral";
        else
            better_Cmax = "Lung";
        end

    else
        if AUCs_oral_mean > AUCs_lung_mean
            better_AUC = "Lung";
        else
            better_AUC = "Oral";
        end

        if Cmaxs_oral_mean > Cmaxs_lung_mean
            better_Cmax = "Lung";
        else
            better_Cmax = "Oral";
        end
    end

    % calculate effect sizes 
    AUCs_effsize  = meanEffectSize(AUCs_oral, AUCs_lung);
    Cmaxs_effsize = meanEffectSize(Cmaxs_oral, Cmaxs_lung);

    % perform hypothesis tests
    [~, AUC_p]  = ttest(AUCs_oral, AUCs_lung);
    [~, Cmax_p] = ttest(Cmaxs_oral, Cmaxs_lung);

    % record results
    % means
    current_cell{2, 3} = round(AUCs_oral_mean, 2);
    current_cell{3, 3} = round(Cmaxs_oral_mean, 2);
    current_cell{4, 3} = round(AUCs_lung_mean, 2);
    current_cell{5, 3} = round(Cmaxs_lung_mean, 2);

    % SDs
    current_cell{2, 4} = round(AUCs_oral_sd, 2);
    current_cell{3, 4} = round(Cmaxs_oral_sd, 2);
    current_cell{4, 4} = round(AUCs_lung_sd, 2);
    current_cell{5, 4} = round(Cmaxs_lung_sd, 2);

    % tests
    current_cell{7, 3}  = append(better_AUC, ", p = ", num2str(AUC_p));
    current_cell{8, 3} = append(better_Cmax, ", p = ", num2str(Cmax_p, 3));

    current_cell{7, 4}  = abs(round(AUCs_effsize.Effect, 2));
    current_cell{8, 4} = abs(round(Cmaxs_effsize.Effect, 2));

    cells_store{compt_idx} = current_cell;

    % store metrics
    AUCs_oral_store{compt_idx} = AUCs_oral;
    AUCs_lung_store{compt_idx} = AUCs_lung;
    Cmaxs_oral_store{compt_idx} = Cmaxs_oral;
    Cmaxs_lung_store{compt_idx} = Cmaxs_lung;

    % plotting
    fig = figure();
    tiledlayout(1, 2);

    % AUC
    nexttile
    hold on;
    current_bar = bar(1:2, [AUCs_oral_mean, AUCs_lung_mean], 0.5);
    current_bar.FaceColor = "flat";
    current_bar.CData(1, :) = [0 0 1];
    current_bar.CData(2, :) = [1 0 0];
        
    bar_heights = [AUCs_oral_mean, AUCs_lung_mean];
    label_text = append("Effect size: ", num2str(-round(AUCs_effsize.Effect, 2)));
    text(2, bar_heights(2) + 1.25, label_text, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    errorbar(1:2, [AUCs_oral_mean, AUCs_lung_mean], ...
                                [AUCs_oral_sd ./ sqrt(n_pts), AUCs_lung_sd ./ sqrt(n_pts)], ...
                                "LineStyle", "none", ...
                                "LineWidth", 0.75, ...
                                "Color", "Black", ...
                                "CapSize", 25);

    lower_AUC = min([AUCs_oral_mean, AUCs_lung_mean]);
    higher_AUC = max([AUCs_oral_mean, AUCs_lung_mean]);
    ylims = [floor(lower_AUC / 5) * 5 - 1, ceil(higher_AUC / 5) * 5 + 1];
    ylim(ylims);
    xticks(1:2)
    xticklabels(["Oral", "Lung"]);
    ylabel("AUC (\mug*h/L)")
    title("AUC", "FontSize", 15);
    hold off;

    % Cmax
    nexttile
    hold on;
    current_bar = bar(1:2, [Cmaxs_oral_mean, Cmaxs_lung_mean], 0.5);
    current_bar.FaceColor = "flat";
    current_bar.CData(1, :) = [0 0 1];
    current_bar.CData(2, :) = [1 0 0];

    bar_heights = [Cmaxs_oral_mean, Cmaxs_lung_mean];
    label_text = append("Effect size: ", num2str(-round(Cmaxs_effsize.Effect, 2)));
    text(2, bar_heights(2) + 1.25, label_text, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    errorbar(1:2, [Cmaxs_oral_mean, Cmaxs_lung_mean], ...
                                [Cmaxs_oral_sd ./ sqrt(n_pts), Cmaxs_lung_sd ./ sqrt(n_pts)], ...
                                "LineStyle", "none", ...
                                "LineWidth", 0.75, ...
                                "Color", "Black", ...
                                "CapSize", 25);

    lower_Cmax = min([Cmaxs_oral_mean, Cmaxs_lung_mean]);
    higher_Cmax = max([Cmaxs_oral_mean, Cmaxs_lung_mean]);
    ylims = [floor(lower_Cmax) - 5, ceil(higher_Cmax) + 5];
    ylim(ylims);
    xticks(1:2)
    xticklabels(["Oral", "Lung"]);
    ylabel("C_{max} (\mug/L)")
    title("C_{max}", "FontSize", 15);
    hold off;

    sgtitle(append("PK Metric Comparison for ", current_compt, " Compartment"), "FontSize", 20);

    saveas(fig, append("Outputs/", drug, "/Figures/", current_compt, "_metric_comparisons.png"))

end


%% Write output file
for compt_idx = 1:length(relevant_compts)
    writecell(cells_store{compt_idx}, ...
        append("Outputs/", drug, "/Tables/popPK_metric_analysis.xlsx"), ...
        "Sheet", relevant_compts{compt_idx});
end
     writecell(cells_store{compt_idx}, append("Outputs/", drug, "/Tables/popPK_analysis_day", num2str(n_days), ".xlsx"), "Sheet", relevant_compts{compt_idx});
end