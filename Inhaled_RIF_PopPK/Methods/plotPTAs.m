function plotPTAs(drug, nontoxic_compts, AUCs_oral, AUCs_lung, Cmaxs_oral, Cmaxs_lung, MIC_dist, AUC_target, Cmax_target, oral_dose, oral_dose_freq, lung_dose, lung_dose_freq, n_pts, n_days)
% PLOTPTAS - Calculates and plots the probability of target attainment (for
% AUC/MIC and Cmax/MIC targets) for an oral vs. lung dose of medication.
%
% DESCRIPTION:
%   This function uses the PK metrics calculated for all patients and a
%   drug-TB MIC distribution in order to calculate and plot the probability 
%   of target attainment (PTA) for each metric in each compartment.
%
% INPUTS:
% - drug (str): An all-caps three-letter identifier of the relevant drug
% - nontoxic compts (str array): Contains the names of all nontoxic
%   compartments for which PTA is to be calculated and plotted
% - AUCs_oral (cell array): Contains AUCs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - AUCs_lung (cell array): Contains AUCs for all patients given a 
%   lung dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_oral (cell array): Contains Cmaxs for all patients given an 
%   oral dose in each cell, with each cell representing a different 
%   compartment
% - Cmaxs_lung (cell array): Contains Cmaxs for all patients given a
%   lung dose in each cell, with each cell representing a different 
%   compartment
% - MIC_dist (cell array): Contains arrays formatted like [MIC, # isolates]
%   describing the drug-TB MIC distribution
% - AUC_target (double): The AUC/MIC target for the drug
% - Cmax_target (double): The Cmax/MIC target for the drug
% - oral_dose (int): The oral dose amt. of the drug (in mg)
% - oral_dose_freq (int): The number of times per day an oral dose is to be
%   administered
% - lung_dose (int): The lung dose amt. of the drug (in mg)
% - lung_dose_freq (int): The number of times per day an inhaled/lung dose
%   is to be administered
% - n_pts (int): The number of patients for which concentration-time
%   courses have been calculated
% - n_days (int): The number of days for which concentration-time courses
%   have been calculated

% unpack MIC dist info
concs = cell2mat(cellfun(@(x) x(:, 1), MIC_dist, "UniformOutput", false));
isolates = cell2mat(cellfun(@(x) x(:, 2), MIC_dist, "UniformOutput", false));

MICs_to_test = 0.001:0.001:25; % start at 0.001 ug/mL

% set up storage
PTA_AUC_storage = cell(1, length(nontoxic_compts));
PTA_Cmax_storage = cell(1, length(nontoxic_compts));
CFR_AUC_storage = cell(1, length(nontoxic_compts));
CFR_Cmax_storage = cell(1, length(nontoxic_compts));

% find the continuous lognormal probabiliity dist for MICs
dummy_data = repelem(concs, isolates);
MIC_PD = fitdist(dummy_data, "Lognormal");

%% Plot MICs and the estimated full dist
log_MICs = log(dummy_data);
MIC_PD_mu = mean(dummy_data);
MIC_PD_sigma = std(log_MICs);

MIC_dist_to_plot = makedist("Lognormal", "mu", MIC_PD_mu, "sigma", MIC_PD_sigma);

fig = figure();
hold on;

scatter(concs, isolates / max(isolates), ...
    "DisplayName", "MIC Distribution (Normalized)", ...
    "MarkerFaceColor", "Black", ...
    "MarkerEdgeColor", "Black");

plot(MICs_to_test, pdf(MIC_dist_to_plot, MICs_to_test) / max(pdf(MIC_dist_to_plot, MICs_to_test)), ...
    "DisplayName", "Estimated Population MIC Dist. (Normalized)", ...
    "LineWidth", 2, ...
    "Color", "Black");

yline(0, "LineWidth", 2, "HandleVisibility", "off");
xlim([0.008, max(MICs_to_test)]);
ylim([-0.05, 1.05]);

set(gca, "XScale", "log");
xlabel("Minimum Inhibitory Concentration (\mug/mL)");
ylabel("Probability of Target Attainment (PTA)");
leg = legend("show");
leg.FontSize = 8;

title(append(drug, "/TB MIC Distribution"), "FontSize", 20);

saveas(fig, append("Outputs/", drug, "/Figures/MIC_dist.png"));

hold off;


%% Calculate PTA dists

parfor compt_idx = 1:length(nontoxic_compts)

    current_AUCs_oral = AUCs_oral{compt_idx};
    current_AUCs_lung = AUCs_lung{compt_idx};
    current_Cmaxs_oral = Cmaxs_oral{compt_idx};
    current_Cmaxs_lung = Cmaxs_lung{compt_idx};

    AUC_PTA_array = zeros(2, length(MICs_to_test));
    Cmax_PTA_array = zeros(2, length(MICs_to_test));

    % iterate through MICs
    for mic_idx = 1:length(MICs_to_test)

        current_mic = MICs_to_test(mic_idx);

        success_pts_AUC_oral = 0;
        success_pts_AUC_lung = 0;
        success_pts_Cmax_oral = 0;
        success_pts_Cmax_lung = 0;

        % iterate through patients
        for pt_idx = 1:n_pts
            pt_AUC_oral = current_AUCs_oral(pt_idx);
            pt_AUC_lung = current_AUCs_lung(pt_idx);
            pt_Cmax_oral = current_Cmaxs_oral(pt_idx);
            pt_Cmax_lung = current_Cmaxs_lung(pt_idx);

            % track AUC/MIC successes
            if pt_AUC_oral/current_mic >= AUC_target
                success_pts_AUC_oral = success_pts_AUC_oral + 1;
            end

            if pt_AUC_lung/current_mic >= AUC_target
                success_pts_AUC_lung = success_pts_AUC_lung + 1;
            end

            % track Cmax/MIC successes
            if pt_Cmax_oral/current_mic >= Cmax_target
                success_pts_Cmax_oral = success_pts_Cmax_oral + 1;
            end

            if pt_Cmax_lung/current_mic >= Cmax_target
                success_pts_Cmax_lung = success_pts_Cmax_lung + 1;
            end

        end

        AUC_PTA_array(1:2, mic_idx) = [success_pts_AUC_oral / n_pts; success_pts_AUC_lung / n_pts];
        Cmax_PTA_array(1:2, mic_idx) = [success_pts_Cmax_oral / n_pts; success_pts_Cmax_lung / n_pts];

    end

    PTA_AUC_storage{compt_idx} = AUC_PTA_array;
    PTA_Cmax_storage{compt_idx} = Cmax_PTA_array;

end


%% Plot AUC/MIC results

for compt_idx = 1:length(PTA_AUC_storage)

    current_PTA_array = PTA_AUC_storage{compt_idx};
    oral_PTAs = current_PTA_array(1, :);
    lung_PTAs = current_PTA_array(2, :);

    fig = figure();
    hold on;

    % MIC distribution
    plot(MICs_to_test, pdf(MIC_PD, MICs_to_test) / max(pdf(MIC_PD, MICs_to_test)), ...
        "DisplayName", "Continuous MIC Dist. Estimate (Normalized)", ...
        "LineWidth", 2, ...
        "Color", "Black");

    % PTAs
    plot(MICs_to_test, oral_PTAs, ...
        "DisplayName", append("Oral PTA; ", num2str(oral_dose), " mg, ", num2str(oral_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Blue");
    plot(MICs_to_test, lung_PTAs, ...
        "DisplayName", append("Lung PTA; ", num2str(lung_dose), " mg, ", num2str(lung_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Red");
    
    yline(0, "LineWidth", 2, "HandleVisibility", "off");
    xlim([min(MICs_to_test), max(MICs_to_test)]);
    ylim([-0.05, 1.05]);
    xlabel("Minimum Inhibitory Concentration (\mug/mL)");
    ylabel("Probability of Target Attainment (PTA)");
    title(append(nontoxic_compts{compt_idx}, " PTA (AUC/MIC)"), "FontSize", 20);

    set(gca, "XScale", "log");
    leg = legend("show");
    leg.FontSize = 8;
    hold off;
    
    % predict prop of isolates inhibited
    oral_inhibition_yvals = oral_PTAs .* pdf(MIC_PD, MICs_to_test);
    lung_inhibition_yvals = lung_PTAs .* pdf(MIC_PD, MICs_to_test);

    CFR_AUC_oral = trapz(MICs_to_test, oral_inhibition_yvals);
    CFR_AUC_lung = trapz(MICs_to_test, lung_inhibition_yvals);

    CFR_AUC_storage{compt_idx} = [CFR_AUC_oral, CFR_AUC_lung];

    % save figure
    saveas(fig, append("Outputs/", drug, "/Figures/", nontoxic_compts{compt_idx}, "_AUCMIC_PTAs.png"));

end


%% Plot Cmax/MIX results

for compt_idx = 1:length(PTA_Cmax_storage)

    current_PTA_array = PTA_Cmax_storage{compt_idx};
    oral_PTAs = current_PTA_array(1, :);
    lung_PTAs = current_PTA_array(2, :);

    fig = figure();
    hold on;

    % MIC distribution
    plot(MICs_to_test, pdf(MIC_PD, MICs_to_test) / max(pdf(MIC_PD, MICs_to_test)), ...
        "DisplayName", "Continuous MIC Dist. Estimate (Normalized)", ...
        "LineWidth", 2, ...
        "Color", "Black");

    % PTAs
    plot(MICs_to_test, oral_PTAs, ...
        "DisplayName", append("Oral PTA; ", num2str(oral_dose), " mg, ", num2str(oral_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Blue");
    plot(MICs_to_test, lung_PTAs, ...
        "DisplayName", append("Lung PTA; ", num2str(lung_dose), " mg, ", num2str(lung_dose_freq), " x/day"), ...
        "LineWidth", 1.5, "Color", "Red");
    
    yline(0, "LineWidth", 2, "HandleVisibility", "off");
    xlim([min(MICs_to_test), max(MICs_to_test)]);
    ylim([-0.05, 1.05]);
    xlabel("Minimum Inhibitory Concentration (\mug/mL)");
    ylabel("Probability of Target Attainment (PTA)");
    title(append(nontoxic_compts{compt_idx}, " PTA (C_{max}/MIC)"), "FontSize", 20);

    set(gca, "XScale", "log");
    leg = legend("show");
    leg.FontSize = 8;
    hold off;
    
    % predict prop of isolates inhibited
    oral_inhibition_yvals = oral_PTAs .* pdf(MIC_PD, MICs_to_test);
    lung_inhibition_yvals = lung_PTAs .* pdf(MIC_PD, MICs_to_test);

    CFR_Cmax_oral = trapz(MICs_to_test, oral_inhibition_yvals);
    CFR_Cmax_lung = trapz(MICs_to_test, lung_inhibition_yvals);

    CFR_Cmax_storage{compt_idx} = [CFR_Cmax_oral, CFR_Cmax_lung];

    % save figure
    saveas(fig, append("Outputs/", drug, "/Figures/", nontoxic_compts{compt_idx}, "_CmaxMIC_PTAs.png"));

end


%% CFR barplots
for compt_idx = 1:length(PTA_AUC_storage)
    current_compt = nontoxic_compts{compt_idx};
    
    fig = figure();
    tlayout = tiledlayout(1, 2);

    % AUC/MIC
    nexttile
    current_bar = bar(1:2, [CFR_AUC_oral * 100, CFR_AUC_lung * 100], 0.5);
    current_bar.FaceColor = "flat";
    current_bar.CData(1, :) = [0 0 1];
    current_bar.CData(2, :) = [1 0 0];

    bar_heights = [CFR_AUC_oral * 100, CFR_AUC_lung * 100];
    label_text1 = append(num2str(round(CFR_AUC_oral * 100, 2)), "%");
    label_text2 = append(num2str(round(CFR_AUC_lung * 100, 2)), "%");
    text(1, bar_heights(1) + 2, label_text1, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
    text(2, bar_heights(2) + 2, label_text2, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    ylim([0, 100]);
    xticks(1:2)
    xticklabels(["Oral", "Lung"]);
    title("CFR for AUC/MIC Target", "FontSize", 15);

    % Cmax/MIC
    nexttile
    current_bar = bar(1:2, [CFR_Cmax_oral * 100, CFR_Cmax_lung * 100], 0.5);
    current_bar.FaceColor = "flat";
    current_bar.CData(1, :) = [0 0 1];
    current_bar.CData(2, :) = [1 0 0];

    bar_heights = [CFR_Cmax_oral * 100, CFR_Cmax_lung * 100];
    label_text1 = append(num2str(round(CFR_Cmax_oral * 100, 2)), "%");
    label_text2 = append(num2str(round(CFR_Cmax_lung * 100, 2)), "%");
    text(1, bar_heights(1) + 2, label_text1, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');
    text(2, bar_heights(2) + 2, label_text2, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 10, ...
        'FontWeight', 'bold');

    ylim([0, 100]);
    xticks(1:2)
    xticklabels(["Oral", "Lung"]);
    title("CFR for C_{max}/MIC Target", "FontSize", 15);

    ylabel(tlayout, "% CFR");
    sgtitle(append(current_compt, " CFRs"), "FontSize", 20);

    saveas(fig, append("Outputs/", drug, "/Figures/", current_compt, "_CFRs_comparison.png"));

end

%% Write output files

tbl_format = cell(3, 3);

% set labels
tbl_format{1, 1} = append("Day ", num2str(n_days), ", n = ", num2str(n_pts));
tbl_format{2, 1} = "Oral Dose";
tbl_format{3, 1} = "Lung Dose";

tbl_format{1, 2} = "% CFR (AUC/MIC)";
tbl_format{1, 3} = "% CFR (Cmax/MIC)";

% iterate through compartments
for compt_idx = 1:length(nontoxic_compts)

    current_cell = tbl_format;
    current_AUC_CFRs = CFR_AUC_storage{compt_idx};
    current_Cmax_CFRs = CFR_Cmax_storage{compt_idx};

    current_cell{2, 2} = round(current_AUC_CFRs(1) * 100, 2);
    current_cell{3, 2} = round(current_AUC_CFRs(2) * 100, 2);
    current_cell{2, 3} = round(current_Cmax_CFRs(1) * 100, 2);
    current_cell{3, 3} = round(current_Cmax_CFRs(2) * 100, 2);

    % write output file
    writecell(current_cell, append("Outputs/", drug, "/Tables/popPK_CFRs.xlsx"), "Sheet", nontoxic_compts{compt_idx});

end


end