
% I cannot compute one tailed BFs in matlab, so simply get the one sided
% t-values and save them for computation of BFs in R, where necessary:

condvect = {'v2', 'fus', 'fis'};
diff_score = nan(3,13);
for icond = 1:length(condvect)
    load(fullfile('H:\dfi_experiment_figures\Paper_figures\iAF\phase\ynt\dp_diff_phase_opp',  ...
        sprintf('figure_data_PO4_O2_PO8_%s.mat', condvect{icond})));
    diff_score(icond,:) = fslide_GA(1) - fslide_GA(2);
    
end

T = array2table(diff_score','VariableNames',condvect);

writetable(T,fullfile('H:\dfi_experiment_figures\Paper_figures\iAF\phase\ynt\dp_diff_phase_opp\tvals_for_R'))



