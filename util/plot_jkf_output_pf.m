
%% Joined 2IFC fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\2IFC\JKF\joined_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission'};

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    load(fullfile(data_dir, sprintf('d%s_2ifc_jkf.mat', subject)));
    
    for icnd = 1:3

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,icnd,1)), squeeze(par_mat(:,1,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,icnd,3)), squeeze(par_mat(:,1,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,icnd,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,icnd,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,icnd,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,icnd,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end





%% Separate 2IFC fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\2IFC\JKF\sep_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission'};

cond_vect = [3,6,9];

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    for icnd = 1:3
        
        load(fullfile(data_dir, sprintf('d%s_2ifc_jkf_cond%i.mat', subject, cond_vect(icnd))));

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,1)), squeeze(par_mat(:,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,3)), squeeze(par_mat(:,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end






%% Joined YN fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\yesno\JKF\joined_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission', '2F2S'};

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    load(fullfile(data_dir, sprintf('d%s_yn_jkf.mat', subject)));
    
    for icnd = 1:4

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,icnd,1)), squeeze(par_mat(:,1,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,icnd,3)), squeeze(par_mat(:,1,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,icnd,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,icnd,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,icnd,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,icnd,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end






%% Separate YN fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\yesno\JKF\sep_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission', '2F2S'};

cond_vect = [3,6,8,9];

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    for icnd = 1:4
        
        load(fullfile(data_dir, sprintf('d%s_yn_jkf_cond%i.mat', subject, cond_vect(icnd))));

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,1)), squeeze(par_mat(:,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,3)), squeeze(par_mat(:,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end






%% Joined YN-pooled fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\yesno_pooled\JKF\joined_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission'};

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    load(fullfile(data_dir, sprintf('d%s_yn_pool_jkf.mat', subject)));
    
    for icnd = 1:3

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,icnd,1)), squeeze(par_mat(:,1,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,icnd,3)), squeeze(par_mat(:,1,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,icnd,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,icnd,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,icnd,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,icnd,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end





%% Separate YN-pooled fit
subjvect = {'701', '702', '703', '704', '705', '706', '708', '709', '712', '714', ...
            '715', '716', '717', '718', '719', '720', '722', '725', '726', '727'};

data_dir = 'Z:\dfi_experiment_data\data\experiment\yesno_pooled\JKF\sep_fit';

labelvect = {'VisualOnly', 'Fusion', 'Fission'};

cond_vect = [3,6,9];

for isubj = 1:20
    
    subject = subjvect{isubj};
    
    for icnd = 1:3
        
        load(fullfile(data_dir, sprintf('d%s_yn_pool_jkf_cond%i.mat', subject, cond_vect(icnd))));

        fh = figure('color', 'w', 'position', [0 0 1450 700]);
        subplot(231)
        hist3([squeeze(par_mat(:,1)), squeeze(par_mat(:,2))])
        xlabel('Threshold'); ylabel('Slope');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Threshold and slope'); colormap(haxby)
        plotspecs
        subplot(234)
        hist3([squeeze(par_mat(:,3)), squeeze(par_mat(:,4))])
        xlabel('Guess-rate'); ylabel('Lapse-rate');
        set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        title('Guess-rate and lapse-rate'); colormap(haxby)
        plotspecs
        subplot(232)
        hist(squeeze(par_mat(:,1)),25)
        xlabel('Threshold'); ylabel('Count');
        title('Threshold marginal distribution'); colormap(haxby)
        plotspecs
        subplot(233)
        hist(squeeze(par_mat(:,2)),25)
        xlabel('Slope'); ylabel('Count');
        title('Slope marginal distribution'); colormap(haxby)
        plotspecs
        subplot(235)
        hist(squeeze(par_mat(:,3)),25)
        xlabel('Guess-rate'); ylabel('Count');
        title('Guess-rate marginal distribution'); colormap(haxby)
        plotspecs
        subplot(236)
        hist(squeeze(par_mat(:,4)),25)
        xlabel('Lapse-rate'); ylabel('Count');
        title('Lapse-rate marginal distribution'); colormap(haxby)
        plotspecs
        suptitle(sprintf('Participant %s, condition %s', subject, labelvect{icnd}));

        saveas(fh, fullfile(data_dir, sprintf('fig_params_cond%s_%s.emf',labelvect{icnd},subject)));
        close all
    
    end

end



% eof






