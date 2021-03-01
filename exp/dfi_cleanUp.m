function dfi_cleanUp(s, data, full, debug, dpractice)
% dfi_cleanUp(s, full) cleans up after a PTB sessions, i.e. closes the screen,
% clears the keyboard queue, returns priviledges to os etc.
% s is the setup structure used throughout the VE experiment and full
% specifies whether or not to clean up the workspace. 
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% created, April 2015
% last update, Februrary, 2016
% 
% Steffen Buergers, sbuergers@gmail.com
% 

% Eyelink edf2asc flags:
% -l or -nr   outputs left-eye data only if binocular data file 
% -r or -nl   outputs right-eye data only if binocular data file
% -sp         outputs sample raw pupil position if present 
% -sh         outputs sample HREF angle data if present 
% -sg         outputs sample GAZE data if present (default) 
% -res        outputs resolution data if present
% -vel        outputs velocity data in samples if possible 
% -s or -ne   outputs sample data only
% -e or -ns   outputs event data only 
% -nse        blocks output of start events 
% -nmsg       blocks message event output
% -neye       outputs only non-eye events (for sample-only files)
% -miss <string> replaces missing data in ASC file with <string>
% -setres <xr> <yr> uses a fixed <xr>,<yr> resolution always 
% -defres <xr> <yr> uses a default <xr>,<yr> resolution if none in file
    

    
    disp('Cleaning up...') 
    
    time = datestr(now, 'HH:MM');

    
    %% *** Clean up EyeLink ***
    
    if s.el.online
        
        % Stop writing to edf
        disp('Stop Eyelink recording...')
        Eyelink('Command', 'set_idle_mode'); 
        WaitSecs(0.5);
        Eyelink('CloseFile'); 
        
        % Copy over edf
        try
            fprintf('Receiving data file ''%s''\n', s.edfFile );
            edf_new_filename = sprintf('subj%ssess%drun%d_%s.edf',s.subj.id, ...
                                       s.subj.session, s.runid, [time(1:2),'h',time(4:5),'m']);
            status=Eyelink('ReceiveFile', s.edfFile, fullfile(s.mypath.data, edf_new_filename));
            s.edfFile = edf_new_filename;
            if status > 0
                fprintf('ReceiveFile status %d\n', status);
            end
            if 2==exist(fullfile(s.mypath.data, edf_new_filename), 'file')
                fprintf('Data file can be found in ''%s''\n', fullfile(s.mypath.data, edf_new_filename));  
            end
        catch %#ok<*CTCH>
            fprintf('Problem receiving data file ''%s''\n', s.edfFile );
        end
        fprintf('\n\nTrying to convert to ASCII...\n')
        shellfb = system(sprintf('%sedf2asc %s/%s %s', s.mypath.root, s.mypath.data, s.edfFile));
        
        % Shut down connection
        Eyelink('Shutdown'); 
    end
    
    
    
    %% *** Clean up screen ***
    
    sca
    if ~debug
        ShowCursor;
        Screen('Preference', 'VisualDebugLevel', s.disp.old.visdeb);
    end
    Screen('Preference', 'Verbosity', s.disp.old.verb);
    
    % return CPU priviledges
    Priority(s.disp.oldPriority);
    
    % Allow keyboard input to matlab
    ListenChar(1); 
    KbQueueRelease(s.kb.subj); % destroy Kb queue
    
    % Clean up sound driver
    PsychPortAudio('Close', s.haudio);
    
    
    
    %% *** Save data ***
    
    try
        disp('Saving data...')
        
        % save
        fname = fullfile(s.mypath.data, sprintf('subj%s_sess%g_run%g_%s.mat', ...
                         s.subj.id, s.subj.session, s.runid, [time(1:2),'h',time(4:5),'m']));
        save(fname, 's', 'data')
        % practice data
        if exist('dpractice', 'var')
            save(fname, 'dpractice', '-append');
        end
        
        % All of these accuracies pertain to key configuration 'AB'
        if strcmp(s.paradigm, 'yesno') 
            % add accuracy column to data
            acc = zeros(size(data(:,1)));
            acc(data.resp == 1 & ismember(data.trlid, [2,  5,  8,    4])) = 1;
            acc(data.resp == 2 & ismember(data.trlid, [3,  6,  9,    7])) = 1;
            data.acc = acc;
            if strcmp(s.cond, 'multiAttAud')
                acc = zeros(size(data(:,1)));
                acc(data.resp == 1 & ismember(data.trlid, [4 5 6])) = 1;
                acc(data.resp == 2 & ismember(data.trlid, [7 8 9])) = 1;
                data.acc = acc;
            end
        elseif strcmp(s.paradigm, '2AFC')
            acc = zeros(size(data(:,1)));
            acc(data.resp == 2 & ismember(data.trlid, [2,  5,  8,    4])) = 1;
            acc(data.resp == 1 & ismember(data.trlid, [3,  6,  9,    7])) = 1;
            data.acc = acc;
        elseif strcmp(s.paradigm, 'YN_threshold')
            acc = zeros(size(data(:,1)));
            acc(data.resp < 5 & data.resp > 0 & ismember(data.trlid, [2,  5,  8,    4])) = 1;
            acc(data.resp > 4 &                 ismember(data.trlid, [3,  6,  9,    7])) = 1;
            data.acc = acc;
            data.conf = data.resp;
            data.conf(data.resp > 4) = abs(5 - (data.conf(data.resp > 4) - 4));
        end
        
        % When key configuration is 'BA' accuracy is simply flipped
        data.acc(strcmp(data.rkeycfg, 'BA')) = ~data.acc(strcmp(data.rkeycfg, 'BA'));

        % add block number and trlid per block
        blid       = repmat((1:s.blocks.N)', 1, s.ntrls)';
        data.block = blid(:);
        data.trl   = repmat(1:s.ntrls,1, s.blocks.N)';
        data.partid= repmat(str2double(s.subj.id), [size(data,1),1]);

        % add session number and run number
        data.sess = repmat(s.subj.session, size(data, 1), 1);
        data.run  = repmat(s.runid, size(data, 1), 1);

        % add post-error, post-correct columns
        posterr = [0; ~data.acc(1:end-1)]; % let all first trials per block
        postcor = [0;  data.acc(1:end-1)]; % be post-correct
        posterr(1:s.ntrls+1:end) = 0;
        postcor(1:s.ntrls+1:end) = 1;
        data.perr = posterr;
        data.pcor = postcor;

        % save
        fname = fullfile(s.mypath.data, sprintf('subj%s_sess%g_run%g_%s.mat', ...
                         s.subj.id, s.subj.session, s.runid, [time(1:2),'h',time(4:5),'m']));
        save(fname, 's', 'data')
        % practice data
        if exist('dpractice', 'var')
            save(fname, 'dpractice', '-append');
        end
    catch
        disp('Unable to save data...')
    end
    
    
    
    %% *** Upload data ***
    
    try
        % turn wifi on
        tic
        system('networksetup -setairportpower en1 on');
        
        % Upload data
        fprintf('\n\nUploading data to google drive .....\n');
        try
            upload_path = fullfile('/Users/sxb1173/Google Drive/SIFI_interim_data', s.subj.id, s.paradigm);
            if ~exist(upload_path, 'dir'), mkdir(upload_path); end;
            if exist('dpractice', 'var')
                save(fullfile(upload_path, sprintf('data_practice_final_subj%s_%s_sess%d_run%d%',s.subj.id, s.paradigm, s.subj.session, s.runid)), 'dpractice'); 
            end
            if exist('data', 'var')
                save(fullfile(upload_path, sprintf('data_final_subj%s_%s_sess%d_run%d%',s.subj.id, s.paradigm, s.subj.session, s.runid)), 'data'); 
            end
            fprintf('.....successful!\n');
        catch
            fprintf('.....unsuccessful!\n');
        end
        toc
        
        % turn wifi off again
        system('networksetup -setairportpower en1 off')
    catch
        fprintf('\n\nUnable to upload data to Google Drive!\n\n');
    end
    
    
    %% *** Clear workspace ***
    
    if full 
        clear functions
        clearvars
        disp('Clearing functions and variables...')
    else
        disp('Leaving workspace variables untouched...')
    end
    disp('Finished.')
    diary off
    
    
end % end dfi_cleanUp

% eof








