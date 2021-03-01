function data = dfi_generate_data(s, ispractice)
% get empty dataset of SIFI paradigm for one participant
%
% USAGE:
%    data = dfi_generate_data(s)
% 
% DETAILS:
%    takes specifications from dfi_run_experiment and creates a stimulus-response
%    dataset object with randomized IV and empty DV vectors of the form:
%    resp    RT    SOA     trialtype
%
%    The optional argument 'ispractice' is a boolean that if provided will
%    use the practice settings specified in 's'. The default is not to do
%    this.
%
% INPUT:
%    s, the structure containing all specifications in this experiment
%
% OUTPUT:
%    data: dataset specified in details.
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, Oct 2015
%

% set defaults
if ~exist('ispractice', 'var')
    ispractice = 0;
end


ntrls = s.ntrls;


%% data is generated differently for YN_threshold (different SOAs for different conditions)
if strcmp(s.paradigm, 'YN_threshold')
    
%                     ntrls   = s.ntrls + 2; % for 2 nil trials per block
    data = zeros(ntrls*s.blocks.N, 6); 
    
    % for now I'm assuming I mix trials within the same block
    if ~s.randwblocks
        error('Having separate blocks for different stimulus conditions is not implemented for YN_threshold!')
    else
        nblocks = s.blocks.N;
        SOA = [s.stim.SOA_FF, s.stim.SOA_Fus, s.stim.SOA_Fis, s.stim.SOA_2av];
                
        % randomize trl condition vector
        trlid   = repmat(s.stim.trlid, [1 (s.ntrls-1)/numel(s.stim.trlid)])';
        soaid   = repmat(SOA,   [1 (s.ntrls-1)/numel(SOA)])';
        
        % add average SOA over conditions for 1 Nil trials per block:
        avg_SOA = round(mean(soaid/(1/s.disp.refr))) * (1/s.disp.refr);
        trlid = [trlid; 1];
        soaid = [sort(soaid); avg_SOA];
        soaid(trlid == 3) = s.stim.SOA_FF;
        soaid(trlid == 6) = s.stim.SOA_Fus;
        soaid(trlid == 8) = s.stim.SOA_Fis;
        soaid(trlid == 9) = s.stim.SOA_2av;
        trlnsoa = [trlid, soaid];
        
        for b=1:s.blocks.N
            trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
            if s.pseudoran
                sTime = GetSecs;
                while badshuffle(trlnsoa, s.nallowrep)
                    trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
                    elapsedTime = GetSecs - sTime;
                    if elapsedTime > 10, error('Stack overflow while randomizing trial sequence.'); sca; end;
                end
            end
            data(1+(b-1)*ntrls:b*ntrls,3:4) = trlnsoa;
        end
    end
else
    
    %% mixed audiovisual, visual and auditory blocks
    if strcmp(s.cond, 'all') && ~s.randwblocks
        error('This code might actually work, but you have to include nil trials and pilot it again before using it for an experiment!')

%         ntrls   = s.ntrls + 1; % for 1 nil trial per block
%         trlid   = s.stim.trlid;
%         
%         data = zeros(ntrls*s.blocks.N, 6);
%         SOA  = s.stim.SOA;
%         id.a  = [4 7 4 7 4 7 4 7];
%         id.v  = [2 3 2 3 2 3 2 3 2 3 2 3];
%         id.av = [5 6 8 9 5 6 8 9 5 6 8 9];
%
%         if length(s.stim.trlid) ~= length(id.v)
%             error('id.v and id.av have to have the same number of items as s.stim.trlid')
%         end
%         
%         % randomize blocks and create corresponding trial sequences
%         cond = Shuffle([repmat('a', [s.blocks.Na  1]); ...
%                         repmat('v', [s.blocks.Nv  1]); ...
%                         repmat('m', [s.blocks.Nav 1])]);
%         
%         for b = 1:numel(cond)
%             % determine trial type of this block
%             switch floor(cond(b))
%                 case floor('a'), trlid = id.a;
%                 case floor('v'), trlid = id.v;
%                 case floor('m'), trlid = id.av;
%             end
%             % randomize trial sequence for block
%             trlid   = repmat(trlid, [1 ntrls/numel(trlid)])';
%             soaid   = repmat(SOA, [1 ntrls/numel(SOA)])';
%             trlnsoa = [trlid, sort(soaid)];
%             trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
%             if s.pseudoran
%                 sTime = GetSecs;
%                 while badshuffle(trlnsoa, s.nallowrep)
%                     trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
%                     elapsedTime = GetSecs - sTime;
%                     if elapsedTime > 10, error('Stack overflow while randomizing trial sequence.');  end;
%                 end
%             end
%             data(1+(b-1)*ntrls:b*ntrls,3:4) = trlnsoa;
%         end
        
    else
        %% or all trials randomized within the same blocks
        
%                         if strcmp(s.paradigm, 'yesno')
%                             ntrls   = s.ntrls + 1; % for 1 nil trial per block
%                         else
%                             ntrls   = s.ntrls;
%                         end
        trlid   = s.stim.trlid;
        
        data = zeros(ntrls*s.blocks.N, 6);
        SOA  = s.stim.SOA;
        
        nblocks = s.blocks.N;
        
        % randomize trl condition vector
        trlid   = repmat(trlid, [1 s.ntrls/numel(trlid)])';
        soaid   = repmat(SOA,   [1 s.ntrls/numel(SOA)])';
        trlnsoa = [trlid, sort(soaid)];
        
        if strcmp(s.paradigm, 'yesno') && ~strcmp(s.cond, 'multiAttAud') && ~strcmp(s.cond,'audio')
            % add average SOA over conditions for 1 Nil trial per block:
            avg_SOA = round(mean(soaid/(1/s.disp.refr))) * (1/s.disp.refr);
            trlid = [trlid; 1];
            soaid = [sort(soaid); avg_SOA];
            trlnsoa = [trlid, soaid];
        end
        
        for b=1:s.blocks.N
            trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
            if s.pseudoran
                sTime = GetSecs;
                %   This is a terribly inefficient way of doing this, but as I have a decent number
                %   of conditions it is unlikely that I will run into problems here.
                while badshuffle(trlnsoa, s.nallowrep)
                    trlnsoa = trlnsoa(randperm(size(trlnsoa,1)),:);
                    elapsedTime = GetSecs - sTime;
                    if elapsedTime > 10, error('Stack overflow while randomizing trial sequence.'); sca; end;
                end
            end
            data(1+(b-1)*ntrls:b*ntrls,3:4) = trlnsoa;
        end
    end
    
end

% Create final dataset array (easier to manipulate than matrix)
data = dataset({[data(:,1) data(:,2) data(:,3) data(:,4) data(:,5) data(:,6)], ...
                   'resp',   'RT',    'trlid',   'soa',  'badRT', 'nbadresp'});
               
% add confidence judgment column
data.conf = nan(size(data,1),1);

% add response key configuration (left 1, right 2 = A, left 2, right 1 = B)
data.rkeycfg = repmat({s.keyswitch}, [size(data,1), 1]);


% % Errors (unbalanced design)
% if mod(ntrls, numel(SOA))>0, error('ERROR: dfi_generate_data.m: number of trials per block is not a multiple of SOA.', UD); end;
% if mod(ntrls, numel(trlid))>0, error('ERROR: dfi_generate_data.m: number of trials per block is not a multiple of trial types.', UD); end;
% freqtable = tabulate(data.trlid);
% imbalance = (freqtable(1,2)/numel(SOA)) - round(freqtable(1,2)/numel(SOA));
% if imbalance, error('ERROR: dfi_generate_data.m: Trial types and SOAs are unbalanced.', UD); end;              


% select only the top part for practice sessions
if ispractice
    data(s.pract.ntrls+1:end,:) = [];
end


%
% ------------------ Nested functions --------------------
%
    function nogoodshuffle = badshuffle(mat, nrep)
        % goes through a vector of numerical items and checks if there are a number
        % of nrep repititions in succession. If there are it gives 1, otherwise 0.
        csum = 1;
        msum = 1;
        vect = sum(mat,2);
        if numel(vect) <= nrep
            nogoodshuffle = 0;
            return; 
        end;
        for i = 2:numel(vect)
            if vect(i) == vect(i-1)
                csum = csum + 1;
                if csum > msum
                    msum = csum;
                end
            else
                csum = 1;
            end
        end % efor
        if msum >= nrep
            nogoodshuffle = 1;
        else
            nogoodshuffle = 0;
        end
    end % efun

end
% eof














