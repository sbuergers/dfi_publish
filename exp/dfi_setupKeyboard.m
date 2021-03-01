function s = dfi_setupKeyboard(s)
% s = dfi_setupKeyboard(s)
% reads in the response keys and prepares and creates a KbQueue to save
% all future keystrokes with high accuracy in addition to KbCheck 
%
% ----------------------------------------
% adapted from Agoston Mihalik, 
% last updated, April 2015
%
    KbName('UnifyKeyNames'); % Match any key code scheme to mac os x
    % or use predefined keys
    s.key.code = KbName([s.key.resp, s.key.quit]);
%     if ~debug, ListenChar(2); end % listening, not writing 
    
    % Initialize KbQueue to retain all response-key presses
    keysOfInterest=zeros(1,256);
	keysOfInterest(KbName({[s.key.resp, s.key.quit]}))=1;
    KbQueueCreate(s.kb.subj, keysOfInterest);
    
end % end dfi_setupKeyboard