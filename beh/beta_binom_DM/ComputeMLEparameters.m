function dataMLE = ComputeMLEparameters(dataPF)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract the relevant info from the fits %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%AV disparity used for this subject
deltaAV = dataPF.data.deltaAV;

%The true emperical PSE and SD from the 'dependent' fit
PSEsEmp = dataPF.fit.paramsFitted(:,1);        
SDsEmp = (1./dataPF.fit.paramsFitted(:,2))/sqrt(2);                         %We use the true SDs, not the JNDs (=1/beta)
lapseR = dataPF.fit.paramsFitted(1,4);
eta = dataPF.fit.eta;

%Number of boostrap samples
BootstrapN = dataPF.settings.BootstrapN;

%'paramsSim' is a [BootstrapN x nStimTypes x 4] matrix
paramsSim = dataPF.Bootstrap.paramsFitted;
PSEsBootstrap = squeeze(paramsSim(:,:,1));
SDsBootstrap = squeeze((1./paramsSim(:,:,2))/sqrt(2));           
lapseRBootstrap = squeeze(paramsSim(:,1,4));
etaBootstrap = dataPF.Bootstrap.eta;

%Merge the emperical and bootstrapped parameters such that matrix dimensions are [nStimTypes x BootstrapN]
PSEsEmp = [PSEsEmp PSEsBootstrap'];          
SDsEmp  = [SDsEmp  SDsBootstrap'];                                          %Note that the true values are the first column, and the bootstrapped values come afterwards
lapseR = [lapseR lapseRBootstrap'];
eta = [eta etaBootstrap'];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute the MLE parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Empirical Auditory weights

%1. Make no assumption: the emperical weights are based upon the unisensory PSEs, including their respective LR-bias
wAemp_NoAssumption_VLAR  = (PSEsEmp(4,:) - (PSEsEmp(2,:) + (deltaAV/2))) ./ ((PSEsEmp(1,:) - (deltaAV/2)) - (PSEsEmp(2,:) + (deltaAV/2)));       
wAemp_NoAssumption_ALVR  = (PSEsEmp(5,:) - (PSEsEmp(2,:) - (deltaAV/2))) ./ ((PSEsEmp(1,:) + (deltaAV/2)) - (PSEsEmp(2,:) - (deltaAV/2)));
wAemp_NoAssumption_AVG = mean([wAemp_NoAssumption_VLAR; wAemp_NoAssumption_ALVR]); %average of the two above

%2. Assume that the unisensory PSEs contain noise, and use a fixed value for both unisensory PSEs instead (either PSE_congr as in Fetsch et al., 2012, or zero as in Alais and Burr, 2004)
wAemp_WithAssumption_VLAR = (PSEsEmp(4,:) - (PSEsEmp(3,:) + (deltaAV/2))) ./ (-deltaAV);
wAemp_WithAssumption_ALVR = (PSEsEmp(5,:) - (PSEsEmp(3,:) - (deltaAV/2))) ./ (deltaAV);
wAemp_WithAssumption_AVG = 0.5 + (PSEsEmp(5,:) - PSEsEmp(4,:)) ./ (2*deltaAV);  %the above simplifies, because the fixed values for the unisensory PSEs are subtracted away. 

%MLE predictions (based on unisensory SDs)

%1. Auditory weights 
wAmle = (SDsEmp(2,:).^2) ./ (SDsEmp(1,:).^2 + SDsEmp(2,:).^2);          
wVmle = (SDsEmp(1,:).^2) ./ (SDsEmp(1,:).^2 + SDsEmp(2,:).^2);          

%2. Bisensory SDs
SDavMLE = sqrt( (SDsEmp(1,:).^2 .* SDsEmp(2,:).^2) ./ (SDsEmp(1,:).^2 + SDsEmp(2,:).^2) );           

%3. Bisensory PSEs
%Note, the modality-specific PSE is predicted to be on the opposite side of the respective modality's true location (because we present the AV conflict in the probe, not in the standard)    
%3.1 Make no assumption: the MLE predicted PSE for bimodal trials is based upon the unimodal PSEs, including their respective LR-bias
PSEsMLE_NoAssumption = nan(3,BootstrapN+1);                          
PSEsMLE_NoAssumption(1,:) = wVmle.*PSEsEmp(2,:)               + wAmle.*PSEsEmp(1,:);                     %Congruent
PSEsMLE_NoAssumption(2,:) = wVmle.*(PSEsEmp(2,:)+(deltaAV/2)) + wAmle.*(PSEsEmp(1,:)-(deltaAV/2));       %Incongruent 1 (V left, A right)
PSEsMLE_NoAssumption(3,:) = wVmle.*(PSEsEmp(2,:)-(deltaAV/2)) + wAmle.*(PSEsEmp(1,:)+(deltaAV/2));       %Incongruent 2 (V right, A left)
%3.2. Assume that the PSEs of the congruent bimodal responses are more accurate estimates than the PSEs of the unimodal responses, as in Fetsch et al., 2012
PSEsMLE_CongrAssumption = nan(3,BootstrapN+1); 
PSEsMLE_CongrAssumption(1,:) = PSEsEmp(3,:);                                                             %Congruent
PSEsMLE_CongrAssumption(2,:) = wVmle.*(PSEsEmp(3,:)+(deltaAV/2)) + wAmle.*(PSEsEmp(3,:)-(deltaAV/2));    %Incongruent 1 (V left + A right)
PSEsMLE_CongrAssumption(3,:) = wVmle.*(PSEsEmp(3,:)-(deltaAV/2)) + wAmle.*(PSEsEmp(3,:)+(deltaAV/2));    %Incongruent 2 (V right + A left)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
%%% Perform some bootstrap tests on the SDs, PSEs, and weights %%%                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define all the comparisons for the bootstrap t-tests                       Last column defines one or two-sided tests: 1st column values are hypothesized to be larger    
ComparisonsToDo = {SDsEmp(1,:),             SDsEmp(2,:),                    '1. SD Aud > OR < than SD Vis',                                                    2 ; ... 
                   mean(SDsEmp(4:5,:)),     SDsEmp(3,:),                    '2. SD Incongruent > SD Congruent',                                                1 ; ...
                   SDsEmp(1,:),             SDsEmp(3,:),                    '3. SD Auditory > SD Congruent',                                                   1 ; ...
                   SDsEmp(2,:),             SDsEmp(3,:),                    '4. SD Visual > SD Congruent',                                                     1 ; ...     
                   SDsEmp(3,:),             SDavMLE(1,:),                   '5. SD Congruent > SD MLE',                                                        1 ; ...
                   mean(SDsEmp(4:5,:)),     SDavMLE(1,:),                   '6. SD Incongruent > SD MLE',                                                      1 ; ...
                   wAmle(1,:),              wAemp_NoAssumption_AVG(1,:),    '7. Auditory Weight Emperical_NoAssumption < Auditory Weight MLE-predicted',       1 ; ...
                   wAmle(1,:),              wAemp_WithAssumption_AVG(1,:),  '8. Auditory Weight Emperical_WithAssumption < Auditory Weight MLE-predicted',     1 };
                     %Larger                  %Smaller     
               
%Perform all the bootstrap ttests following the procedure described in chapter 16 of Efron's and Tibshirani's An Introduction to the bootstrap (page 220-224)
nComparisons = size(ComparisonsToDo,1);
BootstrapComparisonResults = cell(nComparisons,2);
for i=1:nComparisons
    
    %Data (bootstrapped SDs, PSEs or weights)
    x_original = ComparisonsToDo{i,1}(1,1);
    y_original = ComparisonsToDo{i,2}(1,1);
    
    x_bootstraps = ComparisonsToDo{i,1}(1,2:(BootstrapN+1));
    y_bootstraps = ComparisonsToDo{i,2}(1,2:(BootstrapN+1));    
    
    %Compute the differences between the parameters     
    XminusY_original = x_original - y_original;             %scalar   
    XminusY_bootstraps = x_bootstraps - y_bootstraps;       %vector of length n_bootstraps
    
    %Assuming that the original difference is the best estimate of the mean of the bootstrapped differences, 
    %then we can create a null distribution of differences (centred at 0: no difference) by subtracting the original difference from the bootstrapped differences
    nullDistributionOfDifferences = XminusY_bootstraps - XminusY_original;   
    
    %Find the proportion of bootstrapped differences under H0, that are more extreme than the original difference of interest (in either direction: 2 sided, or in one direction only: 1-sided) 
    if ComparisonsToDo{i,4} == 2
        bootstrapP = (1 + sum(abs(nullDistributionOfDifferences) > abs(XminusY_original))) / (BootstrapN+1);        % Two sided test   
    elseif ComparisonsToDo{i,4} == 1
        bootstrapP = (1 + sum(nullDistributionOfDifferences > XminusY_original)) / (BootstrapN+1);                  % One sided test (Note the order of "ComparisonsToDo")
    end
    %Save and collect the results
    BootstrapComparisonResults{i,1} = bootstrapP;           %p-value of bootstrap-test of differences
end
BootstrapComparisonResults(:,2) = ComparisonsToDo(:,3);     %What was tested in words

%%%%%%%%%%%%%%%%%%%%%%%                 
%%% End of function %%%                 
%%%%%%%%%%%%%%%%%%%%%%%                 

%Collect all variables in struct 'dataMLE' for the output of this function
dataMLE.deltaAV = deltaAV;
dataMLE.BootstrapN = BootstrapN;

dataMLE.PSEsEmp = PSEsEmp;      
dataMLE.SDsEmp  = SDsEmp;
dataMLE.lapseR  = lapseR;
dataMLE.eta     = eta;

dataMLE.wAemp_NoAssumption_VLAR   = wAemp_NoAssumption_VLAR;
dataMLE.wAemp_NoAssumption_ALVR   = wAemp_NoAssumption_ALVR;
dataMLE.wAemp_NoAssumption_AVG    = wAemp_NoAssumption_AVG;
dataMLE.wAemp_WithAssumption_VLAR = wAemp_WithAssumption_VLAR;
dataMLE.wAemp_WithAssumption_ALVR = wAemp_WithAssumption_ALVR;
dataMLE.wAemp_WithAssumption_AVG  = wAemp_WithAssumption_AVG;

dataMLE.wAmle   = wAmle;
dataMLE.SDavMLE = SDavMLE;

dataMLE.PSEsMLE_NoAssumption = PSEsMLE_NoAssumption;
dataMLE.PSEsMLE_CongrAssumption = PSEsMLE_CongrAssumption;

dataMLE.BootstrapComparisonResults = BootstrapComparisonResults;

return %[EOF]
