function PlotGroupLevel(AllData,DataPF,DataMLE)

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare the colours %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set some colors for the different modalities
colorAuditory         = [ 0  255  0 ] ./ 255;   %green
colorVisual           = [255  0   0 ] ./ 255;   %red
colorCongruent        = [ 0   0  255] ./ 255;   %blue
colorIncongruent_VLAR = [255  0  255] ./ 255;   %magenta
colorIncongruent_ALVR = [ 0  255 255] ./ 255;   %cyan
colorMLE              = [ 0   0   0 ] ./ 255;   %black
colorCodings = {colorAuditory,colorVisual,colorCongruent,colorIncongruent_VLAR,colorIncongruent_ALVR,colorMLE};
colorCodings2 = {'g', 'r', 'b', 'm', 'c', 'k'};

StimTypesNames = AllData{1,2}.StimTypesNames;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the psychometric functions (+/- SEM) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Find nr of Subjects
nSubjects = size(DataMLE,1);

nTimePoints = 1000;
standardTimeline = linspace(-2,2,nTimePoints);
PF = @PAL_CumulativeNormal;

%Collect the data per person
paramsFitted = nan(nSubjects,5,4);
Locations = nan(nSubjects,13);
deltaAVs = nan(nSubjects,1);
FittedCurvesEMP = nan(nSubjects,5,nTimePoints);
FittedCurvesMLE = nan(nSubjects,3,nTimePoints);
for i=1:nSubjects
    
    %Collect the relevant data for this subject
    paramsFitted(i,:,:) = DataPF{i,2}.fit.paramsFitted;
    Locations(i,:) = AllData{i,2}.Locations;
    deltaAVs(i,1) = DataPF{i,2}.data.deltaAV;
    
    %Create subj-specific timeline and PF curves
    LocationsFineGrain = linspace(-2*deltaAVs(i,1),2*deltaAVs(i,1),nTimePoints);
    for j=1:5
        FittedCurvesEMP(i,j,:) = PF(squeeze(paramsFitted(i,j,:)),LocationsFineGrain);
        
        %Compute bimodal curves from the MLE predicted parameters  
        if ismember(j,3:5)
            FittedCurvesMLE(i,j-2,:) = PF([DataMLE{i,2}.PSEsMLE_CongrAssumption(j-2,1) 1/(sqrt(2)*DataMLE{i,2}.SDavMLE(1,1)) DataMLE{i,2}.lapseR(1,1) DataMLE{i,2}.lapseR(1,1)], LocationsFineGrain);
        end
    end
end

%Plot unisensory versus bisensory congruent (Figure 2A)
figure; set(gcf, 'Position', get(0,'Screensize'));         

hData = [];
LegendTexts = cell(1);

hold on; box on;  
plot([-2 2],[.5 .5],'k--');
plot([0 0],[0 1],'k--');
for j=1:3
    FittedCurvesThisStimtype = squeeze(FittedCurvesEMP(:,j,:));
    meanCurve = mean(FittedCurvesThisStimtype,1);
    semCurve = std(FittedCurvesThisStimtype,[],1)./sqrt(nSubjects);
    [hl, hp] = boundedline(standardTimeline, meanCurve, semCurve, colorCodings2{j},'alpha');
    hData = [hData hl];
    LegendTexts{numel(hData)} = StimTypesNames{j}; 
end

%Plot the MLE-predicted psychometric curve (congruent)
meanMLECurves = squeeze(mean(FittedCurvesMLE,1));
hl = plot(standardTimeline,meanMLECurves(1,:),'--','Color','black','LineWidth',4); %congruent             
hData = [hData hl];
LegendTexts{numel(hData)} = 'MLE Predicted'; 

hold off;
xlim([-1.6 1.6]);
ylim([0 1]);
set(gca,'XTick',-2:0.5:2);
set(gca,'XTickLabel',{'-2*deltaAV','-1.5*deltaAV','-deltaAV','-0.5*deltaAV','0°','0.5*deltaAV','deltaAV','1.5*deltaAV','2*deltaAV'});
title('GroupLevel Normalised Psychometric Functions (+/-SEM)');
xlabel('Probe Location (relative to deltaAV)');
ylabel('P("probe right")');
leg = legend(hData,LegendTexts);
set(leg,'Location','northwest');

%Plot bisensory congruent and incongruent (Figure 2B)
figure; set(gcf, 'Position', get(0,'Screensize'));         

hData = [];
LegendTexts = cell(1);

hold on; box on;  
plot([-2 2],[.5 .5],'k--');
plot([-.5 -.5],[0 1],'k--');
plot([0 0],[0 1],'k--');
plot([.5 .5],[0 1],'k--');
for j=3:5
    FittedCurvesThisStimtype = squeeze(FittedCurvesEMP(:,j,:));
    meanCurve = mean(FittedCurvesThisStimtype,1);
    semCurve = std(FittedCurvesThisStimtype,[],1)./sqrt(nSubjects);
    [hl, hp] = boundedline(standardTimeline, meanCurve, semCurve, colorCodings2{j},'alpha');
    hData = [hData hl];
    LegendTexts{numel(hData)} = StimTypesNames{j}; 
end

%Plot the MLE-predicted psychometric curves
meanMLECurves = squeeze(mean(FittedCurvesMLE,1));     
hl4 = plot(standardTimeline,meanMLECurves(2,:),'-.','Color',colorCodings{4},'LineWidth',4); %incongr1 (VL AR)
hl5 = plot(standardTimeline(10:end),meanMLECurves(3,10:end),'--','Color',colorCodings{5},'LineWidth',4); %incongr2 (AL VR)
hData = [hData hl4 hl5];
LegendTexts{numel(hData)-1} = 'MLE Predicted VL AR'; 
LegendTexts{numel(hData)} = 'MLE Predicted VR AL';

hold off;
xlim([-1.6 1.6]);
ylim([0 1]);
set(gca,'XTick',-2:0.5:2);
set(gca,'XTickLabel',{'-2*deltaAV','-1.5*deltaAV','-deltaAV','-0.5*deltaAV','0°','0.5*deltaAV','deltaAV','1.5*deltaAV','2*deltaAV'});
title('GroupLevel Normalised Psychometric Functions (+/-SEM)');
xlabel('Probe Location (relative to deltaAV)');
ylabel('P("probe right")');
leg = legend(hData,LegendTexts);
set(leg,'Location','northwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare the other parameters %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Collect the SDs and auditory weights
SDaud        = nan(nSubjects,1);
SDvis        = nan(nSubjects,1);
SDcongr      = nan(nSubjects,1);
SDincongr1   = nan(nSubjects,1);
SDincongr2   = nan(nSubjects,1);
SDmle        = nan(nSubjects,1);

wAemp1noass  = nan(nSubjects,1);
wAemp2noass  = nan(nSubjects,1);
wAempwithass = nan(nSubjects,1);
wAmle        = nan(nSubjects,1);
for i=1:nSubjects
    SDaud(i)        = DataMLE{i,2}.SDsEmp(1,1);        
    SDvis(i)        = DataMLE{i,2}.SDsEmp(2,1);
    SDcongr(i)      = DataMLE{i,2}.SDsEmp(3,1);
    SDincongr1(i)   = DataMLE{i,2}.SDsEmp(4,1);
    SDincongr2(i)   = DataMLE{i,2}.SDsEmp(5,1);
    SDmle(i)        = DataMLE{i,2}.SDavMLE(1,1);
    
    wAemp1noass(i)  = DataMLE{i,2}.wAemp_NoAssumption_VLAR(1,1);
    wAemp2noass(i)  = DataMLE{i,2}.wAemp_NoAssumption_ALVR(1,1);
    wAempwithass(i) = DataMLE{i,2}.wAemp_WithAssumption_AVG(1,1);
    wAmle(i)        = DataMLE{i,2}.wAmle(1,1);
end

%Normalize the SDs such that the MLE-predicted SD is set to 1.
SDaud_norm       = SDaud./SDmle;
SDvis_norm       = SDvis./SDmle;
SDcongr_norm     = SDcongr./SDmle;
SDincongr1_norm  = SDincongr1./SDmle;
SDincongr2_norm  = SDincongr2./SDmle;
SDmle_norm       = SDmle./SDmle;    

%Compute the 95% CI (based on the SEMs) of the (normalized) SDs and auditory weights
SDaud_norm_CI      = mean(SDaud_norm) + [-1.96 1.96]*(std(SDaud_norm)/sqrt(nSubjects));
SDvis_norm_CI      = mean(SDvis_norm) + [-1.96 1.96]*(std(SDvis_norm)/sqrt(nSubjects));
SDcongr_norm_CI    = mean(SDcongr_norm) + [-1.96 1.96]*(std(SDcongr_norm)/sqrt(nSubjects));
SDincongr1_norm_CI = mean(SDincongr1_norm) + [-1.96 1.96]*(std(SDincongr1_norm)/sqrt(nSubjects));      
SDincongr2_norm_CI = mean(SDincongr2_norm) + [-1.96 1.96]*(std(SDincongr2_norm)/sqrt(nSubjects));     
SDmle_norm_CI      = mean(SDmle_norm) + [-1.96 1.96]*(std(SDmle_norm)/sqrt(nSubjects));

wAemp1noass_CI     = mean(wAemp1noass) + [-1.96 1.96]*(std(wAemp1noass)/sqrt(nSubjects));  
wAemp2noass_CI     = mean(wAemp2noass) + [-1.96 1.96]*(std(wAemp2noass)/sqrt(nSubjects)); 
wAempwithass_CI    = mean(wAempwithass) + [-1.96 1.96]*(std(wAempwithass)/sqrt(nSubjects)); 
wAmle_CI           = mean(wAmle) + [-1.96 1.96]*(std(wAmle)/sqrt(nSubjects));  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the SDs (Figure 2C) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Open a figure
figure; clf('reset'); set(gcf, 'Position', get(0,'Screensize'));        %open, clear and maximize figure
ax = gca; hold on;

yMax = 0;
yMin = 10;
hData = [];
LegendTexts = cell(1);
barWidth = 0.9;

LegendTextsFirstParts = {'Auditory: SD = '; 'Visual: SD = '; 'AVempCongr: SD = '; 'AVempVLAR: SD = '; 'AVempVRAL: SD = '; 'AVmle: SD = '};
SDdataSets = {SDaud_norm; SDvis_norm; SDcongr_norm; SDincongr1_norm; SDincongr2_norm; SDmle_norm};
SDdataCI = {SDaud_norm_CI; SDvis_norm_CI; SDcongr_norm_CI; SDincongr1_norm_CI; SDincongr2_norm_CI; SDmle_norm_CI};

for k=1:5
    %Plot the emprical bars
    SD_temp = round(mean(SDdataSets{k,1})*100)/100;
    x_temp_bar = k + barWidth*[-0.5 0.5 0.5 -0.5];           %[L R R L]; 
    y_temp_bar = [SD_temp SD_temp 0 0];                      %[T T B B];
    patch(x_temp_bar,y_temp_bar,1,'Parent',ax,'EdgeColor',colorCodings{k},'LineStyle','-','Linewidth',1.5,'FaceColor',colorCodings{k},'FaceAlpha',0.2); 
    yMin = min(yMin,SD_temp);

    %Plot the 95% CI on top
    CI_temp = SDdataCI{k,1};
    plot([k k],[CI_temp(1) CI_temp(2)],'k-','Linewidth',2);
    plot(k+barWidth*[-0.25 0.25],[CI_temp(2) CI_temp(2)],'k-','Linewidth',2);
    plot(k+barWidth*[-0.25 0.25],[CI_temp(1) CI_temp(1)],'k-','Linewidth',2);  
    yMin = min(yMin,CI_temp(1));
    yMax = max(yMax,CI_temp(2));
        
    %Plot dummies for legend and create legend text
    dummyHandle = plot([100 100],[-1 -2],'-','Color',colorCodings{k},'LineWidth',2);
    hData = [hData dummyHandle];
    LegendTexts{numel(hData)} = [LegendTextsFirstParts{k,1} num2str(SD_temp) '°'];              
end

%Plot the mle-predicted bar
SD_temp = round(mean(SDdataSets{(k+1),1})*100)/100;
x_temp_bar = (k+1) + barWidth*[-0.5 0.5 0.5 -0.5];           %[L R R L]; 
y_temp_bar = [SD_temp SD_temp 0 0];                      %[T T B B];
patch(x_temp_bar,y_temp_bar,1,'Parent',ax,'EdgeColor',colorCodings{3},'LineStyle',':','Linewidth',1.5,'FaceColor',colorCodings{3},'FaceAlpha',0.05); 
yMin = min(yMin,SD_temp);

%Plot dummies for MLE-predicted legend and create legend text
dummyHandle = plot([100 100],[-1 -2],':','Color',colorCodings{3},'LineWidth',2);
hData = [hData dummyHandle];
LegendTexts{numel(hData)} = [LegendTextsFirstParts{(k+1),1} num2str(SD_temp) '°'];  

%Finish plot
ylabel('Normalised SD (°)');
%xlabel('Conditions')
set(ax,'XTick',1:(k+1));
set(ax,'XTickLabel',{'Aud','Vis','AVempCongr','AVempVLAR','AVempVRAL','AVmle'});
xlim([0.5 (k+1)+0.5]);
ylim([max(0,yMin-0.5) yMax+0.5]);     
title('GroupLevel Normalised SDs and 95% CI (+/-1.96*SEM)');
leg = legend(ax,hData,LegendTexts);
set(leg,'Location','eastoutside');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the auditory weights (Figure 2D) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Create a colormap going from magenta (0) to cyan (1) via white (0.5)
colormap1 = linspace(colorIncongruent_VLAR(1),colorIncongruent_ALVR(1),100)';
colormap1 = [colormap1 linspace(colorIncongruent_VLAR(2),colorIncongruent_ALVR(2),100)'];
colormap1 = [colormap1 linspace(colorIncongruent_VLAR(3),colorIncongruent_ALVR(3),100)'];

%Open a figure for the auditory weights barplot
figure; clf('reset'); set(gcf, 'Position', get(0,'Screensize'));        
ax = gca; hold on;

yMax = 0;
hData = [];
LegendTexts = cell(1);
barWidth = 0.9;

%NO ASSUMPTION
%Incongr1 VLAR
Aweight_Emp = round(mean(wAemp1noass)*1000)/1000;
x_temp_bar = 1 + barWidth*[-0.5 0.5 0.5 -0.5];                %[L R R L]; 
y_temp_bar = [Aweight_Emp Aweight_Emp 0 0];                   %[T T B B];
colors_bar = [0 0 0 0]';                                  %[LT, RT, RB, LB];
patch(x_temp_bar,y_temp_bar,colors_bar,'Parent',ax,'EdgeColor',colorIncongruent_VLAR,'LineStyle','-','Linewidth',1.5,'FaceColor','interp','FaceAlpha',0.2); 

CI_temp = wAemp1noass_CI;
plot([1 1],[CI_temp(1) CI_temp(2)],'k-','Linewidth',2);
plot(1+barWidth*[-0.25 0.25],[CI_temp(2) CI_temp(2)],'k-','Linewidth',2);
plot(1+barWidth*[-0.25 0.25],[CI_temp(1) CI_temp(1)],'k-','Linewidth',2);  
yMax = max(yMax,CI_temp(2));

dummyHandle = plot([100 100],[-1 -2], '-','LineWidth',2,'Color',colorIncongruent_VLAR);                  
hData = [hData dummyHandle];
LegendTexts{numel(hData)} = ['Emperical Auditory Weight VLAR NoAss = ' num2str(Aweight_Emp)];           

%Incongr1 VRAL
Aweight_Emp = round(mean(wAemp2noass)*1000)/1000;
x_temp_bar = 2 + barWidth*[-0.5 0.5 0.5 -0.5];                %[L R R L]; 
y_temp_bar = [Aweight_Emp Aweight_Emp 0 0];                   %[T T B B];
colors_bar = [.99 .99 .99 .99]';                                  %[LT, RT, RB, LB];
patch(x_temp_bar,y_temp_bar,colors_bar,'Parent',ax,'EdgeColor',colorIncongruent_ALVR,'LineStyle','-','Linewidth',1.5,'FaceColor','interp','FaceAlpha',0.2); 

CI_temp = wAemp2noass_CI;
plot([2 2],[CI_temp(1) CI_temp(2)],'k-','Linewidth',2);
plot(2+barWidth*[-0.25 0.25],[CI_temp(2) CI_temp(2)],'k-','Linewidth',2);
plot(2+barWidth*[-0.25 0.25],[CI_temp(1) CI_temp(1)],'k-','Linewidth',2);  
yMax = max(yMax,CI_temp(2));

dummyHandle = plot([100 100],[-1 -2], '-','LineWidth',2,'Color',colorIncongruent_ALVR);                  
hData = [hData dummyHandle];
LegendTexts{numel(hData)} = ['Emperical Auditory Weight VRAL NoAss= ' num2str(Aweight_Emp)];           

%WITH ASSUMPTION
%Plot the emperical bar
Aweight_Emp = round(mean(wAempwithass)*1000)/1000;
x_temp_bar = 3 + barWidth*[-0.5 0.5 0.5 -0.5];                %[L R R L]; 
y_temp_bar = [Aweight_Emp Aweight_Emp 0 0];                   %[T T B B];
colors_bar = [1 0.5 0 0.5]';                                  %[LT, RT, RB, LB];
patch(x_temp_bar,y_temp_bar,colors_bar,'Parent',ax,'LineStyle','-','Linewidth',1.5,'FaceColor','interp','FaceAlpha',0.5); %'EdgeColor','interp', --> doesn't work together with 'LineStyle'

%Plot the 95% CI on top
CI_temp = wAempwithass_CI;
plot([3 3],[CI_temp(1) CI_temp(2)],'k-','Linewidth',2);
plot(3+barWidth*[-0.25 0.25],[CI_temp(2) CI_temp(2)],'k-','Linewidth',2);
plot(3+barWidth*[-0.25 0.25],[CI_temp(1) CI_temp(1)],'k-','Linewidth',2);  
yMax = max(yMax,CI_temp(2));

%Plot dummy for legend and create legend text
dummyHandle = plot([100 100],[-1 -2], 'k-','LineWidth',2);                  %interpolated patches are not possible in Legend
hData = [hData dummyHandle];
LegendTexts{numel(hData)} = ['Emperical Auditory Weight AV WithAss = ' num2str(Aweight_Emp)];           

%Plot the MLE-predicted bar
Aweight_MLE = round(mean(wAmle)*1000)/1000;
x_temp_bar = 4 + barWidth*[-0.5 0.5 0.5 -0.5];                  %[L R R L]; 
y_temp_bar = [Aweight_MLE Aweight_MLE 0 0];                     %[T T B B];
colors_bar = [1 0.5 0 0.5]';                                    %[LT, RT, RB, LB];  
patch(x_temp_bar,y_temp_bar,colors_bar,'Parent',ax,'LineStyle',':','Linewidth',1.5,'FaceColor','interp','FaceAlpha',0.2); %'EdgeColor','interp', --> doesn't work together with 'LineStyle'

%Plot the 95% CI on top (SEM would be too small to see anything)
CI_temp = wAmle_CI;
plot([4 4],[CI_temp(1) CI_temp(2)],'k-','Linewidth',2);
plot(4+barWidth*[-0.25 0.25],[CI_temp(2) CI_temp(2)],'k-','Linewidth',2);
plot(4+barWidth*[-0.25 0.25],[CI_temp(1) CI_temp(1)],'k-','Linewidth',2);  
yMax = max(yMax,CI_temp(2));

%Plot dummy for legend and create legend text
dummyHandle = plot([100 100],[-1 -2], 'k:','LineWidth',2);                  %interpolated patches are not possible in Legend
hData = [hData dummyHandle];
LegendTexts{numel(hData)} = ['MLE-predicted Auditory Weight = ' num2str(Aweight_MLE)];  

%Finish plot
hold off;
colormap(ax,colormap1)
ylabel('Auditory Weight');
%xlabel('Conditions')
set(ax,'XTick',1:4);
set(ax,'XTickLabel',{'Emp VLAR No-ass', 'Emp VRAL No-ass', 'Empirical Congr-ass','MLE-predicted'});
xlim([0.5 4.5]);
ylim([0 yMax+0.05]);
title('Auditory Weights and 95% CI (+/-1.96*SEM)');
leg = legend(ax,hData,LegendTexts);
set(leg,'Location','eastoutside');

return %[EOF]