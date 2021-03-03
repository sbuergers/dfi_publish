function PlotIndividualResults(DataMLE)

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Collect all data %%%
%%%%%%%%%%%%%%%%%%%%%%%%

%Find nr of Subjects
nSubjects = size(DataMLE,1);

%Initialize
Aud_SDs = nan(nSubjects,1);
Vis_SDs = nan(nSubjects,1);
AV_SDs_Emp = nan(nSubjects,1);
AV_SDs_MLE = nan(nSubjects,1);
BestUni_SDs = nan(nSubjects,1);
BestUni_AudOrVis = nan(nSubjects,1);

Aud_Weight_Emp = nan(nSubjects,1);
Aud_Weight_MLE = nan(nSubjects,1);

Test_SDaud_vs_SDvis = nan(nSubjects,1);
Test_SDunibest_vs_SDav = nan(nSubjects,1);          %Bootstrapped tests
Test_SDavEMP_vs_SDavMLE = nan(nSubjects,1);
Test_AweightEMP_vs_AweightMLE = nan(nSubjects,1);

CI_SDaud = nan(nSubjects,2);
CI_SDvis = nan(nSubjects,2);
CI_SDunibest = nan(nSubjects,2);                    %Confidence intervals (low - high)
CI_SDavEMP = nan(nSubjects,2);
CI_SDavMLE = nan(nSubjects,2);
CI_AweightEMP = nan(nSubjects,2);
CI_AweightMLE = nan(nSubjects,2);

%Collect the data per person
for i=1:nSubjects
    
    Aud_SDs(i,1) = DataMLE{i,2}.SDsEmp(1,1);
    Vis_SDs(i,1) = DataMLE{i,2}.SDsEmp(2,1);
    AV_SDs_Emp(i,1) = DataMLE{i,2}.SDsEmp(3,1);
    AV_SDs_MLE(i,1) = DataMLE{i,2}.SDavMLE(1,1);

    Aud_Weight_Emp(i,1) = DataMLE{i,2}.wAemp_WithAssumption_AVG(1,1);
    Aud_Weight_MLE(i,1) = DataMLE{i,2}.wAmle(1,1);
    
    if Aud_SDs(i,1) < Vis_SDs(i,1)
        BestUni_AudOrVis(i,1)           = 1;
        Test_SDunibest_vs_SDav(i,1)     = DataMLE{i,2}.BootstrapComparisonResults{3,1};     %AudSD_vs_avCongrSD
        BestUni_SDs(i,1)                = Aud_SDs(i,1);
        CI_SDunibest(i,1)               = prctile(DataMLE{i,2}.SDsEmp(1,2:end),2.5);
        CI_SDunibest(i,2)               = prctile(DataMLE{i,2}.SDsEmp(1,2:end),97.5);
    else
        BestUni_AudOrVis(i,1)           = 2;
        Test_SDunibest_vs_SDav(i,1)     = DataMLE{i,2}.BootstrapComparisonResults{4,1};     %VisSD_vs_avCongrSD
        BestUni_SDs(i,1)                = Vis_SDs(i,1);
        CI_SDunibest(i,1)               = prctile(DataMLE{i,2}.SDsEmp(2,2:end),2.5);
        CI_SDunibest(i,2)               = prctile(DataMLE{i,2}.SDsEmp(2,2:end),97.5);
    end
    Test_SDaud_vs_SDvis(i,1)            = DataMLE{i,2}.BootstrapComparisonResults{1,1};     %AudSD_vs_VisSD
    Test_SDavEMP_vs_SDavMLE(i,1)        = DataMLE{i,2}.BootstrapComparisonResults{5,1};     %AVempCongrSD_vs_AVmleSD
    Test_AweightEMP_vs_AweightMLE(i,1)  = DataMLE{i,2}.BootstrapComparisonResults{8,1};     %AweightEMP_vs_AweightMLE
    
    CI_SDaud(i,:)      = [prctile(DataMLE{i,2}.SDsEmp(1,2:end),2.5),                    prctile(DataMLE{i,2}.SDsEmp(1,2:end),97.5)];
    CI_SDvis(i,:)      = [prctile(DataMLE{i,2}.SDsEmp(2,2:end),2.5),                    prctile(DataMLE{i,2}.SDsEmp(2,2:end),97.5)];
    CI_SDavEMP(i,:)    = [prctile(DataMLE{i,2}.SDsEmp(3,2:end),2.5),                    prctile(DataMLE{i,2}.SDsEmp(3,2:end),97.5)];
    CI_SDavMLE(i,:)    = [prctile(DataMLE{i,2}.SDavMLE(1,2:end),2.5),                   prctile(DataMLE{i,2}.SDavMLE(1,2:end),97.5)];
    CI_AweightEMP(i,:) = [prctile(DataMLE{i,2}.wAemp_WithAssumption_AVG(1,2:end),2.5),  prctile(DataMLE{i,2}.wAemp_WithAssumption_AVG(1,2:end),97.5)];
    CI_AweightMLE(i,:) = [prctile(DataMLE{i,2}.wAmle(1,2:end),2.5),                     prctile(DataMLE{i,2}.wAmle(1,2:end),97.5)];
end

%Define colours
colorEdges = {[0 1 0],[1 0 0],[0 0 1],[1 0 1],[0 1 1]};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the Two Unisensory SDs against each other %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the optimal participants for when AV_SD_Emp was smaller than the smallest unisensory modality    
figure; set(gcf, 'Position', get(0,'Screensize'));         %open and maximize figure
min_SD = min(min(CI_SDaud(:,1)),min(CI_SDvis(:,1)));
max_SD = max(max(CI_SDaud(:,2)),max(CI_SDvis(:,2)));
hold on; box on;
h1 = plot(-1,-1,'kd','MarkerFaceColor','w','MarkerSize',10); %dummy for legend
h2 = plot(-1,-1,'kd','MarkerFaceColor','k','MarkerSize',10); %dummy for legend
plot([min_SD max_SD],[min_SD max_SD],'k--');
%Plot two lines that mark when the variances differ too much from one another   
plot([min_SD max_SD],sqrt(2)*[min_SD max_SD],'b:');
plot(sqrt(2)*[min_SD max_SD],[min_SD max_SD],'b:');
%First plot all confidence intervals (so that the markers are always on top)    
for j=1:nSubjects
    plot([CI_SDaud(j,1) CI_SDaud(j,2)], [Vis_SDs(j,1) Vis_SDs(j,1)],'k-');      %Aud at x-axis
    plot([Aud_SDs(j,1) Aud_SDs(j,1)],   [CI_SDvis(j,1) CI_SDvis(j,2)],'k-');    %Vis at y-axis
end
%Then plot the markers (grey filled is suboptimal)
for j=1:nSubjects
    xyPoint = [Aud_SDs(j,1),Vis_SDs(j,1)];
    diamondPoints = 0.04*[-1.5 0; 0 -2.5; 1.5 0; 0 2.5]+repmat(xyPoint,[4 1]);  %LTRB
    if Test_SDaud_vs_SDvis(j,1) < 0.05
        plot(xyPoint(1),xyPoint(2),'d','MarkerSize',18,'MarkerEdgeColor',[.55 .55 .55],'MarkerFaceColor',[.55 .55 .55]); 
        plot(diamondPoints([4 1 2],1),diamondPoints([4 1 2],2),'Color',colorEdges{2},'LineWidth',2); %left
        plot(diamondPoints([4 3 2],1),diamondPoints([4 3 2],2),'Color',colorEdges{1},'LineWidth',2); %right
    else
        plot(xyPoint(1),xyPoint(2),'d','MarkerSize',18,'MarkerEdgeColor','w','MarkerFaceColor','w');
        plot(diamondPoints([4 1 2],1),diamondPoints([4 1 2],2),'Color',colorEdges{2},'LineWidth',2); %left
        plot(diamondPoints([4 3 2],1),diamondPoints([4 3 2],2),'Color',colorEdges{1},'LineWidth',2); %right
    end
end
hold off; xlim([min_SD max_SD]); ylim([min_SD max_SD]);
title(['Aud SD vs. Vis SD (N = ' num2str(nSubjects) ')']);
xlabel('Aud SD')
ylabel('Vis SD')
legend([h2 h1],{'Significant Difference','No Significant Difference'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the Best Unisensory vs Bisensory SDs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the optimal participants for when AV_SD_Emp was smaller than the smallest unisensory modality    
figure; set(gcf, 'Position', get(0,'Screensize'));         %open and maximize figure
min_SD = min(min(CI_SDunibest(:,1)),min(CI_SDavEMP(:,1)));
max_SD = max(max(CI_SDunibest(:,2)),max(CI_SDavEMP(:,2)));
hold on; box on;
h1 = plot(-1,-1,'kd','MarkerFaceColor','w','MarkerSize',10); %dummy for legend
h2 = plot(-1,-1,'kd','MarkerFaceColor','k','MarkerSize',10); %dummy for legend
plot([min_SD max_SD],[min_SD max_SD],'k--');
%First plot all confidence intervals (so that the markers are always on top)    
for j=1:nSubjects
    plot([CI_SDunibest(j,1) CI_SDunibest(j,2)], [AV_SDs_Emp(j,1) AV_SDs_Emp(j,1)],'k-');    %Uni at x-axis
    plot([BestUni_SDs(j,1) BestUni_SDs(j,1)],   [CI_SDavEMP(j,1) CI_SDavEMP(j,2)],'k-');    %AV at y-axis
end
%Then plot the markers (grey filled is suboptimal)
for j=1:nSubjects
    if Test_SDunibest_vs_SDav(j,1) < 0.05
        plot(BestUni_SDs(j,1),AV_SDs_Emp(j,1),'d','MarkerSize',15,'MarkerEdgeColor',colorEdges{BestUni_AudOrVis(j,1)},'MarkerFaceColor',[.55 .55 .55],'Linewidth',2);  
    else
        %faceColorTemp = interp1([0 1], [1 1 1; colorEdges{BestUni_AudOrVis(j,1)}], 0.2);    %I got this from "boundedline.m" - last value is the "transparency" --> use instead of 'w'
        plot(BestUni_SDs(j,1),AV_SDs_Emp(j,1),'d','MarkerSize',15,'MarkerEdgeColor',colorEdges{BestUni_AudOrVis(j,1)},'MarkerFaceColor','w','Linewidth',2);
    end
end
hold off; xlim([min_SD max_SD]); ylim([min_SD max_SD]);
title(['Best unisensory SD vs. empirical AV SD (N = ' num2str(nSubjects) ')']);
xlabel('Best unisensory SD')
ylabel('AV SD (empirical)')
legend([h2 h1],{'Significant Variance Reduction','Suboptimal'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the Emprical vs. MLE-predicted SDs %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the optimal participants for when AV_SD_Emp was similar to the AV_SD_MLE 
figure; set(gcf, 'Position', get(0,'Screensize'));         %open and maximize figure
min_SD = min(min(CI_SDavMLE(:,1)),min(CI_SDavEMP(:,1)));
max_SD = max(max(CI_SDavMLE(:,2)),max(CI_SDavEMP(:,2)));
hold on; box on;
h1 = plot(-1,-1,'kd','MarkerFaceColor','w','MarkerSize',10); %dummy for legend
h2 = plot(-1,-1,'kd','MarkerFaceColor','k','MarkerSize',10); %dummy for legend
plot([min_SD max_SD],[min_SD max_SD],'k--');
%First plot all confidence intervals (so that the markers are always on top)  
for j=1:nSubjects
    plot([CI_SDavMLE(j,1) CI_SDavMLE(j,2)], [AV_SDs_Emp(j,1) AV_SDs_Emp(j,1)],'k-');    %MLE AV SD at x-axis
    plot([AV_SDs_MLE(j,1) AV_SDs_MLE(j,1)], [CI_SDavEMP(j,1) CI_SDavEMP(j,2)],'k-');    %EMP AV SD at y-axis
end
%Then plot the markers (grey filled is suboptimal)
for j=1:nSubjects
    if Test_SDavEMP_vs_SDavMLE(j,1) < 0.05
        plot(AV_SDs_MLE(j,1),AV_SDs_Emp(j,1),'d','MarkerSize',15,'MarkerEdgeColor',colorEdges{3},'MarkerFaceColor',[.55 .55 .55],'Linewidth',2);  
    else
        %faceColorTemp = interp1([0 1], [1 1 1; colorEdges{BestUni_AudOrVis(j,1)}], 0.2);    %I got this from "boundedline.m" - last value is the "transparency" --> use instead of 'w'
        plot(AV_SDs_MLE(j,1),AV_SDs_Emp(j,1),'d','MarkerSize',15,'MarkerEdgeColor',colorEdges{3},'MarkerFaceColor','w','Linewidth',2);          
    end
end
hold off; xlim([min_SD max_SD]); ylim([min_SD max_SD]);
title(['Empirical vs. MLE-predicted AV SDs (N = ' num2str(nSubjects) ')']);
xlabel('AV SD (mle)')
ylabel('AV SD (empirical)')
legend([h2 h1],{'Significant: EMP > MLE','(Super-) Optimal'},'Location','northwest')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Plot the Emprical vs. MLE-predicted Auditory Weights %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Plot the optimal participants for when Aud_weight_Emp was was similar to the Aud_weight_MLE    
figure; set(gcf, 'Position', get(0,'Screensize'));         %open and maximize figure
min_Weight = min(min(CI_AweightMLE(:,1)),min(CI_AweightEMP(:,1)));
max_Weight = max(max(CI_AweightMLE(:,2)),max(CI_AweightEMP(:,2)));
hold on; box on;
h1 = plot(-1,-1,'kd','MarkerFaceColor','w','MarkerSize',10); %dummy for legend
h2 = plot(-1,-1,'kd','MarkerFaceColor','k','MarkerSize',10); %dummy for legend
plot([0 max_Weight],[0 max_Weight],'k--');
%First plot all confidence intervals (so that the markers are always on top)  
for j=1:nSubjects
    plot([CI_AweightMLE(j,1) CI_AweightMLE(j,2)], [Aud_Weight_Emp(j,1) Aud_Weight_Emp(j,1)],'k-');    %Uni at x-axis
    plot([Aud_Weight_MLE(j,1) Aud_Weight_MLE(j,1)],   [CI_AweightEMP(j,1) CI_AweightEMP(j,2)],'k-');    %AV at y-axis
end
%Then plot the markers (grey filled is suboptimal)
for j=1:nSubjects
    xyPoint = [Aud_Weight_MLE(j,1),Aud_Weight_Emp(j,1)];
    diamondPoints = 0.007*[-1.5 0; 0 -2.5; 1.5 0; 0 2.5]+repmat(xyPoint,[4 1]);  %LTRB
    if Test_AweightEMP_vs_AweightMLE(j,1) < 0.05
        plot(xyPoint(1),xyPoint(2),'d','MarkerSize',18,'MarkerEdgeColor',[.55 .55 .55],'MarkerFaceColor',[.55 .55 .55]); 
        plot(diamondPoints([4 1 2],1),diamondPoints([4 1 2],2),'Color',colorEdges{5},'LineWidth',2); %left
        plot(diamondPoints([4 3 2],1),diamondPoints([4 3 2],2),'Color',colorEdges{4},'LineWidth',2); %right
    else
        plot(xyPoint(1),xyPoint(2),'d','MarkerSize',18,'MarkerEdgeColor','w','MarkerFaceColor','w');
        plot(diamondPoints([4 1 2],1),diamondPoints([4 1 2],2),'Color',colorEdges{5},'LineWidth',2); %left
        plot(diamondPoints([4 3 2],1),diamondPoints([4 3 2],2),'Color',colorEdges{4},'LineWidth',2); %right
    end
end
hold off; xlim([0 max_Weight]); ylim([0 max_Weight]);
title(['Empirical vs. MLE-predicted Aud Weight (N = ' num2str(nSubjects) ')']);
xlabel('A weight (mle)')
ylabel('A weight (empirical)')
legend([h2 h1],{'Significant: EMP < MLE','No (or Aud-) Overweighting'},'Location','northwest')

%[EOF]
