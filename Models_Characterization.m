initCobraToolbox (0)
changeCobraSolver('gurobi') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Characterize the Baseline model   %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: %Random sampling of hpmecgen model a level 3 sized model with level 4a constraints !!!!!!I would say that this and the above three lines are all the 1st step random sampling
%model = readCbModel('iEC2997_preclinicalmodel_annsurg');%%%% The Ccat-model should be used for the clinical models 
load('EC-GEM3006')
model = ExpandedModelNewBoundaries
FBA = optimizeCbModel(model);%%%%This FBA is part of the random sampling protocol it is used to establish the set the boundary of the possible flux predictions to consider
model.lb(find(model.c)) = 0.5*FBA.f;%%%%Here we find half the values of the maximal predicted fluxes
[sampleMetaOutC, mixedFraction] = gpSampler(model, length(model.rxns), [], 8*3600,10000)%,4);%%%%Here we take the half maximal model and set up the random sampling
save mixedFraction_EndoRecon1NewBounds mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_EndoRecon1NewBounds FBA
pointsC = sampleMetaOutC.points;
save sampleMetaOut_EndoRecon1NewBounds.mat sampleMetaOutC
% STEP 2: List of each metabolits after random sampling 
%%% Loaded the model after the random sampling 
model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 
%%%Set up a place to put the next bit
FluxesPreclinical = zeros(length(model.rxns),7);
%%%% Get the stats about the reaction fluxes the 25th and 75th percentiles
%%%% (columns 6 and 7) would be enough but the rest may be interesting
for idx=1:length(model.rxns)    
    flux = model_points(idx,:);        
    minf = min(flux);
    meanf = mean(flux);
    maxf = max(flux);
    stdf = std(flux);
    UPstd = meanf + (stdf*2);
    LOWstd = meanf - (stdf*2);
    percent25 = prctile(flux,25);
    percent75 = prctile(flux,75);    
    FluxesPreclinical(idx,1) = minf;
    FluxesPreclinical(idx,2) = meanf;
    FluxesPreclinical(idx,3) = maxf;
    FluxesPreclinical(idx,4) = LOWstd;
    FluxesPreclinical(idx,5) = UPstd;
    FluxesPreclinical(idx,6) = percent25;
    FluxesPreclinical(idx,7) = percent75;        
end
clear('ans' , 'idx' , 'minf' , 'meanf' , 'mean(flux)' , 'maxf' , 'stdf' , 'UPstd' , 'LOWstd' , 'percent25' , 'percent75');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   Define Boundaries in Trauma models    %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MODEL 141
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'A1:B51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME141=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME141); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME141_pre model_TME141

%% MODEL 142
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'C1:D51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME142=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME142); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME142_pre model_TME142


%% MODEL 143
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'E1:F51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME143=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME143); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME143_pre model_TME143


%% MODEL 144
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'G1:H51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME144=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME144); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME144_pre model_TME144

%% MODEL 145
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'I1:J51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME145=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME145); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME145_pre model_TME145


%% MODEL 146
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'K1:L51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME146=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME146); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME146_pre model_TME146


%% MODEL 147
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'M1:N51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME147=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME147); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME147_pre model_TME147


%% MODEL 148
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'O1:P51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME148=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME148); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME148_pre model_TME148


%% MODEL 149
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'Q1:R51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME149=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME149); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME149_pre model_TME149


%% MODEL 151
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'S1:T51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME151=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME151); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME151_pre model_TME151


%% MODEL 152
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'U1:V51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME152=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME152); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME152_pre model_TME152

%% MODEL 153
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'W1:X51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME153=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME153); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME153_pre model_TME153


%% MODEL 154
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'Y1:Z51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME154=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME154); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME154_pre model_TME154


%% MODEL 155
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AA1:AB51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME155=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME155); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME155_pre model_TME155


%% MODEL 156
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AC1:AD51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME156=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME156); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME156_pre model_TME156


%% MODEL 157
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AE1:AF51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME157=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME157); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME157_pre model_TME157


%% MODEL 158
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AG1:AH51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME158=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME158); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME158_pre model_TME158


%% MODEL 159
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AI1:AJ51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME159=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME159); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME159_pre model_TME159


%% MODEL 160
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AK1:AL51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME160=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME160); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME160_pre model_TME160


%% MODEL 161
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AM1:AN51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME161=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME161); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME161_pre model_TME161


%% MODEL 162
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AO1:AP51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME162=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME162); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME162_pre model_TME162


%% MODEL 163
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AQ1:AR51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME163=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME163); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME163_pre model_TME163


%% MODEL 164
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AS1:AT51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME164=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME164); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME164_pre model_TME164


%% MODEL 165
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AU1:AV51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME165=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME165); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME165_pre model_TME165


%% MODEL 166
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AW1:AX51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME166=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME166); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME166_pre model_TME166


%% MODEL 167
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'AY1:AZ51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME167=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME167); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME167_pre model_TME167


%% MODEL 168
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BA1:BB51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME168=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME168); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME168_pre model_TME168



%% MODEL 169
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BC1:BD51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME169=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME169); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME169_pre model_TME169

%% MODEL 170
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BE1:BF51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME170=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME170); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME170_pre model_TME170


%% MODEL 171
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BG1:BH51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME171=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME171); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME171_pre model_TME171


%% MODEL 172
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BI1:BJ51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME172=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME172); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME172_pre model_TME172

%% MODEL 173
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BK1:BL51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME173=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME173); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME173_pre model_TME173

%% MODEL 174
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BM1:BN51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME174=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME174); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME174_pre model_TME174


%% MODEL 175
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BO1:BP51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME175=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME175); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME175_pre model_TME175


%% MODEL 176
% STEP 3: Adjusting rxns bounds to 25/75 intervals - this part will be done for each model (A, B, C, D) i.e. I will end up with four models
% Read in the FCs
% MODEL 1
[NUMUrxn1,TXTUrxn1,RAWUrxn1] = xlsread('FC til model _FINAL.xlsx', 1, 'BQ1:BR51'); %% FC for the mean_patient_group_X / Mean_controls
 for i = 1:length(TXTUrxn1) 
     idx = find(ismember(model.rxns,TXTUrxn1(i)));
     meanf = FluxesPreclinical(idx,2);
     dir1 = FluxesPreclinical(idx,6); %Direction of the flux in the 25 percentile
     dir2 = FluxesPreclinical(idx,7); %Direction of the flux in the 75 percentile
     if abs(dir2)+abs(dir1)==0 % if Control UB = LB = 0
         dir1 = -10;
         dir2 = 10;
%         newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution  
%         newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));  %
%         Only if the Ex fluxes have normal distribution      
     end
%     else %Only if the Ex fluxes have normal distribution
     newb(i,1) = dir1;
     newb(i,2) = dir2;
     newb(i,3) = dir1 + (dir1/abs(dir1))*(dir1*(NUMUrxn1(i)-1));   
     newb(i,4) = dir2 + (dir2/abs(dir2))*(dir2*(NUMUrxn1(i)-1));
     %        newb(i,3) = dir1 + (dir1/abs(dir1))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution
     %        newb(i,4) = dir2 + (dir2/abs(dir2))*(meanf*(NUMUrxn1(i)-1)); %
     %        Only if the Ex fluxes have normal distribution 
     if newb(i,1)<0 && newb(i,3)>0
         newb(i,3) = dir1/NUMUrxn1(i);
     end
     if newb(i,2)>0 && newb(i,4)<0
         newb(i,4) = dir2*NUMUrxn1(i);
     end   
 end
% end %Only if the Ex fluxes have normal distribution
 
% STEP 4: Adjust the list from step 3   
modelTrauma = model;
Int=(cell2mat(arrayfun(@(s,e) (s:e), 1, 20, 'UniformOutput', false))/10);
for i = 1:length(TXTUrxn1)
    modelTraumaIni=modelTrauma;
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,3),'l');
    modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),newb(i,4),'u');
    TestTmpModel = optimizeCbModel(modelTrauma);
    if ~isnan(TestTmpModel.f)==1
        modelTraumaIni=modelTrauma;
    else
        modelTrauma=modelTraumaIni;
        flag=0;
        j=1;
        while flag==0
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,3)-(newb(i,3)/abs(newb(i,3)))*(newb(i,3)*Int(j))),'l');
            modelTrauma = changeRxnBounds(modelTrauma,TXTUrxn1(i),(newb(i,4)+(newb(i,4)/abs(newb(i,4)))*(newb(i,4)*Int(j))),'u');                   
            TestTmpModel = optimizeCbModel(modelTrauma);
%           for testing purposes, otherwise comment 
%           i
%           j
%           TestTmpModel.f              
            if ~isnan(TestTmpModel.f)==1
                modelTraumaIni=modelTrauma;
                flag=1;
            elseif j==length(Int)
                modelTrauma=modelTraumaIni;
                flag=1;  
            else
                modelTrauma=modelTraumaIni;
                j=j+1;
            end
%           for testing purposes, otherwise comment 
%           flag    
        end 
    end
end
% STEP 5: Use FBA to see if the model with the trauma constraints works
model_TME176=modelTrauma;
testFBA_model1_I = optimizeCbModel(model_TME176); %%%if this FAILS relax reactions to produce a feasible model if this works proceed to step 7 to random sample the newly defined clinical model
save model_TME176_pre model_TME176

%%RANDOM SAMPLING

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME141);
model.lb(find(model_TME141.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME141, length(model_TME141.rxns), [], 8*3600,10000);
save mixedFraction_model_TME141 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME141 fba
pointsC = sampleMetaOutC.points;
save points_model_TME141 pointsC
save Model_TME141.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME142);
model.lb(find(model_TME142.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME142, length(model_TME142.rxns), [], 8*3600,10000);
save mixedFraction_model_TME142 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME142 fba
pointsC = sampleMetaOutC.points;
save points_model_TME142 pointsC
save Model_TME142.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME143);
model.lb(find(model_TME143.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME143, length(model_TME143.rxns), [], 8*3600,10000);
save mixedFraction_model_TME143 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME143 fba
pointsC = sampleMetaOutC.points;
save points_model_TME143 pointsC
save Model_TME143.mat sampleMetaOutC




% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME144);
model.lb(find(model_TME144.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME144, length(model_TME144.rxns), [], 8*3600,10000);
save mixedFraction_model_TME144 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME144 fba
pointsC = sampleMetaOutC.points;
save points_model_TME144 pointsC
save Model_TME144.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME145);
model.lb(find(model_TME145.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME145, length(model_TME145.rxns), [], 8*3600,10000);
save mixedFraction_model_TME145 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME145 fba
pointsC = sampleMetaOutC.points;
save points_model_TME145 pointsC
save Model_TME145.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME146);
model.lb(find(model_TME146.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME146, length(model_TME146.rxns), [], 8*3600,10000);
save mixedFraction_model_TME146 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME146 fba
pointsC = sampleMetaOutC.points;
save points_model_TME146 pointsC
save Model_TME146.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME147);
model.lb(find(model_TME147.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME147, length(model_TME147.rxns), [], 8*3600,10000);
save mixedFraction_model_TME147 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME147 fba
pointsC = sampleMetaOutC.points;
save points_model_TME147 pointsC
save Model_TME147.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME148);
model.lb(find(model_TME148.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME148, length(model_TME148.rxns), [], 8*3600,10000);
save mixedFraction_model_TME148 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME148 fba
pointsC = sampleMetaOutC.points;
save points_model_TME148 pointsC
save Model_TME148.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME149);
model.lb(find(model_TME149.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME149, length(model_TME149.rxns), [], 8*3600,10000);
save mixedFraction_model_TME149 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME149 fba
pointsC = sampleMetaOutC.points;
save points_model_TME149 pointsC
save Model_TME149.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME151);
model.lb(find(model_TME151.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME151, length(model_TME151.rxns), [], 8*3600,10000);
save mixedFraction_model_TME151 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME151 fba
pointsC = sampleMetaOutC.points;
save points_model_TME151 pointsC
save Model_TME151.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME152);
model.lb(find(model_TME152.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME152, length(model_TME152.rxns), [], 8*3600,10000);
save mixedFraction_model_TME152 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME152 fba
pointsC = sampleMetaOutC.points;
save points_model_TME152 pointsC
save Model_TME152.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME153);
model.lb(find(model_TME153.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME153, length(model_TME153.rxns), [], 8*3600,10000);
save mixedFraction_model_TME153 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME153 fba
pointsC = sampleMetaOutC.points;
save points_model_TME153 pointsC
save Model_TME153.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME154);
model.lb(find(model_TME154.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME154, length(model_TME154.rxns), [], 8*3600,10000);
save mixedFraction_model_TME154 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME154 fba
pointsC = sampleMetaOutC.points;
save points_model_TME154 pointsC
save Model_TME154.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME155);
model.lb(find(model_TME155.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME155, length(model_TME155.rxns), [], 8*3600,10000);
save mixedFraction_model_TME155 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME155 fba
pointsC = sampleMetaOutC.points;
save points_model_TME155 pointsC
save Model_TME155.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME156);
model.lb(find(model_TME156.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME156, length(model_TME156.rxns), [], 8*3600,10000);
save mixedFraction_model_TME156 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME156 fba
pointsC = sampleMetaOutC.points;
save points_model_TME156 pointsC
save Model_TME156.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME157);
model.lb(find(model_TME157.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME157, length(model_TME157.rxns), [], 8*3600,10000);
save mixedFraction_model_TME157 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME157 fba
pointsC = sampleMetaOutC.points;
save points_model_TME157 pointsC
save Model_TME157.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME158);
model.lb(find(model_TME158.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME158, length(model_TME158.rxns), [], 8*3600,10000);
save mixedFraction_model_TME158 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME158 fba
pointsC = sampleMetaOutC.points;
save points_model_TME158 pointsC
save Model_TME158.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME159);
model.lb(find(model_TME159.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME159, length(model_TME159.rxns), [], 8*3600,10000);
save mixedFraction_model_TME159 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME159 fba
pointsC = sampleMetaOutC.points;
save points_model_TME159 pointsC
save Model_TME159.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME160);
model.lb(find(model_TME160.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME160, length(model_TME160.rxns), [], 8*3600,10000);
save mixedFraction_model_TME160 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME160 fba
pointsC = sampleMetaOutC.points;
save points_model_TME160 pointsC
save Model_TME160.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME161);
model.lb(find(model_TME161.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME161, length(model_TME161.rxns), [], 8*3600,10000);
save mixedFraction_model_TME161 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME161 fba
pointsC = sampleMetaOutC.points;
save points_model_TME161 pointsC
save Model_TME161.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME162);
model.lb(find(model_TME162.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME162, length(model_TME162.rxns), [], 8*3600,10000);
save mixedFraction_model_TME162 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME162 fba
pointsC = sampleMetaOutC.points;
save points_model_TME162 pointsC
save Model_TME162.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME163);
model.lb(find(model_TME163.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME163, length(model_TME163.rxns), [], 8*3600,10000);
save mixedFraction_model_TME163 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME163 fba
pointsC = sampleMetaOutC.points;
save points_model_TME163 pointsC
save Model_TME163.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME164);
model.lb(find(model_TME164.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME164, length(model_TME164.rxns), [], 8*3600,10000);
save mixedFraction_model_TME164 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME164 fba
pointsC = sampleMetaOutC.points;
save points_model_TME164 pointsC
save Model_TME164.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME165);
model.lb(find(model_TME165.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME165, length(model_TME165.rxns), [], 8*3600,10000);
save mixedFraction_model_TME165 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME165 fba
pointsC = sampleMetaOutC.points;
save points_model_TME165 pointsC
save Model_TME165.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME166);
model.lb(find(model_TME166.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME166, length(model_TME166.rxns), [], 8*3600,10000);
save mixedFraction_model_TME166 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME166 fba
pointsC = sampleMetaOutC.points;
save points_model_TME166 pointsC
save Model_TME166.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME167);
model.lb(find(model_TME167.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME167, length(model_TME167.rxns), [], 8*3600,10000);
save mixedFraction_model_TME167 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME167 fba
pointsC = sampleMetaOutC.points;
save points_model_TME167 pointsC
save Model_TME167.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME168);
model.lb(find(model_TME168.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME168, length(model_TME168.rxns), [], 8*3600,10000);
save mixedFraction_model_TME168 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME168 fba
pointsC = sampleMetaOutC.points;
save points_model_TME168 pointsC
save Model_TME168.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME169);
model.lb(find(model_TME169.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME169, length(model_TME169.rxns), [], 8*3600,10000);
save mixedFraction_model_TME169 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME169 fba
pointsC = sampleMetaOutC.points;
save points_model_TME169 pointsC
save Model_TME169.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME170);
model.lb(find(model_TME170.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME170, length(model_TME170.rxns), [], 8*3600,10000);
save mixedFraction_model_TME170 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME170 fba
pointsC = sampleMetaOutC.points;
save points_model_TME170 pointsC
save Model_TME170.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME171);
model.lb(find(model_TME171.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME171, length(model_TME171.rxns), [], 8*3600,10000);
save mixedFraction_model_TME171 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME171 fba
pointsC = sampleMetaOutC.points;
save points_model_TME171 pointsC
save Model_TME171.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME172);
model.lb(find(model_TME172.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME172, length(model_TME172.rxns), [], 8*3600,10000);
save mixedFraction_model_TME172 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME172 fba
pointsC = sampleMetaOutC.points;
save points_model_TME172 pointsC
save Model_TME172.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME173);
model.lb(find(model_TME173.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME173, length(model_TME173.rxns), [], 8*3600,10000);
save mixedFraction_model_TME173 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME173 fba
pointsC = sampleMetaOutC.points;
save points_model_TME173 pointsC
save Model_TME173.mat sampleMetaOutC



% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME174);
model.lb(find(model_TME174.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME174, length(model_TME174.rxns), [], 8*3600,10000);
save mixedFraction_model_TME174 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME174 fba
pointsC = sampleMetaOutC.points;
save points_model_TME174 pointsC
save Model_TME174.mat sampleMetaOutC

% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME175);
model.lb(find(model_TME175.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME175, length(model_TME175.rxns), [], 8*3600,10000);
save mixedFraction_model_TME175 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME175 fba
pointsC = sampleMetaOutC.points;
save points_model_TME175 pointsC
save Model_TME175.mat sampleMetaOutC


% STEP 7: Random sampling to characterize Trauma models
fba = optimizeCbModel(model_TME176);
model.lb(find(model_TME176.c)) = 0.5*fba.f;
[sampleMetaOutC, mixedFraction] = gpSampler(model_TME176, length(model_TME176.rxns), [], 8*3600,10000);
save mixedFraction_model_TME176 mixedFraction%%%%Some sort of naming convention is needed but can be anything here is the model and month
save FBA_model_TME176 fba
pointsC = sampleMetaOutC.points;
save points_model_TME176 pointsC
save Model_TME176.mat sampleMetaOutC



