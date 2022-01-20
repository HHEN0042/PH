%initCobraToolbox (0)
%changeCobraSolver('gurobi') 

model = readCbModel('sampleMetaOut_EndoRecon1NewBounds.mat');
model_points = model.points; 

Model_141 = readCbModel('Model_TME141');
Model_142 = readCbModel('Model_TME142');
Model_143 = readCbModel('Model_TME143');
Model_144 = readCbModel('Model_TME144');
Model_145 = readCbModel('Model_TME145');
Model_146 = readCbModel('Model_TME146');
Model_147 = readCbModel('Model_TME147');
Model_148 = readCbModel('Model_TME148');
Model_149 = readCbModel('Model_TME149');
Model_151 = readCbModel('Model_TME151');
Model_152 = readCbModel('Model_TME152');
Model_153 = readCbModel('Model_TME153');
Model_154 = readCbModel('Model_TME154');
Model_155 = readCbModel('Model_TME155');
Model_156 = readCbModel('Model_TME156');
Model_157 = readCbModel('Model_TME157');
Model_158 = readCbModel('Model_TME158');
Model_159 = readCbModel('Model_TME159');
Model_160 = readCbModel('Model_TME160');
Model_161 = readCbModel('Model_TME161');
Model_162 = readCbModel('Model_TME162');
Model_163 = readCbModel('Model_TME163');
Model_164 = readCbModel('Model_TME164');
Model_165 = readCbModel('Model_TME165');
Model_166 = readCbModel('Model_TME166');
Model_167 = readCbModel('Model_TME167');
Model_168 = readCbModel('Model_TME168');
Model_169 = readCbModel('Model_TME169');
Model_170 = readCbModel('Model_TME170');
Model_171 = readCbModel('Model_TME171');
Model_172 = readCbModel('Model_TME172');
Model_173 = readCbModel('Model_TME173');
Model_174 = readCbModel('Model_TME174');
Model_175 = readCbModel('Model_TME175');
Model_176 = readCbModel('Model_TME176');

TasksTable=readtable('Supp_Tables_2_ENERGY METABOLISM2.xls','sheet','Supp_Table_1_Recon1_Filtered');


Tasks=unique(TasksTable(:,1));
% Table, Headers, and colums:
headers=split(unique(join(TasksTable{:,2:4},'_'),'stable'),'_');
data = cell(height(Tasks)+2,length(model.rxns)+3);
for i=1:length(headers)
    data(i+2,1)=headers(i,1);
    data(i+2,2)=headers(i,2);
    data(i+2,3)=headers(i,3);
end
for i=1:length(model.rxns)
    data(1,i+3)=model.subSystems(i);
    data(2,i+3)=model.rxns(i);
end
FunctionSummaryTrauma141=transpose(data);
FunctionSummaryTrauma142=transpose(data);
FunctionSummaryTrauma143=transpose(data);
FunctionSummaryTrauma144=transpose(data);
FunctionSummaryTrauma145=transpose(data);
FunctionSummaryTrauma146=transpose(data);
FunctionSummaryTrauma147=transpose(data);
FunctionSummaryTrauma148=transpose(data);
FunctionSummaryTrauma149=transpose(data);
FunctionSummaryTrauma151=transpose(data);
FunctionSummaryTrauma152=transpose(data);
FunctionSummaryTrauma153=transpose(data);
FunctionSummaryTrauma154=transpose(data);
FunctionSummaryTrauma155=transpose(data);
FunctionSummaryTrauma156=transpose(data);
FunctionSummaryTrauma157=transpose(data);
FunctionSummaryTrauma158=transpose(data);
FunctionSummaryTrauma159=transpose(data);
FunctionSummaryTrauma160=transpose(data);
FunctionSummaryTrauma161=transpose(data);
FunctionSummaryTrauma162=transpose(data);
FunctionSummaryTrauma163=transpose(data);
FunctionSummaryTrauma164=transpose(data);
FunctionSummaryTrauma165=transpose(data);
FunctionSummaryTrauma166=transpose(data);
FunctionSummaryTrauma167=transpose(data);
FunctionSummaryTrauma168=transpose(data);
FunctionSummaryTrauma169=transpose(data);
FunctionSummaryTrauma170=transpose(data);
FunctionSummaryTrauma171=transpose(data);
FunctionSummaryTrauma172=transpose(data);
FunctionSummaryTrauma173=transpose(data);
FunctionSummaryTrauma174=transpose(data);
FunctionSummaryTrauma175=transpose(data);
FunctionSummaryTrauma176=transpose(data);


for i=1:height(Tasks)
    %Ith Function
    IthIndex=find(ismember(TasksTable(:,1),Tasks(i,1)));
    IthTask=join(string(TasksTable{IthIndex(1),2:4}),newline);
%   IthTaskShort=TasksTable{IthIndex(1),4};
    IthTaskShort=join(string(TasksTable{IthIndex(1),2:4}),'_');    
    IthSubstrate=table2cell(TasksTable(IthIndex,5));
    IthProduct=table2cell(TasksTable(IthIndex,8));
    %Testing trauma groups
    [Flux_141, FBAsolution_141,~]=testPathway(Model_141,IthSubstrate,IthProduct);
    [Flux_142, FBAsolution_142,~]=testPathway(Model_142,IthSubstrate,IthProduct);
    [Flux_143, FBAsolution_143,~]=testPathway(Model_143,IthSubstrate,IthProduct);
    [Flux_144, FBAsolution_144,~]=testPathway(Model_144,IthSubstrate,IthProduct);
    [Flux_145, FBAsolution_145,~]=testPathway(Model_145,IthSubstrate,IthProduct);
    [Flux_146, FBAsolution_146,~]=testPathway(Model_146,IthSubstrate,IthProduct);
    [Flux_147, FBAsolution_147,~]=testPathway(Model_147,IthSubstrate,IthProduct);
    [Flux_148, FBAsolution_148,~]=testPathway(Model_148,IthSubstrate,IthProduct);
    [Flux_149, FBAsolution_149,~]=testPathway(Model_149,IthSubstrate,IthProduct);
    [Flux_151, FBAsolution_151,~]=testPathway(Model_151,IthSubstrate,IthProduct);
    [Flux_152, FBAsolution_152,~]=testPathway(Model_152,IthSubstrate,IthProduct);
    [Flux_153, FBAsolution_153,~]=testPathway(Model_153,IthSubstrate,IthProduct);
    [Flux_154, FBAsolution_154,~]=testPathway(Model_154,IthSubstrate,IthProduct);
    [Flux_155, FBAsolution_155,~]=testPathway(Model_155,IthSubstrate,IthProduct);
    [Flux_156, FBAsolution_156,~]=testPathway(Model_156,IthSubstrate,IthProduct);
    [Flux_157, FBAsolution_157,~]=testPathway(Model_157,IthSubstrate,IthProduct);
    [Flux_158, FBAsolution_158,~]=testPathway(Model_158,IthSubstrate,IthProduct);
    [Flux_159, FBAsolution_159,~]=testPathway(Model_159,IthSubstrate,IthProduct);
    [Flux_160, FBAsolution_160,~]=testPathway(Model_160,IthSubstrate,IthProduct);
    [Flux_161, FBAsolution_161,~]=testPathway(Model_161,IthSubstrate,IthProduct);
    [Flux_162, FBAsolution_162,~]=testPathway(Model_162,IthSubstrate,IthProduct);
    [Flux_163, FBAsolution_163,~]=testPathway(Model_163,IthSubstrate,IthProduct);
    [Flux_164, FBAsolution_164,~]=testPathway(Model_164,IthSubstrate,IthProduct);
    [Flux_165, FBAsolution_165,~]=testPathway(Model_165,IthSubstrate,IthProduct);
    [Flux_166, FBAsolution_166,~]=testPathway(Model_166,IthSubstrate,IthProduct);
    [Flux_167, FBAsolution_167,~]=testPathway(Model_167,IthSubstrate,IthProduct);
    [Flux_168, FBAsolution_168,~]=testPathway(Model_168,IthSubstrate,IthProduct);
    [Flux_169, FBAsolution_169,~]=testPathway(Model_169,IthSubstrate,IthProduct);
    [Flux_170, FBAsolution_170,~]=testPathway(Model_170,IthSubstrate,IthProduct);
    [Flux_171, FBAsolution_171,~]=testPathway(Model_171,IthSubstrate,IthProduct);
    [Flux_172, FBAsolution_172,~]=testPathway(Model_172,IthSubstrate,IthProduct);
    [Flux_173, FBAsolution_173,~]=testPathway(Model_173,IthSubstrate,IthProduct);
    [Flux_174, FBAsolution_174,~]=testPathway(Model_174,IthSubstrate,IthProduct);
    [Flux_175, FBAsolution_175,~]=testPathway(Model_175,IthSubstrate,IthProduct);
    [Flux_176, FBAsolution_176,~]=testPathway(Model_176,IthSubstrate,IthProduct);
    
    
    
        for j=1:length(Model1_I.rxns)
        FunctionSummaryTrauma141{j+3,i+2}=FBAsolution_141.x(j);
        FunctionSummaryTrauma142{j+3,i+2}=FBAsolution_142.x(j);
        FunctionSummaryTrauma143{j+3,i+2}=FBAsolution_143.x(j);
        FunctionSummaryTrauma144{j+3,i+2}=FBAsolution_144.x(j);
        FunctionSummaryTrauma145{j+3,i+2}=FBAsolution_145.x(j);
        FunctionSummaryTrauma146{j+3,i+2}=FBAsolution_146.x(j);
        FunctionSummaryTrauma147{j+3,i+2}=FBAsolution_147.x(j);
        FunctionSummaryTrauma148{j+3,i+2}=FBAsolution_148.x(j);
        FunctionSummaryTrauma149{j+3,i+2}=FBAsolution_149.x(j);
        FunctionSummaryTrauma151{j+3,i+2}=FBAsolution_151.x(j);
        FunctionSummaryTrauma152{j+3,i+2}=FBAsolution_152.x(j);
        FunctionSummaryTrauma153{j+3,i+2}=FBAsolution_153.x(j);
        FunctionSummaryTrauma154{j+3,i+2}=FBAsolution_154.x(j);
        FunctionSummaryTrauma155{j+3,i+2}=FBAsolution_155.x(j);
        FunctionSummaryTrauma156{j+3,i+2}=FBAsolution_156.x(j);
        FunctionSummaryTrauma157{j+3,i+2}=FBAsolution_157.x(j);
        FunctionSummaryTrauma158{j+3,i+2}=FBAsolution_158.x(j);
        FunctionSummaryTrauma159{j+3,i+2}=FBAsolution_159.x(j);
        FunctionSummaryTrauma160{j+3,i+2}=FBAsolution_160.x(j);
        FunctionSummaryTrauma161{j+3,i+2}=FBAsolution_161.x(j);
        FunctionSummaryTrauma162{j+3,i+2}=FBAsolution_162.x(j);
        FunctionSummaryTrauma163{j+3,i+2}=FBAsolution_163.x(j);
        FunctionSummaryTrauma164{j+3,i+2}=FBAsolution_164.x(j);
        FunctionSummaryTrauma165{j+3,i+2}=FBAsolution_165.x(j);
        FunctionSummaryTrauma166{j+3,i+2}=FBAsolution_166.x(j);
        FunctionSummaryTrauma167{j+3,i+2}=FBAsolution_167.x(j);
        FunctionSummaryTrauma168{j+3,i+2}=FBAsolution_168.x(j);
        FunctionSummaryTrauma169{j+3,i+2}=FBAsolution_169.x(j);
        FunctionSummaryTrauma170{j+3,i+2}=FBAsolution_170.x(j);
        FunctionSummaryTrauma171{j+3,i+2}=FBAsolution_171.x(j);
        FunctionSummaryTrauma172{j+3,i+2}=FBAsolution_172.x(j);
        FunctionSummaryTrauma173{j+3,i+2}=FBAsolution_173.x(j);
        FunctionSummaryTrauma174{j+3,i+2}=FBAsolution_174.x(j);
        FunctionSummaryTrauma175{j+3,i+2}=FBAsolution_175.x(j);
        FunctionSummaryTrauma176{j+3,i+2}=FBAsolution_176.x(j);
              
        
        end
    subplot(ceil(sqrt(height(Tasks))),floor(sqrt(height(Tasks))),i)
    hold on;    
    t=bar([Flux_141,Flux_142,Flux_143,Flux_144,Flux_145,Flux_146,Flux_147,Flux_148,Flux_149
        Flux_151,Flux_152,Flux_153,Flux_154,Flux_155,Flux_156,Flux_157,Flux_158,Flux_159
        Flux_160,Flux_160,Flux_160,Flux_160,Flux_160,Flux_160,Flux_160,Flux_160,Flux_160 
        Flux_170,Flux_171,Flux_172,Flux_173,Flux_174,Flux_175,Flux_176
        ]);
    ylabel('Flux Size');
    xlabel('Groups');
    title(IthTask);
    ax=gca;
    %save(subplot); 
end