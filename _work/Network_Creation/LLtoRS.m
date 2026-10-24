function LLtoRS(file,type)

load rsMap.mat
load dendriteMap.mat
load activatedCells.mat
fileName = file;

% outputFile = fopen(fileName,'w');
    % Assign a file Id. to fileName
% fprintf(outputFile, 'digraph "LL_to_RS" {\n\n');
%     % Initialize dot file
% fclose(outputFile);
%     % Close outputFile for writing.
    
    numConsPerCell = [16 16 19 19 10 10 7 7 8 8 8 10 10 6 6 12 12 11 11 11]';
    

if type == 1
    
    side = 'L';
    % 1st Section
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeM_',side],activatedCells{1,2,2},6);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeM_',side],activatedCells{1,2,2},6);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLr_',side],activatedCells{1,3,2},1);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLr_',side],activatedCells{1,3,2},6);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLc_',side],activatedCells{1,4,2},2);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLc_',side],activatedCells{1,4,2},2);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLm_',side],activatedCells{1,5,2},1);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLm_',side],activatedCells{1,5,2},1);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_MeM_',side],activatedCells{1,2,2},3);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_MeM_',side],activatedCells{1,2,2},3);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_MeLc_',side],activatedCells{1,4,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_MeLc_',side],activatedCells{1,4,2},1);
    
    % 2nd Section
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM1r_',side],activatedCells{2,2,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM1r_',side],activatedCells{2,2,2},1);
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM1c_',side],activatedCells{2,3,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM1c_',side],activatedCells{2,3,2},1);
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoL1_',side],activatedCells{2,4,2},12);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoL1_',side],activatedCells{2,4,2},12);
    
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    
    % 3rd Section
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    
    % 4th Section
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    
    % 5th Section
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    % 6th Section
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiR2_',side],activatedCells{6,2,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiR2_',side],activatedCells{6,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiV2_',side],activatedCells{6,3,2},5);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiV2_',side],activatedCells{6,3,2},5);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD2cm_',side],activatedCells{6,4,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD2cm_',side],activatedCells{6,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD2cl_',side],activatedCells{6,5,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD2cl_',side],activatedCells{6,5,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    
    % 7th Section
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3cm_',side],activatedCells{7,3,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3cm_',side],activatedCells{7,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    
    % 8th Section
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'18'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'19'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'20'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    % 9th Section 
    
    connectOneWay(fileName,['Sens_',side,'18'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'19'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'20'],1,['Des_T_',side],activatedCells{9,2,2},10);
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_T_',side],activatedCells{9,2,2},10);
    
    
    side = 'R';
    
     % 1st Section
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeM_',side],activatedCells{1,2,2},6);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeM_',side],activatedCells{1,2,2},6);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLr_',side],activatedCells{1,3,2},1);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLr_',side],activatedCells{1,3,2},6);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLc_',side],activatedCells{1,4,2},2);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLc_',side],activatedCells{1,4,2},2);
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_MeLm_',side],activatedCells{1,5,2},1);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_MeLm_',side],activatedCells{1,5,2},1);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_MeM_',side],activatedCells{1,2,2},3);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_MeM_',side],activatedCells{1,2,2},3);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_MeLc_',side],activatedCells{1,4,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_MeLc_',side],activatedCells{1,4,2},1);
    
    % 2nd Section
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM1r_',side],activatedCells{2,2,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM1r_',side],activatedCells{2,2,2},1);
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM1c_',side],activatedCells{2,3,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM1c_',side],activatedCells{2,3,2},1);
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoL1_',side],activatedCells{2,4,2},12);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoL1_',side],activatedCells{2,4,2},12);
    
    connectOneWay(fileName,['Sens_',side,'1'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    connectOneWay(fileName,['Sens_',side,'2'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoL1_',side],activatedCells{2,4,2},6);
    
    % 3rd Section
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'3'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    connectOneWay(fileName,['Sens_',side,'4'],1,['Des_RoM2l_',side],activatedCells{3,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoM2m_',side],activatedCells{3,3,2},1);
    
    % 4th Section
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'5'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    connectOneWay(fileName,['Sens_',side,'6'],1,['Des_RoM3m_',side],activatedCells{4,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_RoV3_',side],activatedCells{4,2,2},3);
    
    % 5th Section
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiR1_',side],activatedCells{5,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiV1_',side],activatedCells{5,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'7'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'8'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_Mauthner_',side],activatedCells{5,4,2},1);
    
    % 6th Section
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiR2_',side],activatedCells{6,2,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiR2_',side],activatedCells{6,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiV2_',side],activatedCells{6,3,2},5);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiV2_',side],activatedCells{6,3,2},5);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD2cm_',side],activatedCells{6,4,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD2cm_',side],activatedCells{6,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD2cl_',side],activatedCells{6,5,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD2cl_',side],activatedCells{6,5,2},1);
    
    connectOneWay(fileName,['Sens_',side,'9'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'10'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'11'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiV2_',side],activatedCells{6,3,2},2);
    
    % 7th Section
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3cm_',side],activatedCells{7,3,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3cm_',side],activatedCells{7,3,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    
    connectOneWay(fileName,['Sens_',side,'12'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    connectOneWay(fileName,['Sens_',side,'13'],1,['Des_MiD3i_',side],activatedCells{7,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_MiD3cl_',side],activatedCells{7,4,2},1);
    
    % 8th Section
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'14'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'15'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    connectOneWay(fileName,['Sens_',side,'18'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'19'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    connectOneWay(fileName,['Sens_',side,'20'],1,['Des_CaD_',side],activatedCells{8,2,2},1);
    
    % 9th Section 
    
    connectOneWay(fileName,['Sens_',side,'18'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'19'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'20'],1,['Des_T_',side],activatedCells{9,2,2},10);
    
    connectOneWay(fileName,['Sens_',side,'16'],1,['Des_T_',side],activatedCells{9,2,2},10);
    connectOneWay(fileName,['Sens_',side,'17'],1,['Des_T_',side],activatedCells{9,2,2},10);
    
    
    
    
    
    
   
    
    
else
%     
%     side = 'L';
%      somaPerSection = [10 14 2 5 3 6 3 1 10]';
%      consPerSection = round(217*somaPerSection./sum(somaPerSection));
%      if sum(consPerSection) < 217
%          consPerSection(10) = consPerSection + 217 - sum(consPerSection);
%      end
%    
%     
%     for i = 1:9
%         for j = 1:consPerSection(i)
%             
%             
%             
%             
%         end
%     end
 
 
 
 
%     
%     for j = 1:2
%         if j == 1
%             side = 'L';
%         else
%             side = 'R';
%         end
%         
%         
%         for i = 1:20
%        
%             k = 0;
%             connections = 0;
%             prevCon = zeros(20,3);
%             
%             while connections <= numConsPerCell(i)
%                 
%             
%                 tmp = activatedCells{mod(k,9)+1,1,1};
%                 
%                 randInt = ceil(tmp*rand(1));
%                 
%                 tempNeuron = prevCon((find(prevCon(:,1) == (mod(k,9)+1))),:);
%                 tempNeuron = tempNeuron((find(tempNeuron(:,2) == (randInt+1))),:);
%                 
%                 badNeurons = tempNeuron(:,3);
%                 
%                 if length(badNeurons) <  activatedCells{mod(k,9)+1,randInt+1,2}             
%                    
%                  nodeNumber = connectOneWay(fileName,['Sens_',side,num2str(i)],1,['Des_',activatedCells{mod(k,9)+1,randInt+1,1},'_',side],activatedCells{mod(k,9)+1,randInt+1,2},1,badNeurons);
%                  prevCon(k+1,1) = mod(k,9)+1; 
%                  prevCon(k+1,2) = randInt+1; 
%                  prevCon(k+1,3) = nodeNumber; 
%                  connections = connections + 1;
%                  
%                 end
%                 
%                 k = k + 1;
%                 
%             end
%     
%         end
%         
%     end
%     
% end



load activeCellsList.mat

outputFile = fopen(file,'a');

fprintf(outputFile,'############################################################################################################################################\n');
fprintf(outputFile,'############################################################################################################################################\n');
fprintf(outputFile,'####################################################     Distributed LL to RS Connect    ###################################################\n');
fprintf(outputFile,'############################################################################################################################################\n');
fprintf(outputFile,'############################################################################################################################################\n\n\n');

for s = 1:2
    if s == 1
        side = 'L';
        tempArray = cell(56);
        tempArray = activeCellsLeft;
        
    else
        side = 'R';
        tempArray = cell(56);
        tempArray = activeCellsRight;
    end
    
    
    
    
    for i=1:20
    
      tempArray = tempArray(randperm(56));
    
    
    
    
    for j = 1:numConsPerCell(i)
    
        fprintf(outputFile,'%s  ',['"Sens_',side,num2str(i),'" -> "',tempArray{j},'";']);
    
        if mod(j,3) == 0
            fprintf(outputFile,'\n');
        end
        
    end
    
    end
end

end

% fprintf(outputFile,'\n\n\n############################################################################################################################################\n')
% fprintf(outputFile,'############################################################################################################################################\n')
% fprintf(outputFile,'##################################################     End of Distributed LL to RS Connect     #############################################\n')
% fprintf(outputFile,'############################################################################################################################################\n')
% fprintf(outputFile,'############################################################################################################################################\n\n\n\n')



















% Create end of Dot file
% outputFile = fopen(fileName,'a');
% fprintf(outputFile, '\n\n } \n\n');    
% fclose(outputFile);  
end



function outputNum = connectOneWay(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon,NotAllowed)
%
% connectOneWay(file,sec1,num1,sec2,num2,numCon)
% 
% Inputs:
%   file - string containing filename to write to
%   sec1 - string of the first section to connect
%   num1 - number of nodes for the first section
%   sec2 - string of the second section
%   num2 - number of nodes for the second section
%   numCon - number of connections between the two sections

if nargin < 7; NotAllowed = 0; end

overflow = 0;
outputNum = 0;


outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//****************** One Way Connect ',sec1,' to ',sec2,' ****************** ']);


sizeOut = nodeNumber1+nodeNumber2; 

tmp = ~eye(sizeOut); % makes a matrix of size sizeOut (both sections concatenated)
tmp(1:nodeNumber1,1:nodeNumber1) = 0; % disallows connections within a single section
tmp(1:sizeOut,(nodeNumber1+1):sizeOut) = 0;
tmp(NotAllowed+nodeNumber1,1:nodeNumber1) = 0;
i = find(tmp); % find all allowed connections
rp = randperm(length(i)); % shuffle connection order to randomize it
irp = i(rp);
if sum(sum(tmp)) < numCon
    overflow = numCon - sum(sum(tmp));
    numCon = sum(sum(tmp));
end
connec = irp(1:numCon); % only take the number of connections desired (numCon)

[I J] = ind2sub([sizeOut, sizeOut], connec); %Convert linear indices to subscripts


% MAIN LOOP
for i = 1:numCon
    fprintf(outputFile,'%s  ',['"',sec1,'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";']);
    if mod(i,3) == 0
        fprintf(outputFile,'\n');
    end
end % End of MAIN LOOP

fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',sec1,' to ',sec2,' **********************']);
fclose(outputFile);

outputNum = I(i)-nodeNumber1;

end