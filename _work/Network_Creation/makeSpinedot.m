function makeSpinedot(file)



%*************************** Initilize Parameters *************************



% Create beginning of Dot file
fileName = file;
outputFile = fopen(fileName,'a');
fprintf(outputFile, '\n\n subgraph cluster_spine {\n\n');    
fclose(outputFile);  


% The number of spinal Segments in the model
number_of_Segments = 25;

% The number of Interneurons per side
number_of_interneurons = 14;

% The number of primary motor neurons per side
number_of_primary = 3;

% The number of secondary motor neurons per side
number_of_secondary = 10;



% The total number of connections allotted to the spine
TotalSpineCon = 50000;


% The number of connections given to segment clustering connections
cons4SpineClust = round(.4*TotalSpineCon);

% The number of connections given to interneuron to interneuron traverse connections
% cons4IntTraverse = round(.4*TotalSpineCon);
% 
% % The number of connections given to interneuron to motor neurons
% %cons4IntMot = round(.2*TotalSpineCon);
% 
% %The number of connections given to interneuron to motor traverse connections
% cons4IntMotTraverse = round(.2*TotalSpineCon);



%**************************************************************************
%**************************************************************************
%**************************** True File contents **************************
%**************************************************************************
%**************************************************************************
 
%##########################################################################

%**************************************************************************
% Start of spinal segment clustering
%**************************************************************************
makeSectionTypeBreakStart(fileName,'Segment Clustering')

for i = 1:number_of_Segments
    defineSegment(fileName,i,number_of_interneurons,number_of_primary,number_of_secondary)
end

numConInt = round(.6*cons4SpineClust/number_of_Segments);

for i = 1:number_of_Segments
    sec1 = [num2str(i),'In_L'];
    sec2 = [num2str(i),'In_R'];    
    connect2sections(fileName,sec1,number_of_interneurons,sec2,number_of_interneurons,numConInt,true)
end

numConSec = round(.15*cons4SpineClust/number_of_Segments);

for i = 1:number_of_Segments
    sec1 = [num2str(i),'In_L'];
    sec2 = [num2str(i),'Sec_L'];
    connectOneWay(fileName,sec1,number_of_interneurons,sec2,number_of_secondary,numConSec)
end
for i = 1:number_of_Segments
    sec1 = [num2str(i),'In_R'];
    sec2 = [num2str(i),'Sec_R'];
    connectOneWay(fileName,sec1,number_of_interneurons,sec2,number_of_secondary,numConSec)
end

numConPri = round(.05*cons4SpineClust/number_of_Segments);

for i = 1:number_of_Segments
    sec1 = [num2str(i),'In_L'];
    sec2 = [num2str(i),'Pri_L'];
    connectOneWay(fileName,sec1,number_of_interneurons,sec2,number_of_primary,numConPri)
end
for i = 1:number_of_Segments
    sec1 = [num2str(i),'In_R'];
    sec2 = [num2str(i),'Pri_R'];
    connectOneWay(fileName,sec1,number_of_interneurons,sec2,number_of_primary,numConPri)
end


outputFile = fopen(fileName,'a');
fprintf(outputFile, '\n}\n');    
fclose(outputFile);


makeSectionTypeBreakEnd(fileName,'Segment Clustering')
%**************************************************************************
% End of spinal segment clustering
%**************************************************************************
 
%##########################################################################

%**************************************************************************
% Start of spinal Interneuron Traverse Connections
%**************************************************************************
makeSectionTypeBreakStart(fileName,'Interneuron Traverse Connections')

numConInt = round(cons4SpineClust/number_of_Segments);

% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','In_R_1',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','In_L_1',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','In_R_2',10)
connectAscending(fileName,number_of_Segments,'In_L_3','In_R_3',10)
connectAscending(fileName,number_of_Segments,'In_R_2','In_L_2',10)
connectAscending(fileName,number_of_Segments,'In_R_3','In_L_3',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','In_R_4',5)
connectAscending(fileName,number_of_Segments,'In_R_4','In_L_4',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','In_R_5',13)
connectDescending(fileName,number_of_Segments,'In_L_6','In_R_6',13)
connectDescending(fileName,number_of_Segments,'In_R_5','In_L_5',13)
connectDescending(fileName,number_of_Segments,'In_R_6','In_L_6',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','In_R_7',13)
connectDescending(fileName,number_of_Segments,'In_L_8','In_R_8',13)
connectDescending(fileName,number_of_Segments,'In_R_7','In_L_7',13)
connectDescending(fileName,number_of_Segments,'In_R_8','In_L_8',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','In_R_9',14)
connectDescending(fileName,number_of_Segments,'In_L_10','In_R_10',14)
connectDescending(fileName,number_of_Segments,'In_R_9','In_L_9',14)
connectDescending(fileName,number_of_Segments,'In_R_10','In_L_10',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','In_L_11',9)
connectDescending(fileName,number_of_Segments,'In_L_12','In_L_12',9)
connectDescending(fileName,number_of_Segments,'In_R_11','In_R_11',9)
connectDescending(fileName,number_of_Segments,'In_R_12','In_R_12',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','In_L_13',4)
connectAscending(fileName,number_of_Segments,'In_R_13','In_L_13',4)
connectDescending(fileName,number_of_Segments,'In_L_13','In_R_13',4)
connectAscending(fileName,number_of_Segments,'In_L_13','In_R_13',4)
connectDescending(fileName,number_of_Segments,'In_R_14','In_L_14',4)
connectAscending(fileName,number_of_Segments,'In_R_14','In_L_14',4)
connectDescending(fileName,number_of_Segments,'In_L_14','In_R_14',4)
connectAscending(fileName,number_of_Segments,'In_L_14','In_R_14',4)


makeSectionTypeBreakEnd(fileName,'Interneuron Traverse Connections')
%**************************************************************************
% End of spinal Interneuron Traverse Connections
%**************************************************************************
 
%##########################################################################

%**************************************************************************
% Start of spinal Interneuron to Motor Traverse Connections
%**************************************************************************
makeSectionTypeBreakStart(fileName,'Interneuron to Motor Traverse Connections')

% Interneuron to Primary motor neurons

% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Pri_R_1',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Pri_L_1',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Pri_R_2',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Pri_R_2',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Pri_L_2',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Pri_L_2',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Pri_R_3',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Pri_L_3',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Pri_R_1',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Pri_R_1',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Pri_L_1',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Pri_L_1',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Pri_R_2',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Pri_R_2',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Pri_L_2',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Pri_L_2',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Pri_R_3',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Pri_R_3',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Pri_L_3',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Pri_L_3',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Pri_L_1',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Pri_L_1',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Pri_R_1',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Pri_R_1',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Pri_L_2',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Pri_L_2',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Pri_R_2',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Pri_R_2',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Pri_L_2',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Pri_L_2',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Pri_R_2',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Pri_R_2',4)

% 2nd set **********************************************************
% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Pri_R_2',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Pri_L_2',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Pri_R_3',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Pri_R_3',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Pri_L_3',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Pri_L_3',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Pri_R_1',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Pri_L_1',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Pri_R_2',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Pri_R_2',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Pri_L_2',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Pri_L_2',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Pri_R_3',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Pri_R_3',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Pri_L_3',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Pri_L_3',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Pri_R_1',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Pri_R_1',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Pri_L_1',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Pri_L_1',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Pri_L_2',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Pri_L_2',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Pri_R_2',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Pri_R_2',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Pri_L_3',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Pri_L_3',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Pri_R_3',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Pri_R_3',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Pri_L_3',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Pri_L_3',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Pri_R_3',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Pri_R_3',4)

%***********************************************************************
% Interneuron to Secondary motor neurons

% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Sec_R_1',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Sec_L_1',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Sec_R_2',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Sec_R_2',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Sec_L_2',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Sec_L_2',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Sec_R_3',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Sec_L_3',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Sec_R_4',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Sec_R_4',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Sec_L_4',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Sec_L_4',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Sec_R_5',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Sec_R_5',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Sec_L_5',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Sec_L_5',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Sec_R_6',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Sec_R_6',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Sec_L_6',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Sec_L_6',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Sec_L_7',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Sec_L_7',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Sec_R_7',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Sec_R_7',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Sec_L_8',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Sec_L_8',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Sec_R_8',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Sec_R_8',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Sec_L_8',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Sec_L_8',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Sec_R_8',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Sec_R_8',4)

% 2nd set ************************************************************
% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Sec_R_2',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Sec_L_2',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Sec_R_3',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Sec_R_3',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Sec_L_3',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Sec_L_3',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Sec_R_4',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Sec_L_4',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Sec_R_5',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Sec_R_5',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Sec_L_5',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Sec_L_5',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Sec_R_6',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Sec_R_6',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Sec_L_6',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Sec_L_6',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Sec_R_7',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Sec_R_7',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Sec_L_7',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Sec_L_7',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Sec_L_8',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Sec_L_8',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Sec_R_8',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Sec_R_8',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Sec_L_9',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Sec_L_9',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Sec_R_9',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Sec_R_9',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Sec_L_9',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Sec_L_9',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Sec_R_9',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Sec_R_9',4)

% 3rd set ************************************************************

% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Sec_R_3',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Sec_L_3',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Sec_R_4',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Sec_R_4',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Sec_L_4',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Sec_L_4',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Sec_R_5',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Sec_L_5',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Sec_R_6',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Sec_R_6',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Sec_L_6',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Sec_L_6',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Sec_R_7',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Sec_R_7',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Sec_L_7',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Sec_L_7',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Sec_R_8',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Sec_R_8',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Sec_L_8',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Sec_L_8',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Sec_L_9',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Sec_L_9',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Sec_R_9',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Sec_R_9',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Sec_L_10',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Sec_L_10',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Sec_R_10',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Sec_R_10',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Sec_L_10',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Sec_L_10',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Sec_R_10',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Sec_R_10',4)

% 4th set *************************************************************

% The CoPA interneuron
connectAscending(fileName,number_of_Segments,'In_L_1','Sec_R_4',number_of_Segments)
connectAscending(fileName,number_of_Segments,'In_R_1','Sec_L_4',number_of_Segments)

% The CoSA interneurons
connectAscending(fileName,number_of_Segments,'In_L_2','Sec_R_5',10)
connectAscending(fileName,number_of_Segments,'In_L_3','Sec_R_5',10)
connectAscending(fileName,number_of_Segments,'In_R_2','Sec_L_5',10)
connectAscending(fileName,number_of_Segments,'In_R_3','Sec_L_5',10)

% The CoLA interneuron
connectAscending(fileName,number_of_Segments,'In_L_4','Sec_R_6',5)
connectAscending(fileName,number_of_Segments,'In_R_4','Sec_L_6',5)

% The CiD interneurons
connectDescending(fileName,number_of_Segments,'In_L_5','Sec_R_7',13)
connectDescending(fileName,number_of_Segments,'In_L_6','Sec_R_7',13)
connectDescending(fileName,number_of_Segments,'In_R_5','Sec_L_7',13)
connectDescending(fileName,number_of_Segments,'In_R_6','Sec_L_7',13)

% The MCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_7','Sec_R_8',13)
connectDescending(fileName,number_of_Segments,'In_L_8','Sec_R_8',13)
connectDescending(fileName,number_of_Segments,'In_R_7','Sec_L_8',13)
connectDescending(fileName,number_of_Segments,'In_R_8','Sec_L_8',13)

% The UCoD interneurons
connectDescending(fileName,number_of_Segments,'In_L_9','Sec_R_9',14)
connectDescending(fileName,number_of_Segments,'In_L_10','Sec_R_9',14)
connectDescending(fileName,number_of_Segments,'In_R_9','Sec_L_9',14)
connectDescending(fileName,number_of_Segments,'In_R_10','Sec_L_9',14)

% The VeMe interneurons
connectDescending(fileName,number_of_Segments,'In_L_11','Sec_L_10',9)
connectDescending(fileName,number_of_Segments,'In_L_12','Sec_L_10',9)
connectDescending(fileName,number_of_Segments,'In_R_11','Sec_R_10',9)
connectDescending(fileName,number_of_Segments,'In_R_12','Sec_R_10',9)

% The CoBL interneurons
connectDescending(fileName,number_of_Segments,'In_R_13','Sec_L_1',4)
connectAscending(fileName,number_of_Segments,'In_R_13','Sec_L_1',4)
connectDescending(fileName,number_of_Segments,'In_L_13','Sec_R_1',4)
connectAscending(fileName,number_of_Segments,'In_L_13','Sec_R_1',4)
connectDescending(fileName,number_of_Segments,'In_R_14','Sec_L_1',4)
connectAscending(fileName,number_of_Segments,'In_R_14','Sec_L_1',4)
connectDescending(fileName,number_of_Segments,'In_L_14','Sec_R_1',4)
connectAscending(fileName,number_of_Segments,'In_L_14','Sec_R_1',4)


makeSectionTypeBreakEnd(fileName,'Interneuron to Motor Traverse Connections')
%**************************************************************************
% End of spinal Interneuron to Motor Traverse Connections
%**************************************************************************
 
%##########################################################################

%**************************************************************************
%**************************************************************************
%************************* End True File contents *************************
%**************************************************************************
%**************************************************************************


 %Create end of Dot file
 outputFile = fopen(fileName,'a');
 fprintf(outputFile, '\n\n } \n\n');    
 fclose(outputFile);  


end % END of makeRS function



function makeSectionTypeBreakStart(file,section)

outputFile = fopen(file,'a');

fprintf(outputFile,'\n//*************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'%s\n',['//                   Starting Section ',section,'']);
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n\n');

fclose(outputFile);

end


function makeSectionTypeBreakEnd(file,section)

outputFile = fopen(file,'a');

fprintf(outputFile,'\n//*************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'%s\n',['//                      Ending Section ',section,'']);
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n\n\n\n');

fclose(outputFile);

end



function clusterInt(file,seg,numInt,numCon)
% This is a process adapted from Olaf Sporns makerandCIJ_dir.m from the
% brain connectivity toolbox.

% open file and print header
outputFile = fopen(file,'a');
fprintf(outputFile,'\n\n %s \n\n',['//************************ Clustering for ',seg,' ********************** ']);

% randomly create the connections we are going to use
tmp = ~eye(numInt);
i = find(tmp);
rp = randperm(length(i));
irp = i(rp);
i = irp(1:numCon);

% convert linear indices to subscripts
[I J] = ind2sub([numInt, numInt], i);

% MAIN LOOP
for k=1:length(I)
    
    fprintf(outputFile, '%s', ['"', seg, '_', num2str(J(k)), '" -> "', seg, '_', num2str(I(k)), '";  ']);
    
    if mod(k,3)==0
        fprintf(outputFile, '\n');
    end
        
end % End of MAIN LOOP

% print out comment and close file
fprintf(outputFile,'\n\n %s \n\n',['//************************ END of section ',seg,' **********************']);
fclose(outputFile);

end % END of clusterInt function




function defineSegment(file,segNum,INT,PRI,SEC)

% open file and print header
outputFile = fopen(file,'a');
fprintf(outputFile,'\n\n %s \n\n',['//******************* Defintitions for Segment #',num2str(segNum),' ******************* ']);

counter = 0;

% Define all neurons to be within a spinal segment
fprintf(outputFile,'%s \n',['subgraph cluster_SpinalSegment_',num2str(segNum),'  { ']);
fprintf(outputFile,'%s \n',['label = "SpinalSegment_',num2str(segNum),'";  ']);

% Left Side of spinal cord
fprintf(outputFile,'%s \n',['subgraph cluster_Left_',num2str(segNum),'  { ']);
fprintf(outputFile,'%s \n',['label = "Left_',num2str(segNum),'" ;  ']);
% Interneurons
for i = 1:INT
        fprintf(outputFile,'%s',['"',num2str(segNum),'In_L_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
% Primary Motor Neurons
for i = 1:PRI
        fprintf(outputFile,'%s',['"',num2str(segNum),'Pri_L_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
% Secondary Motor Neurons
for i = 1:SEC
        fprintf(outputFile,'%s',['"',num2str(segNum),'Sec_L_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
fprintf(outputFile,'\n                         } \n\n');

% Right Side of spinal cord
fprintf(outputFile,'%s \n',['subgraph cluster_Right_',num2str(segNum),'  { ']);
fprintf(outputFile,'%s \n',['label = "Right_',num2str(segNum),'"; ']);
% Interneurons
for i = 1:INT
        fprintf(outputFile,'%s',['"',num2str(segNum),'In_R_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
% Primary Motor Neurons
for i = 1:PRI
        fprintf(outputFile,'%s',['"',num2str(segNum),'Pri_R_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
% Secondary Motor Neurons
for i = 1:SEC
        fprintf(outputFile,'%s',['"',num2str(segNum),'Sec_R_',num2str(i),'";  ']);
        counter = counter + 1;
        if mod(counter,6) == 0
            fprintf(outputFile,'\n');
        end
end
fprintf(outputFile,'\n                         } \n');
fprintf(outputFile,'\n          } \n\n');

end




function connect2sections(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon,selfConnects)
%
% connect2sections(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon)
% 
% Inputs:
%   file - string containing filename to write to
%   sec1 - string of the first section to connect
%   nodeNumber1 - number of nodes for the first section
%   sec2 - string of the second section
%   nodeNumber2 - number of nodes for the second section
%   numCon - number of connections between the two sections


outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//***************** Two Way Connect ',sec1,' to ',sec2,' ****************** ']);


sizeOut = nodeNumber1+nodeNumber2; 

tmp = ~eye(sizeOut); % makes a matrix of size sizeOut (both sections concatenated)
if selfConnects == false % disallows connections within a single section if false
    tmp(1:nodeNumber1,1:nodeNumber1) = 0; 
    tmp((nodeNumber1+1):sizeOut,(nodeNumber1+1):sizeOut) = 0;
end
i = find(tmp); % find all allowed connections
rp = randperm(length(i)); % shuffle connection order to randomize it
irp = i(rp);
connec = irp(1:numCon); % only take the number of connections desired (numCon)

[I J] = ind2sub([sizeOut, sizeOut], connec); %Convert linear indices to subscripts


% MAIN LOOP
for i = 1:numCon
    if I(i) <= nodeNumber1
        if J(i) <= nodeNumber1
            fprintf(outputFile,'%s  ',['"',sec1,'_',num2str(J(i)),'" -> "',sec1,'_',num2str(I(i)),'";  ']);
        else
            fprintf(outputFile,'%s  ',['"',sec2,'_',num2str(J(i)-nodeNumber1),'" -> "',sec1,'_',num2str(I(i)),'";  ']);
        end
    else
        if J(i) <= nodeNumber1
            fprintf(outputFile,'%s  ',['"',sec1,'_',num2str(J(i)),'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";  ']);
        else
            fprintf(outputFile,'%s  ',['"',sec2,'_',num2str(J(i)-nodeNumber1),'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";  ']);
        end
    end
    
    if mod(i,3) == 0
        fprintf(outputFile,'\n');
    end
end % End of MAIN LOOP

fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',sec1,' to ',sec2,' **********************']);
fclose(outputFile);

end


function connectOneWay(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon)

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


outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//****************** One Way Connect ',sec1,' to ',sec2,' ****************** ']);


sizeOut = nodeNumber1+nodeNumber2; 

tmp = ~eye(sizeOut); % makes a matrix of size sizeOut (both sections concatenated)
tmp(1:nodeNumber1,1:nodeNumber1) = 0; % disallows connections within a single section
tmp(1:sizeOut,(nodeNumber1+1):sizeOut) = 0;
i = find(tmp); % find all allowed connections
rp = randperm(length(i)); % shuffle connection order to randomize it
irp = i(rp);
connec = irp(1:numCon); % only take the number of connections desired (numCon)

[I J] = ind2sub([sizeOut, sizeOut], connec); %Convert linear indices to subscripts


% MAIN LOOP
for i = 1:numCon
    fprintf(outputFile,'%s  ',['"',sec1,'_',num2str(J(i)),'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";']);
    if mod(i,3) == 0
        fprintf(outputFile,'\n');
    end
end % End of MAIN LOOP

fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',sec1,' to ',sec2,' **********************']);
fclose(outputFile);

end




function connectAscending(file,numSeg,neuron1,neuron2,conLength)
%
% connectAscending(file,neuron1,neuron2,conLength)

outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//************ Ascending connect ',neuron1,' to ',neuron2,' for ',num2str(conLength),' sections *************** ']);
      
counter = 0;

  for sub = numSeg:-1:1

    if sub > conLength
      start = sub - 1;
      stop = sub - conLength;
    else
      start = sub - 1;
      stop = 1;
    end
    
    
    
    for i=start:-1:stop
          fprintf(outputFile,'  %s  ',['"',num2str(sub),neuron1,'" -> "',num2str(i),neuron2,'";']);
          
          counter = counter + 1;
          
          if mod(counter,4) == 0
            fprintf(outputFile,'\n');
          end    
    end

  end
  
fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',neuron1,' to ',neuron2,' **********************']);
fclose(outputFile);  % Close the file for printing

end




function connectDescending(file,numSeg,neuron1,neuron2,conLength)
%
% connectDescending(file,neuron1,neuron2,conLength)

outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//************ Descending connect ',neuron1,' to ',neuron2,' for ',num2str(conLength),' sections *************** ']);

  counter = 0; 
           
  for sub = 1:numSeg

    if (sub + conLength) < numSeg
      start = sub + 1;
      stop = sub + conLength;
    else
      start = sub + 1;
      stop = numSeg;
    end
    
 
    
    for i=start:stop
          fprintf(outputFile,'  %s  ',['"',num2str(sub),neuron1,'" -> "',num2str(i),neuron2,'";']);
          
          counter = counter + 1;
          
          if mod(counter,4) == 0
            fprintf(outputFile,'\n');
          end    
    end

  end
  
fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',neuron1,' to ',neuron2,' **********************']);
fclose(outputFile);  % Close the file for printing

end