function [adjacency, cellIndex] = makeRSdot2(file_to_save)


% Load Descending Neuron layout arrays
load('rsMap.mat')
load('dendriteMap.mat')


% Create beginning of Dot file

    %fileName = 'swtest_RS.dot';
    fileName = file_to_save;
        % File name for writing dot output
    outputFile = fopen(fileName,'a');
        % Assign a file Id. to fileName
    %fprintf(outputFile, 'subgraph cluster_RS {\n\n');
        % Initialize dot file
    %fclose(outputFile);
        % Close outputFile for writing.
        
% END beginning of Dot file


%*************************** Initialize Parameters *************************

tNC = 75000;
    % Total number of connections within RS
    
intRatio = 3;
    % ratio of interneurons to descending neurons
    
rsMap(:,:,3) = num2cell(intRatio*cell2mat(rsMap(:,:,2)));
    % Catalogues # of interneurons in each brain region. 
    for i=1:9
        for j=1:4
            if rsMap{i,j,3} == 0; rsMap{i,j,3} = 2 * intRatio;end
            % If no desc. cells exist in a section then interneurons equivalent
            % to 2 cell bodies still exist in that section.
        end
    end
    
intConDist = [ 0.20, 0.12, 0.14, 0.16, 0.18, 0.18];
    % intConDist is the distribution of connections among types of connections
    %      between mid/hindbrain interneurons.
    % intConDist = [Clustering, Mediolateral, Rostrocaudal, Contralateral,
    %      mixed Contralat./Rostroca., mixed Mediolat./Rostroca.]


totalInt = sum(sum(cell2mat(rsMap(:,:,3))));
    % total # of interneurons


totalDes = sum(sum(cell2mat(rsMap(:,:,2))));
    % Matrix of the number of Descending Neurons
    

% rsConDist is the distribution of connections within the
%      rs based on number of soma.
rsConDist = [ totalInt*2;
              (totalInt + totalDes);
              (totalDes + totalInt);
              totalDes;];
rsConDist = rsConDist ./ sum(rsConDist);


tIC = floor(rsConDist(1) * tNC); 
    % Interneuron interconnections
    
tI2DC = floor(rsConDist(2) * tNC); 
    % Interneuron to Descending Neuron connections
    
tD2IC = floor(rsConDist(3) * tNC); 
    % Descending Neuron to Interneuron connections
    
tDC = floor(rsConDist(4) * tNC); 
    % Descending Neuron interconnections

% The following ensures that connections are not lost during distribution.
%      Rounded off connections are assigned to interneuron connections.
distCon = sum([tIC, tI2DC, tD2IC, tDC]);
if (distCon < tNC); tIC = tIC + (tNC - distCon);end


    % Matlab Command Window Output --------------------------------------------
    fprintf(1, '\nThe total mid/hindbrain connections:  %g\n', tNC);
    fprintf(1, '     Connections among Interneurons:  %g\n', tIC);
    fprintf(1, '     Connections from Interneurons to Descending Neurons:  %g\n', tI2DC);
    fprintf(1, '     Connections from Descending Neurons to Interneurons:  %g\n', tD2IC);
    fprintf(1, '     Connections among Descending Neurons:  %g\n\n', tDC);
    fprintf(1, 'Check of total connections after distribution:  %g\n\n', sum([tIC, tI2DC, tD2IC, tDC]));
    % END Matlab CW Output ----------------------------------------------------


descendingNeurons = cell2mat(rsMap(:,:,2));

numSegments = cell2mat(dendriteMap(:,1,2));

numSoma = cell2mat(dendriteMap(:,1,3));  % Number of cell bodies each cell type has

%  It is possible for a call to ask for more connections than can
%  physically be accomodated.  In this case we overflow connections, which
%  later are added into interneuron clustering in the mid/hindbrain.
%  totalOverFlow and totalAssigned keep track of connections involved
%  during assignment proceedures.

% Number of overflow connections (defining variable)

totalOverFlow = 0;
%OverFlow = 0;
overflowed = 0;


% Number of connections trying to be assigned

totalAssigned = 0;
%connections = 0;
connected = 0;

totalRoundOff = 0;

rounded = 0;


%**************************************************************************
%**************************************************************************
%**************************** dot File contents **************************
%**************************************************************************
%**************************************************************************
%##########################################################################

% Create matrix for network
    
cellIndex = defineRS(fileName, rsMap, dendriteMap)';
    % Creates vector containing all cell bodies in the created Network
    % defineRS clusters RS neurons in the dot file
    
adjacency = sparse(length(cellIndex), length(cellIndex));
    % Initializes the adjacency matrix for the Network created


%**************************************************************************
% Specific Axon to Descending cell body connections
%**************************************************************************
pCon = 0.8;
    % Probability of connection to cell bodies of a given type/segment
    % References used: Segmental Homologies by Metcalfe 1986; Brain Neurons
    % which project to the spine by Kimmel 1982; T Reticular Interneurons
    % by Kimmel 1985
    
    makeSectionTypeBreakStart(fileName,'Axon Projections in Mid/Hindbrain');

% MiT axon connections in RS
[cond ovrd] = specCons2Cell(fileName, rsMap, dendriteMap,'Des_MiT_L', [{'Des_MiV1_L'}, {'Des_MiV2_L'}, {'Des_RoM3m_L'}, {'Des_RoM3l_L'}, {'Des_Mauthner_L'}, {'Des_CaV_L'}, {'Des_CaD_L'}, {'Des_T_L'}],pCon);
     connected = connected + cond;
     overflowed = overflowed + ovrd;
[cond ovrd] = specCons2Cell(fileName, rsMap, dendriteMap, 'Des_MiT_R', [{'Des_MiV1_R'}, {'Des_MiV2_R'}, {'Des_RoM3m_R'}, {'Des_RoM3l_R'}, {'Des_Mauthner_R'}, {'Des_CaV_R'}, {'Des_CaD_R'}, {'Des_T_R'}],pCon);
     connected = connected + cond;
     overflowed = overflowed + ovrd;
% T-cell axon connections in RS
[cond ovrd] = specCons2Sect(fileName, rsMap, dendriteMap, 'Des_T_L', [{'Mid_MR'},{'Ro1_MR'},{'Ro2_MR'},{'Ro3_MR'},{'Mi1_MR'},{'Mi2_MR'},{'Mi3_MR'},{'Ca_MR'}],pCon);
     connected = connected + cond;
     overflowed = overflowed + ovrd;
[cond ovrd] = specCons2Sect(fileName, rsMap, dendriteMap, 'Des_T_R', [{'Mid_ML'},{'Ro1_ML'},{'Ro2_ML'},{'Ro3_ML'},{'Mi1_ML'},{'Mi2_ML'},{'Mi3_ML'},{'Ca_ML'}],pCon);
     connected = connected + cond;
     overflowed = overflowed + ovrd;


    makeSectionTypeBreakStart(fileName,'Axon Projections in Mid/Hindbrain');


totalOverFlow = totalOverFlow + overflowed;
totalAssigned = totalAssigned + connected;

tDC = tDC - connected;

fprintf(1, 'Axon projections within the brain:\n');
fprintf(1, '     Assigned:  %g\n', connected);
fprintf(1, '     Overflowed:  %g\n', overflowed);
fprintf(1, '     Remaining connections among Des. Neurons:  %g\n\n', tDC);

%**************************************************************************
% END - Specific Axon to Descending cell body connections
%**************************************************************************
%##########################################################################

%##########################################################################
%**************************************************************************
% Start Descending neuron direct connections
%**************************************************************************

makeSectionTypeBreakStart(fileName,'Direct Descending Neurons')

overflowed = 0;
connected = 0;

numSegments = cell2mat(dendriteMap(:,1,2)); % The number of segments each 
                                            % descending type connects to

availCon = zeros(length(dendriteMap(:,1,2)), max(numSegments));
% Makes an empty matrix just as big as numSegments: number of descending 
% types by largest number of dendrite segment projections

    totalAvail = 0; % Available connections to be made in area

    % Map out the number of available connections in each segment and keep a
    % running total of the number of possible connections
    for i=1:length(numSoma)
        for k=1:numSegments(i)

            tmpSecName = [dendriteMap{i,k+1,1}, 'L']; % Using 'L' because RS
            tmpSecPos = dendriteMap{i,k+1,2};  % is symeterical

            [I J K] = ind2sub([9,4,4],find(strcmp(tmpSecName,rsMap)));

            if strcmp(tmpSecPos, 'ic')
                availCon(i, k) = (numSoma(i) * rsMap{I,J,2});
            else
                availCon(i, k) = (numSoma(i) * rsMap{I,J,2});
            end

            totalAvail = totalAvail + availCon(i, k);

        end
    end
    
    DDWeight = availCon;
    
%     totalAvail = sum(sum(availCon));
%     if totalAvail > tDC
%         weightCon = availCon/totalAvail;
%         weightCon = round(weightCon * tDC);
%     else
%         weightCon = availCon;
%     end
    
        pCon = 0.1;    
    
for i = 1:length(numSegments)
    
    for j = 1:numSegments(i)
        
        [I J K] = ind2sub([9,4,4],find(strcmp([dendriteMap{i,j+1,1},'L'],rsMap))); % Finds Section name in rsMap
        side = 'L';
        tmpNames = rsMap{I,J,4}; % Loads the desc. neuron types in the above section
        
        if strcmp(tmpNames(1), ''); break; end % Breaks if no desc. neurons in section

        [cond ovrd] = Sect2Soma(fileName,rsMap, dendriteMap, rsMap(I,J,1), dendriteMap{i,1,1}, side, availCon(i, j), pCon);
        connected = connected + cond;
        overflowed = overflowed + ovrd;
        
   end % j
end % i

for i = 1:length(numSegments)
    
    for j = 1:numSegments(i)
        
        [I J K] = ind2sub([9,4,4],find(strcmp([dendriteMap{i,j+1,1},'R'],rsMap))); % Finds Section name in rsMap
        side = 'R';
        tmpNames = rsMap{I,J,4}; % Loads the desc. neuron types in the above section
        
        if strcmp(tmpNames(1), ''); break; end % Breaks if no desc. neurons in section

        [cond ovrd] = Sect2Soma(fileName,rsMap, dendriteMap, rsMap(I,J,1), dendriteMap{i,1,1}, side, availCon(i, j), pCon);
        connected = connected + cond;
        overflowed = overflowed + ovrd;
        
   end % j
end % i 


totalOverFlow = totalOverFlow + overflowed;            

totalAssigned = totalAssigned + connected;

tDC = tDC - connected;

fprintf(1, 'Descending to Descending connections:\n');
fprintf(1, '     Assigned:  %g\n', connected);
fprintf(1, '     Overflowed:  %g\n', overflowed);
fprintf(1, '     Remaining connections among Des. Neurons:  %g\n\n', tDC);


makeSectionTypeBreakEnd(fileName,'Direct Descending Neurons')

%************************************************************
% End Descending neuron direct connections
%************************************************************


%**************************************************************************
% Start Interneruons connections to Descending Neurons
%**************************************************************************

makeSectionTypeBreakStart(fileName,'Interneurons to Descending Neurons');

overflowed = 0;
connected = 0;

availCon = zeros(length(dendriteMap(:,1,2)), max(numSegments));

totalAvail = 0; % Available connections to be made in area

% Map out the number of available connections in each segment and keep a
% running total of the number of possible connections
for i=1:length(numSoma)
    for k=1:numSegments(i)
        
        tmpSecName = [dendriteMap{i,k+1,1}, 'L'];
        tmpSecPos = dendriteMap{i,k+1,2};
        
        [I J K] = ind2sub([9,4,4],find(strcmp(tmpSecName,rsMap)));
        
        if strcmp(tmpSecPos, 'ic')
            availCon(i, k) = (numSoma(i) * rsMap{I,J,3})*4;
        else
            availCon(i, k) = (numSoma(i) * rsMap{I,J,3})*2;
        end

        totalAvail = totalAvail + availCon(i, k);

    end
end

% Scale number of connections and ensure no connections are lost
if totalAvail > tI2DC
   availCon = floor((availCon ./ sum(sum(availCon))) * tI2DC);
   diff = tI2DC - sum(sum(availCon));
    
    while diff > 0
        i = randi(length(numSoma));
        j = randi(max(numSegments));
        if availCon(i,j) ~= 0
            availCon(i,j) = availCon(i,j) + 1;
            diff = diff - 1;
        end
    end
end

I2DWeight = availCon;


for i = 1:length(availCon(:,1))
    
    for j =1:2
        
        numSec = numSegments(i);
        
        for k = 1:numSec
            
            tmpSecName = dendriteMap{i,k+1,1};
            tmpSecPos = dendriteMap{i,k+1,2};
            
            if j == 1
                    cellName = ['Des_',dendriteMap{i,1,1},'_L'];
                    if strcmp(tmpSecPos,'i')
                        secName = [tmpSecName,'L'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/2);
                        roundOff = availCon(i,k) - numCon*2;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                            
                    elseif strcmp(tmpSecPos,'c')
                        secName = [tmpSecName,'R'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/2);
                        roundOff = availCon(i,k) - numCon*2;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                    else
                        secName = [tmpSecName,'L'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/4);
                        roundOff = availCon(i,k) - numCon*4;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                            
                        secName = [tmpSecName,'R'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/4);
                        roundOff = availCon(i,k) - numCon*4;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                    end

            else
                    cellName = ['Des_',dendriteMap{i,1,1},'_R'];
                    if strcmp(tmpSecPos,'i')
                        secName = [tmpSecName,'R'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/2);
                        roundOff = availCon(i,k) - numCon*2;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                            
                    elseif strcmp(tmpSecPos,'c')
                        secName = [tmpSecName,'L'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/2);
                        roundOff = availCon(i,k) - numCon*2;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                            
                    else                        
                        secName = [tmpSecName,'L'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/4);
                        roundOff = availCon(i,k) - numCon*4;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                            
                        secName = [tmpSecName,'R'];
                        [I J K] = ind2sub([9,4,4],find(strcmp(secName,rsMap)));
                        numCon = floor(availCon(i,k)/4);
                        roundOff = availCon(i,k) - numCon*4;
                        [cond ovrd] = connectOneWay(fileName,secName,rsMap{I,J,3},cellName,dendriteMap{i,1,3},numCon);
                            connected = connected + cond;
                            overflowed = overflowed + ovrd;
                            totalRoundOff = totalRoundOff + roundOff;
                    end
            end
            
            %fprintf(1,'So far we have made %g connections.  We have made %g connections for %s.\n', [connected, (connected-lastSoFar), cellName]);
            %lastSoFar = connected;
        end
    end
end
                

totalOverFlow = totalOverFlow + overflowed;            

totalAssigned = totalAssigned + connected;

carryForward = sum(sum(availCon)) - connected;

tI2DC = tI2DC - connected;

fprintf(1, 'Interneuron to Descending Neuron connections:\n');
fprintf(1, '     Assigned:  %g\n', connected);
fprintf(1, '     Overflowed:  %g\n', overflowed);
fprintf(1, '     Rounded:  %g\n', carryForward);
fprintf(1, '     Remaining connections Int 2 Desc:  %g\n\n', tI2DC);


makeSectionTypeBreakEnd(fileName,'Interneurons to Descending Neurons');

%************************************************************
% End Interneruons to Descending Neurons connections
%************************************************************




%##########################################################################




%************************************************************
% Start Descending to Interneurons
%************************************************************
makeSectionTypeBreakStart(fileName,'Descending to Interneuron')

overflowed = 0;
connected = 0;

Des2IntWeight = round(descendingNeurons./sum(sum(descendingNeurons)) * tD2IC);

for i = 1:9
    for j = 1:4
        typesDescen = rsMap{i,j,4};
        
        for k = 1:length(typesDescen)
            
            if strcmp('',typesDescen{k}) == false
                
                if j == 1 || j == 2; cName = ['Des_',typesDescen{k},'_L']; end
                if j == 3 || j == 4; cName = ['Des_',typesDescen{k},'_R']; end
                
                cTypeNum = find(strcmp(typesDescen{k},dendriteMap));
                cTypeNum = dendriteMap{cTypeNum,1,3};
                
                furtherWeight = cTypeNum/rsMap{i,j,2};
                
                numConn = round(furtherWeight * Des2IntWeight(i,j));
                
                [cond ovrd] = connectOneWay(fileName,cName,cTypeNum,rsMap{i,j,1},rsMap{i,j,3},numConn);
                    connected = connected + cond;
                    overflowed = overflowed + ovrd;
                    
            end % If statement
        end % k loop
    end % j loop
end % i loop


fprintf(1,'Descending Neuron to Interneuron connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;            

totalAssigned = totalAssigned + connected;


makeSectionTypeBreakEnd(fileName,'Descending to Interneuron')
%************************************************************
% End Descending to Interneurons
%************************************************************

%##########################################################################
%     INTERNEURON CONNECTIONS
%##########################################################################

aIC = tIC + totalOverFlow + carryForward; % Ajusted total Interneuron connections

cons4Clust = floor(intConDist(1) * aIC);
    % Number of connections just for clustering
    
cons4MedLat = floor(intConDist(2) * aIC);
    % Number of connections just for Medial-Lateral Connections
    
cons4RoCad = floor(intConDist(3) * aIC);
    % Number of connections just for Rostral-Caudal Connections
    
cons4Contra = floor(intConDist(4) * aIC);
    % Number of connections just for Contralateral Connections
    
cons4ContraCross = floor(intConDist(5) * aIC);
    % Number of connections just for Contralateral Cross Connections
    
cons4Cross = floor(intConDist(6) * aIC);
    % Number of connections just for Cross Connections

% The following ensures that connections are not lost during distribution.
%      Rounded off connections are assigned to interneuron clustering.
distIntCon = sum([cons4Clust, cons4MedLat, cons4RoCad, cons4Contra, cons4ContraCross, cons4Cross]);
if (distIntCon < tIC); cons4Clust = cons4Clust + (tIC - distIntCon);end
distIntCon = sum([cons4Clust, cons4MedLat, cons4RoCad, cons4Contra, cons4ContraCross, cons4Cross]);

    % Matlab Command Window Output --------------------------------------------
    fprintf(1, 'Interneuron connection breakdown:\n');
    fprintf(1, '     Connections among Interneurons:  %g\n', tIC);
    fprintf(1, '     Connections unassigne in RS:  %g\n', (totalOverFlow + carryForward));
    fprintf(1, '     Adjusted Connections among Interneurons:  %g\n', aIC);
    fprintf(1, '     Clustering Connections:  %g\n', cons4Clust);
    fprintf(1, '     Mediolateral Connections:  %g\n', cons4MedLat);
    fprintf(1, '     Rostrocaudal Connections:  %g\n', cons4RoCad);
    fprintf(1, '     Contralateral Connections:  %g\n', cons4Contra);
    fprintf(1, '     Mixed RoCad/Contra Connections:  %g\n', cons4ContraCross);
    fprintf(1, '     Mixed MedLat/Contra Connections:  %g\n\n', cons4Cross);
    fprintf(1, 'Total interneuron connections after distribution:  %g\n\n', distIntCon);
    % END Matlab CW Output ----------------------------------------------------

    
    
%************************************************************
% Start of Medial-Lateral Connections
%************************************************************
makeSectionTypeBreakStart(fileName,'Medial-Lateral Connections');

overflowed = 0;
connected = 0;


% find the proper weight matrix
%MedLatWeight = [27 27; 15 15; 6 6; 7 7; 24 24; 9 9; 0 0; 0 0; 0 0];
MedLatWeight = zeros(9,2);
for i = 1:9
    j = 1;
    for k = 1:2
        %if descendingNeurons(i,j) == 0 || descendingNeurons(i,j+1) == 0
            %MedLatWeight(i,k) = 0;
        %else
            MedLatWeight(i,k) = cell2mat(rsMap(i,j,3))+cell2mat(rsMap(i,j+1,3));

        %end
        j = 3;
    end
end
MedLatSize = floor(MedLatWeight./sum(sum(MedLatWeight))*cons4MedLat);

% Basic check for consistent number of connections
if ( sum(sum(MedLatSize)) < cons4MedLat )
    diff = ( cons4MedLat - sum(sum(MedLatSize)) );
    MedLatSize(1,1) = MedLatSize(1,1) + floor(diff/2);
    MedLatSize(1,2) = MedLatSize(1,2) + ceil(diff/2);
end

for i = 1:9
    %connections = connections + MedLatSize(i,1) + MedLatSize(i,2);
        [cond ovrd] = connect2sections(fileName,rsMap{i,1,1},rsMap{i,1,3},rsMap{i,2,1},rsMap{i,2,3},MedLatSize(i,1));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
        [cond ovrd] = connect2sections(fileName,rsMap{i,3,1},rsMap{i,3,3},rsMap{i,4,1},rsMap{i,4,3},MedLatSize(i,2));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
end

fprintf(1,'Medial-Lateral connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;

totalAssigned = totalAssigned + connected;

makeSectionTypeBreakEnd(fileName,'Medial-Lateral Connections');


%************************************************************
% End of Medial-Lateral Connections
%************************************************************





%##########################################################################





%************************************************************
% Start Rostral-Caudal Connections
%************************************************************
makeSectionTypeBreakStart(fileName,'Rostral-Caudal Connections');

overflowed = 0;
connected = 0;

% Make proper weight matrix
%RoCadWeight = [27 15 15 27; 17 4 4 17; 5 8 8 5; 16 15 15 16; 17 16 16 17; 0 10 10 0; 0 6 6 0; 0 53 53 0];
RoCadWeight = zeros(8,4);
for i = 1:8
    for j = 1:4
       % if descendingNeurons(i,j)==0 || descendingNeurons(i+1,j)==0
          %  RoCadWeight(i,j) = 0;
       % else
            RoCadWeight(i,j) = cell2mat(rsMap(i,j,3))+cell2mat(rsMap(i+1,j,3));
       % end
    end
end
RoCadSize = round(RoCadWeight./sum(sum(RoCadWeight))*cons4RoCad);

% Basic check for consistent number of connections
if ( sum(sum(RoCadSize)) < cons4RoCad )
    diff = ( cons4RoCad - sum(sum(RoCadSize)) );
    RoCadSize(1,1) = RoCadSize(1,1) + ceil(floor(diff/2)/2);
    RoCadSize(1,2) = RoCadSize(1,2) + floor(floor(diff/2)/2);

    RoCadSize(1,3) = RoCadSize(1,3) + floor(ceil(diff/2)/2);
    RoCadSize(1,4) = RoCadSize(1,4) + ceil(ceil(diff/2)/2);
end

for i = 1:8
    for j = 1:4
        %connections = connections + RoCadSize(i,j);
        [cond ovrd] = connect2sections(fileName,rsMap{i,j,1},rsMap{i,j,3},rsMap{i+1,j,1},rsMap{i+1,j,3},RoCadSize(i,j));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
    end
end

fprintf(1,'Rostral-Caudal connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;
totalAssigned = totalAssigned + connected;

makeSectionTypeBreakEnd(fileName,'Rostral-Caudal Connections');


%************************************************************
% End Rostral-Caudal Connections
%************************************************************
%##########################################################################
%************************************************************
% Start Contralateral Connections
%************************************************************
makeSectionTypeBreakStart(fileName,'Contralateral Connections');

overflowed = 0;
connected = 0;

% Make proper weight matrix
%ContraWeight = [26; 4; 4; 12; 18; 14; 6; 6; 100];
ContraWeight = cell2mat(rsMap(:,2,3))+cell2mat(rsMap(:,3,3));
ContraSize = round(ContraWeight./sum(sum(ContraWeight))*cons4Contra);

%
% !!!!!!!!!!! ADD Check to Verify number of connections correct!!!!!!!
%

for i = 1:9
    %connections = connections + ContraSize(i);
    [cond ovrd] = connect2sections(fileName,rsMap{i,2,1},rsMap{i,2,3},rsMap{i,3,1},rsMap{i,3,3},ContraSize(i));
        connected = connected + cond;
        overflowed = overflowed + ovrd;
end

fprintf(1,'Contralateral connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;
totalAssigned = totalAssigned + connected;

makeSectionTypeBreakEnd(fileName,'Contralateral Connections');


%************************************************************
% End contralateral Connections
%************************************************************
%##########################################################################
%************************************************************
% Start Contralateral Cross Connections
%************************************************************
makeSectionTypeBreakStart(fileName,'Contralateral Cross Connections');

overflowed = 0;
connected = 0;

%ContraCrossWeight = [15 15; 4 4; 8 8; 15 15; 16 16; 10 10; 6 6; 53 53];
ContraCrossWeight = zeros(8,2);
for i = 1:8
    ContraCrossWeight(i,1) = cell2mat(rsMap(i,2,3))+cell2mat(rsMap(i+1,3,3));
    ContraCrossWeight(i,2) = cell2mat(rsMap(i,3,3))+cell2mat(rsMap(i+1,2,3));
end
ContraCrossSize = round(ContraCrossWeight./sum(sum(ContraCrossWeight))*cons4ContraCross);

%
% !!!!!!!!!!! ADD Check to Verify number of connections correct!!!!!!!
%

for i = 1:8
    %connections = connections + ContraCrossSize(i,1) + ContraCrossSize(i,2);
        [cond ovrd] = connect2sections(fileName,rsMap{i,2,1},rsMap{i,2,3},rsMap{i+1,3,1},rsMap{i+1,3,3},ContraCrossSize(i,1));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
        [cond ovrd] = connect2sections(fileName,rsMap{i,3,1},rsMap{i,3,3},rsMap{i+1,2,1},rsMap{i+1,2,3},ContraCrossSize(i,2));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
end

fprintf(1,'Contralateral Cross connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;

totalAssigned = totalAssigned + connected;

makeSectionTypeBreakEnd(fileName,'Contralateral Cross Connections');

%************************************************************
% End contralateral Cross Connections
%************************************************************
%##########################################################################
%************************************************************
% Start Contralateral Rostral-Caudal/Medial-Lateral (Cross)
%************************************************************
makeSectionTypeBreakStart(fileName,'Cross Connections');

overflowed = 0;
connected = 0;

% Make proper weight matrix
%CrossWeight = [16 26 26 16; 15 6 6 15; 10 3 3 10; 10 21 21 10; 22 11 11 22; 5 0 0 5; 0 0 0 0; 0 0 0 0];

CrossWeight = zeros(8,4);
for i = 1:8
    for j = 1:4
        if mod(j,2) == 1
            %if descendingNeurons(i,j) == 0 || descendingNeurons(i+1,j+1) == 0
            %    CrossWeight(i,j) = 0;
            %else
                CrossWeight(i,j) = cell2mat(rsMap(i,j,3))+cell2mat(rsMap(i,j+1,3));
            %end
        else
            %if descendingNeurons(i,j) == 0 || descendingNeurons(i+1,j-1) == 0
            %    CrossWeight(i,j) = 0;
            %else
                CrossWeight(i,j) = cell2mat(rsMap(i,j,3))+cell2mat(rsMap(i,j-1,3));
            %end
        end
    end
end
CrossSize = round(CrossWeight./sum(sum(CrossWeight))*cons4Cross);

%
% !!!!!!!!!!! ADD Check to Verify number of connections correct!!!!!!!
%

for i = 1:8
    for j = 1:4
        if mod(j,2) == 1
            %connections = connections + CrossSize(i,j);
            [cond ovrd] = connect2sections(fileName,rsMap{i,j,1},rsMap{i,j,3},rsMap{i+1,j+1,1},rsMap{i+1,j+1,3},CrossSize(i,j));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
        else
            %connections = connections + CrossSize(i,j);
            [cond ovrd] = connect2sections(fileName,rsMap{i,j,1},rsMap{i,j,3},rsMap{i+1,j-1,1},rsMap{i+1,j-1,3},CrossSize(i,j));
            connected = connected + cond;
            overflowed = overflowed + ovrd;
        end
    end
end

fprintf(1,'RoCad-MedLat Cross connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[overflowed connected])

totalOverFlow = totalOverFlow + overflowed;

totalAssigned = totalAssigned + connected;

makeSectionTypeBreakEnd(fileName,'Cross Connections');
%************************************************************
% End Contralateral Rostral-Caudal/Medial-Lateral (Cross)
%************************************************************
%##########################################################################
%*************************************************************
% Cluster interneuron in 36 brain regions
%*************************************************************
makeSectionTypeBreakStart(fileName,'Interneuron Clustering');

overflowed = 0;
connected = 0;

fprintf(1,'Total Connections for clustering: %g \n',[cons4Clust])

OverFlow = zeros(9,4);
    
clusterSize = round((cell2mat(rsMap(:,:,3)) ./ totalInt) * cons4Clust);
    % Matrix of number of connections within each brain segment
    
    %
    % !!!!! ADD check to verify number of connections correct!
    %

for i=1:4
    for j=1:9
        %connections = connections + clusterSize(j,i);
        [cond OverFlow(j,i)] = clusterInt(fileName,rsMap{j,i,1},rsMap{j,i,3}, clusterSize(j,i));
            connected = connected + cond;
            overflowed = overflowed + OverFlow(j,i);
    end
end

totalAssigned = totalAssigned + connected;
totalOverFlow = totalOverFlow + overflowed;

fprintf(1,'Clustering connections:\n     Overflowed:  %g\n     Assigned:  %g\n\n',[sum(sum(OverFlow)) connected]);
fprintf(1,'Average overflow for clustering: %g\n\n',mean(mean(OverFlow)));

makeSectionTypeBreakEnd(fileName,'Interneuron Clustering');
%************************************************************
% End of Clustering
%************************************************************

fprintf(1,'TOTAL CONNECTIONS ATTEMPTING TO BE ASSIGNED:  %g\n',tNC);
fprintf(1,'TOTAL CONNECTIONS SUCCESSFULLY ASSIGNED:  %g\n', totalAssigned);
fprintf(1,'TOTAL CONNECTIONS OVERFLOWED:  %g\n', overflowed);
fprintf(1,'TOTAL CONNECTIONS MISSING:  %g\n\n', (tNC - (totalAssigned + overflowed)));


%##########################################################################
%**************************************************************************
%**************************************************************************
%************************* End True File contents *************************
%**************************************************************************
%**************************************************************************


% Create end of Dot file
%outputFile = fopen(fileName,'a');
%fprintf(outputFile, '\n\n } \n\n');    
%fclose(outputFile);


end % END of makeRSdot function



%**************************************************************************
%
% ***************** BEGIN AUXILARY FUNCTIONS
%
%**************************************************************************



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



function cellIndex = defineRS(file,rmap, dmap)

outputFile = fopen(file,'a');  % Open the file for printing and write the header

cellIndex = cell(1,1);

counter = 1;

for i = 1:9
    switch i
        case 1; seg = 'Midbrain';
        case 2; seg = 'Ro1';
        case 3; seg = 'Ro2';
        case 4; seg = 'Ro3';
        case 5; seg = 'Mi1';
        case 6; seg = 'Mi2';
        case 7; seg = 'Mi3';
        case 8; seg = 'Ca';
        case 9; seg = 'T';
    end
    
    fprintf(outputFile, '%s \n\n', ['subgraph cluster_', seg,'  {']);
    
    for j = 1:4

        % Start dot file subgraph cluster
        if (rmap{i,j,2} + rmap{i,j,3}) ~= 0
            fprintf(outputFile,'%s \n\n',['subgraph cluster_',rmap{i,j,1},'  {']);
        end
        cellsInSeg = rmap{i,j,4}; % List of cell types in (i,j) segment
        
        % Descending neurons
        for k = 1:length(rmap{i,j,4})
            
            cType = cellsInSeg{k};
            
            if strcmp(char(cType), ''); break; end
          
            %  Find the number of cell bodies of cType
            for l = 1:length(dmap(:,1,3))
                if strcmp(cType, dmap{l,1,1})
                    cCount = dmap{l,1,3}; 
                    break
                end
            end
            
            % Declaration of Descending neurons
            for m = 1:cCount 
                if ismember(j, [1,2]); side = 'L';
                else side = 'R';
                end                
                fprintf(outputFile,'%s',['"Des_',char(cType),'_',side,'_',num2str(m),'";  ']);
                
                cellIndex(counter) = {['Des_',char(cType),'_',side,'_',num2str(m)]};
                
                if mod(counter,6) == 0
                    fprintf(outputFile,'\n');
                end
                
                counter = counter + 1;
            end % End of declaration of Descending neurons (m)
            
        end % End Descending neuron loop (k)
        
        
        % Declaration of Interneurons
        for n = 1:rmap{i,j,3}
            
            seg = rmap{i,j,1};
            
            fprintf(outputFile,'%s',['"',seg,'_',num2str(n),'";  ']);
            
            cellIndex(counter) = {[seg,'_',num2str(n)]};
            
            if mod(counter,6) == 0
                    fprintf(outputFile,'\n');
            end
            
            counter = counter + 1;
            
        end % End of Interneuron loop (n)
        
        if (rmap{i,j,2} + rmap{i,j,3}) ~= 0
            fprintf(outputFile,'\n       } \n\n');
        end
        
    end % End of column loop (j)
    
    fprintf(outputFile,'\n       } \n\n');
  
end % End of row loop (i)

fclose(outputFile);

end



%**************************************************************************
%
% ***************** CONNECTION ASSIGNMENT
%
%**************************************************************************




function [connected ovrflw] = clusterInt(file,seg,numInt,numCon)
% This is a process adapted from Olaf Sporns makerandCIJ_dir.m from the
% brain connectivity toolbox.

ovrflw = 0;
connected = 0;

% open file and print header
outputFile = fopen(file,'a');

fprintf(outputFile,'\n\n %s \n\n',['//************************ Interneuron clustering for ',seg,' ********************** ']);

% randomly create the connections we are going to use
tmp = ~eye(numInt);
i = find(tmp);
rp = randperm(length(i));
irp = i(rp);
if sum(sum(tmp)) < numCon
    ovrflw = numCon - sum(sum(tmp));
    connected = sum(sum(tmp));
else
    connected = numCon;
end
i = irp(1:connected);

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





function [connected ovrflw] = connect2sections(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon)
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

ovrflw = 0;
connected = 0;

outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//************************ ',sec1,' to ',sec2,' ********************** ']);


sizeOut = nodeNumber1+nodeNumber2; 

tmp = ~eye(sizeOut); % makes a matrix of size sizeOut (both sections concatenated)
tmp(1:nodeNumber1,1:nodeNumber1) = 0; % disallows connections within a single section
tmp((nodeNumber1+1):sizeOut,(nodeNumber1+1):sizeOut) = 0;
i = find(tmp); % find all allowed connections
rp = randperm(length(i)); % shuffle connection order to randomize it
irp = i(rp);
if sum(sum(tmp)) < numCon
    ovrflw = numCon - sum(sum(tmp));
    connected = sum(sum(tmp));
else
    connected = numCon;
end
connec = irp(1:connected); % only take the number of connections desired (numCon)

[I J] = ind2sub([sizeOut, sizeOut], connec); %Convert linear indices to subscripts


% MAIN LOOP
for i = 1:connected
    if J(i) <= nodeNumber1
        fprintf(outputFile,'%s  ',['"',sec1,'_',num2str(J(i)),'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";']);
    else
        fprintf(outputFile,'%s  ',['"',sec2,'_',num2str(J(i)-nodeNumber1),'" -> "',sec1,'_',num2str(I(i)),'";']);
    end
    
    if mod(i,3) == 0
        fprintf(outputFile,'\n');
    end
end % End of MAIN LOOP

fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',sec1,' to ',sec2,' **********************']);
fclose(outputFile);

end


function [connected ovrflw] = connectOneWay(file,sec1,nodeNumber1,sec2,nodeNumber2,numCon)

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

ovrflw = 0;
connected = 0;

outputFile = fopen(file,'a');  % Open the file for printing and write the header
fprintf(outputFile,'\n\n %s \n\n',['//****************** One Way Connect ',sec1,' to ',sec2,' ****************** ']);


sizeOut = nodeNumber1+nodeNumber2; 

tmp = ~eye(sizeOut); % makes a matrix of size sizeOut (both sections concatenated)
tmp(1:nodeNumber1,1:nodeNumber1) = 0; % disallows connections within a single section
tmp(1:sizeOut,(nodeNumber1+1):sizeOut) = 0;
i = find(tmp); % find all allowed connections
rp = randperm(length(i)); % shuffle connection order to randomize it
irp = i(rp);
if sum(sum(tmp)) < numCon
    ovrflw = numCon - sum(sum(tmp));
    connected = sum(sum(tmp));
else
    connected = numCon;
end
connec = irp(1:connected); % only take the number of connections desired that do not overflow

[I J] = ind2sub([sizeOut, sizeOut], connec); %Convert linear indices to subscripts


% MAIN LOOP
for i = 1:connected
    fprintf(outputFile,'%s  ',['"',sec1,'_',num2str(J(i)),'" -> "',sec2,'_',num2str(I(i)-nodeNumber1),'";']);
    if mod(i,3) == 0
        fprintf(outputFile,'\n');
    end
end % End of MAIN LOOP

fprintf(outputFile,'\n\n %s \n\n',['//************************ END of ',sec1,' to ',sec2,' **********************']);
fclose(outputFile);

end



function [connected ovrflw] = specCons2Cell(file,rMap, dMap, fromCell, to, Kden)

    ovrflw = 0;
    connected = 0;
    
    fCell = fromCell;
    fromCell = regexp(fromCell, '_', 'split');
    fromCell = fromCell{2};
    
    numSomaF = getNumSoma(dMap, fromCell);
    
    for i = 1:length(to)
     
        numSomaT = getNumSoma(dMap, to{i});
        
        avgSoma = (numSomaF + numSomaT) / 2;
        numCon = round( Kden * avgSoma);
        if numCon == 0; numCon = 1; end
    
        [cond ovrd] = connectOneWay(file, fCell, numSomaF, to{i}, numSomaT, numCon);
        
        connected = connected + cond;
        ovrflw = ovrflw + ovrd;
            
    end
        
    
end


function [connected ovrflw] = specCons2Sect(file,rMap, dMap, fromCell, toSect, Kden)

    ovrflw = 0;
    connected = 0;
    
    fCell = fromCell;
    fromCell = regexp(fromCell, '_', 'split');
    fromCell = fromCell{2};
    
    numSomaF = getNumSoma(dMap, fromCell);
    
    
    for i = 1:length(toSect)
            
        %neurons = getNeurInSect(rMap, toSect{i});
        
                toSectLoc = regexp(toSect{i}, '_', 'split');
                toSectLoc = toSectLoc(2);
                
                [x y] = ind2sub([9 4],find(strcmp( toSect{i}, rMap(:,:,1))));
                
                types = rMap{x,y,4};
                
            for j=1:length(types)
                
                cType = types{j};    
                if strcmp(toSectLoc,'LL') | strcmp(toSectLoc,'ML')
                    cType = ['Des_',cType,'_L'];
                else
                    cType = ['Des_',cType,'_R'];
                end
                    
                    
                numSomaT = getNumSoma(dMap, cType);
                avgSoma = (numSomaF + numSomaT) / 2;    
                    
                numCon = round( Kden * avgSoma);
                if numCon == 0; numCon = 1; end
                
                [cond ovrd] = connectOneWay(file, fCell, numSomaF, cType, numSomaT, numCon);
                
                connected = connected + cond;
                ovrflw = ovrflw + ovrd;
            end
    end
        
    
end


function [connected ovrflw] = Sect2Soma(file,rMap, dMap, fromSect, toCell, side, numCon, pCon)

    ovrflw = 0;
    connected = 0;
    
[tx ty] = ind2sub(size(dMap(:,:,1)),find(strcmp( toCell, dMap(:,:,1))));
toCell = ['Des_', toCell, '_', side];


    
    numSomaT = getNumSoma(dMap, toCell);
    
            
        %neurons = getNeurInSect(rMap, toSect{i});
        
%                 fromSectLoc = regexp(fromSect{1}, '_', 'split');
%                 fromSectLoc = fromSectLoc{2};
           
                [x y] = ind2sub([9 4],find(strcmp( fromSect, rMap(:,:,1))));
                
                types = rMap{x,y,4};
                somas = zeros(size(types));
                
                for t=1:length(types)
                    [m n] = ind2sub(size(dMap(:,:,1)),find(strcmp( types{1}, dMap(:,:,1))));
                    somas(t) = dMap{m, n, 3};
                end
                
            for j=1:length(types)

                if strcmp(dMap{tx, j+1, 2}, 'c')
                    if side == 'R'
                        side = 'L';
                    else
                        side = 'R';
                    end
                elseif strcmp(dMap{tx, j+1, 2}, 'ic')
                    side = 'B';
                end
                
                cType = types{j};
                numSomaF = getNumSoma(dMap, cType);
                if side == 'B'
                    tmpCtype = ['Des_',cType,'_L'];
                    if ~strcmp(tmpCtype, toCell)
                    [cond ovrd] = connectOneWay(file, tmpCtype, numSomaF, toCell, numSomaT, round(pCon * numCon));
                    connected = connected + cond;
                    ovrflw = ovrflw + ovrd;
                    end
                    
                    cType = ['Des_',cType,'_R'];
                    if ~strcmp(tmpCtype, toCell)
                    [cond ovrd] = connectOneWay(file, cType, numSomaF, toCell, numSomaT, round(pCon * numCon));
                    connected = connected + cond;
                    ovrflw = ovrflw + ovrd;
                    end
                else
                    cType = ['Des_',cType,'_',side];
                    if ~strcmp(cType, toCell)
                    [cond ovrd] = connectOneWay(file, cType, numSomaF, toCell, numSomaT, round(pCon * numCon));
                    
                    connected = connected + cond;
                    ovrflw = ovrflw + ovrd;
                    end
                end

            end
        
end


%**************************************************************************
%
% ***************** INFORMATION RETRIEVAL
%
%**************************************************************************


function dMapLoc = getDMapLoc(dMap, cType)

    cType = regexp(cType, '_', 'split');
    if length(cType) >2;    cType = cType{2};
    else cType = cType{1};
    end
        % Ignores irrelevant '_L/R' distinction and 'Des_'
    
    dMapLoc = find(strcmp(dMap(:,1,1), cType));
        % Finds cType in dendriteMap registry

end

function neurInSect = getNeurInSect(rMap, sect)

rMapLoc = find(strcmp(rMap(:,:,1), sect));

[i j] = ind2sub([9 4], rMapLoc);

neurInSect = rMap(i,j,2);

end

function numSoma = getNumSoma(dMap, cType)

    dMapLoc = getDMapLoc(dMap, cType);
        % Finds position in dendriteMap
    numSoma = dMap{dMapLoc, 1, 3};
        % Returns associated number of soma for the cell type
end


function sectInfo = getSectInfo(dMap, cType)

    dMapLoc = getDMapLoc(dMap, cType);

    sectInfo = cell(3, dMap{dMapLoc, 1, 2});

for i = 1:dMap{dMapLoc, 1, 2};

    sectInfo(1,i) = {dMap{dMapLoc, i+1, 1}};
    sectInfo(2,i) = {dMap{dMapLoc, i+1, 2}};
    sectInfo(3,i) = {dMap{dMapLoc, i+1, 3}};

end

end


