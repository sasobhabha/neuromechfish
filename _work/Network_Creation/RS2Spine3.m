function RS2Spine3(file,dedORdist)

load dendriteMap
load SpineList

MasterConnecCounter = 0;
%MasterKCounter = 0;
kcounter = 1;


fileName = file;
% File name for writing dot output
outputFile = fopen(fileName,'a');
% Assign a file Id. to fileName
fprintf(outputFile,'\n//*************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'//                   Starting Section RS to Spine connections \n');
fprintf(outputFile,'//**************************************************************************************************************************\n');
fprintf(outputFile,'//**************************************************************************************************************************\n\n');
% Initialize dot file



DesNeuronList = dendriteMap(:,1,1);
DesNeuronList(:,2) = dendriteMap(:,1,3);
DesNeuronList(:,5) = dendriteMap(:,2,4);
DesNeuronList(:,6) = dendriteMap(:,3,4);

%Take care of the long descending neurons
temp = dendriteMap(:,1,4);
temp = find(cell2mat(temp)>22);

for i = 1:length(temp)
    DesNeuronList(temp(i),3) = {'U'};
    DesNeuronList(temp(i),4) = {25};
end

DesNeuronList = DesNeuronList(randperm(length(DesNeuronList)),:);


if dedORdist == 1 %Dedicated model
    
    counter1 = 0; %Count the number of cell bodies
    counter2 = 1; %Count the number of cell types
    while counter1 < 29   % This is going to pick which cells are to be dist type
        temp = DesNeuronList{counter2,2};
        if (temp + counter1) < 30 & isempty(DesNeuronList{counter2,4})
            counter1 = counter1 + temp; %Cell bodies will fit and are thus assigned
            DesNeuronList(counter2,3) = {'U'}; %Assigning "uniform distribution"
            
            randDiceRoll = rand(1);
            if randDiceRoll > .5
                DesNeuronList(counter2,4) = {10};  %Assign length of connections
            elseif randDiceRoll > .15
                DesNeuronList(counter2,4) = {15};  %Assign length of connections
            elseif randDiceRoll > .05 
                DesNeuronList(counter2,4) = {20};  %Assign length of connections
            else
                DesNeuronList(counter2,4) = {25};  %Assign length of connections
            end
            
        elseif isempty(DesNeuronList{counter2,4})
            DesNeuronList(counter2,3) = {'T'}; %Assigning "Tapered distribution"            
        end
        counter2 = counter2 + 1; %Cell bodies will NOT fit and are thus skipped
    end
    
    for i = counter2:length(DesNeuronList)
        if isempty(DesNeuronList{i,4})
            DesNeuronList(i,3) = {'T'};
        end
    end
    
    
else
    
    
    counter1 = 0; %Count the number of cell bodies
    counter2 = 1; %Count the number of cell types
    while counter1 < 29   % This is going to pick which cells are to be dist type
            temp = DesNeuronList{counter2,2};
        if (temp + counter1) < 30 & isempty(DesNeuronList{counter2,4})
            counter1 = counter1 + temp; %Cell bodies will fit and are thus assigned
            DesNeuronList(counter2,3) = {'T'}; %Assigning "taper distribution"
        elseif isempty(DesNeuronList{counter2,4})
            DesNeuronList(counter2,3) = {'U'}; %Assigning "Uniform distribution"
            
            randDiceRoll = rand(1);
            if randDiceRoll > .5
                DesNeuronList(counter2,4) = {10};  %Assign length of connections
            elseif randDiceRoll > .15
                DesNeuronList(counter2,4) = {15};  %Assign length of connections
            elseif randDiceRoll > .05 
                DesNeuronList(counter2,4) = {20};  %Assign length of connections
            else
                DesNeuronList(counter2,4) = {25};  %Assign length of connections
            end
            
        end
        counter2 = counter2 + 1; %Cell bodies will NOT fit and are thus skipped
    end
    
    for i = counter2:length(DesNeuronList)
        if isempty(DesNeuronList{i,4})
            DesNeuronList(i,3) = {'U'}; %Assigning "Uniform distribution"
            
            randDiceRoll = rand(1);
            if randDiceRoll > .5
                DesNeuronList(i,4) = {10};  %Assign length of connections
            elseif randDiceRoll > .15
                DesNeuronList(i,4) = {15};  %Assign length of connections
            elseif randDiceRoll > .05 
                DesNeuronList(i,4) = {20};  %Assign length of connections
            else
                DesNeuronList(i,4) = {25};  %Assign length of connections
            end
            
            
        end
    end
    
end


% MAIN LOOPS


for s = 1:2

    if s == 1
        NeuronSide = 'L';
    else
        NeuronSide = 'R';
    end
    
    
for i = 1:length(DesNeuronList)
    for j = 1:DesNeuronList{i,2}
        
        
        % #####################  Modified From Original ###################
        % -New RS2Spine, everything starts at 1
        % -Uniform connection has different possible lengths
        %##################################################################
        % Getting the connections for every segment that the neuron
        % projects to in the spine
        if DesNeuronList{i,3} == 'U'  % the descending neuron is uniform
            segs = DesNeuronList{i,4};
            if segs == 25 % The descending neuron goes all the way
                connec = 4*ones(25,1);
                temp = zeros(25,1);
                x = randperm(25);
                temp(x(1:20)) = 1;
                connec = connec + temp;
                startplace = 1;
            else % the neuron only goes 10 segments
                %connec = 12*ones(10,1);
                %startplace = randi(16);
                startplace = 1;
                
                if segs == 20
                    connec = 6*ones(20,1);
                elseif segs == 15
                    connec = 8*ones(15,1);
                else
                    connec = 12*ones(10,1);
                end
                
            end
            
        else  % the descending neuron is tapered
            connec = getGamma;
            while max(connec) > 27 | length(connec) > 25
                connec = getGamma;
            end
            segs = length(connec);
            %startplace = randi(25-segs+1);
            startplace = 1;
        end

    
        
        
        
        
        if DesNeuronList{i,5} == 'c'
            switch NeuronSide
                case 'L'
                    axonSide = 'R';
                case 'R'
                    axonSide = 'L';
            end
        else
            switch NeuronSide
                case 'R'
                    axonSide = 'R';
                case 'L'
                    axonSide = 'L';
            end
        end
        
        if     (DesNeuronList{i,6} == 'c' & axonSide == 'L')
            colatSide = 'R';
        elseif (DesNeuronList{i,6} == 'i' & axonSide == 'L')
            colatSide = 'L';
        elseif (DesNeuronList{i,6} == 'c' & axonSide == 'R')
            colatSide = 'L';
        elseif (DesNeuronList{i,6} == 'i' & axonSide == 'R')
            colatSide = 'R';
        else
            colatSide = 'B';
        end
        
%         if (MasterKCounter == MasterKCounter + kcounter -1 & MasterKCounter ~= 0)
%             error('MasterKCounter didn''t change!');
%         end
        
        %MasterKCounter = MasterKCounter + kcounter - 1;
        kcounter = 1;
        %disp(MasterKCounter)
        
        for k = startplace:(startplace+segs-1)  % Go through all segments this neuron connects to in spine
            
            if k == startplace
                fprintf(outputFile,'\n\n\n %s \n\n\n',['//Connecting the ',DesNeuronList{i,1},':  Cell Body# ',num2str(j)]);
            end
            
            if colatSide == 'R'
                tempSpine = spineListRight(randperm(27));
            elseif colatSide == 'L'
                tempSpine = spineListLeft(randperm(27));
            else
                tempSpine = [spineListRight;spineListLeft];
                tempSpine = tempSpine(randperm(length(tempSpine)));
            end
        
             %disp(sum(connec))
             for m=1:connec(kcounter) % This is for each individual segment connections
                 
                 fprintf(outputFile,'%s  ',['"Des_',DesNeuronList{i,1},'_',NeuronSide,'_',num2str(j),'" -> "',num2str(k),tempSpine{m},'";']);
                 MasterConnecCounter = MasterConnecCounter + 1;
                 if mod(m,3) ==0
                     fprintf(outputFile,'\n');
                 end
             end
             
             kcounter = kcounter + 1;
        end
    end
end
                 
end

            
            
            
   fprintf(outputFile,'\n\n\n\n');
   MasterConnecCounter;         
           
  end      
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
    
    