function adj = updateMatrixDed(adj,names)


numNeuroConnec = [1 1 1 1 1 2 2 2 2 3 3 3 4 4 5];



sensToRS =   [   1     2     3     4     0     0     0;
                 1     2     3     4     5     6     0;
                 3     4     5     6     7     8     0;
                 5     6     7     8     9    10    11;
                 7     8     9    10    11    12    13;
                 9    10    11    12    13    14    15;
                12    13    14    15    16    17     0;
                14    15    16    17    18    19    20;
                16    17    18    19    20     0     0];

 

rsSec =     {'Mid_M',    'Mid_L';
             'Ro1_M',    'Ro1_L';
             'Ro2_M',    'Ro2_L';
             'Ro3_M',    'Ro3_L';
             'Mi1_M',    'Mi1_L';
             'Mi2_M',    'Mi2_L';
             'Mi3_M',    'Mi3_L';
             'Ca_M',     'Ca_L' ;
             'T_M',      'T_L'  };
            
            
 numConsPerCell = [16 16 19 19 10 10 7 7 8 8 8 10 10 6 6 12 12 11 11 11]';
 
 newMaxConnec = 25;
 
 
            
            
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%      Sensory Cells       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           

for s=1:2
    
    
    if s == 1
        side = 'L';
    else
        side = 'R';
    end
    
    
    for sens = 1:20
        
        [X Y] = find(sens == sensToRS);
        
        rsIntInd = [];
        
        counter = 0;
        
        for i = 1:length(X)
            
            tmpInd = [];
            
            rsIntTypes = [rsSec{X(i),1},side];
            tmpLengths = length(rsIntTypes);
            tmpInd = find(strncmp(names,rsIntTypes,tmpLengths));
            
            currentLength = length(rsIntInd)+1;
            
            rsIntInd(currentLength:(currentLength+length(tmpInd)-1)) = tmpInd;
            
            tmpInd = [];
            
            rsIntTypes = [rsSec{X(i),2},side];
            tmpLengths = length(rsIntTypes);
            tmpInd = find(strncmp(names,rsIntTypes,tmpLengths));
            
            currentLength = length(rsIntInd)+1;
            
            rsIntInd(currentLength:(currentLength+length(tmpInd)-1)) = tmpInd;
            
        end
        
        
        rsIntInd = rsIntInd(randperm(length(rsIntInd)));
        
        
        numConnec = newMaxConnec - numConsPerCell(sens);
        
        tmpSens = find(strcmp(names,['Sens_',side,num2str(sens)]));
        
        for j = 1:numConnec

            adj(tmpSens,rsIntInd(j)) = 1;
        
        end
        
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%      Neuromasts Only       %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

for s = 1:2
    
    if s == 1
        side = 'L';
    else
        side = 'R';
    end
    
    
    numNeuroConnec = numNeuroConnec(randperm(length(numNeuroConnec)));
    indSens = find(strncmpi(names,['Sens_',side],6));
    
    
    
    for i = 1:15
        ind = [];
        ind = find(strcmp(['Neuromast_',side,num2str(i)],names));
        
        tempInd = find(adj(ind,:) ==1);
        
        adj(ind,tempInd) = 0;
        
        num = numNeuroConnec(i);        
        
        indSensPerm = indSens(randperm(length(indSens)));
        
        adj(ind,indSensPerm(1:num)) = 1;
        
        
    end
    
    
    
end

%}


adj(1:70,1:70) = 0;  % Zero out the Neuromast to Sensory connections

for s = 1:2
    
    if s == 1
        side = 'L';
    else
        side = 'R';
    end
    
    
    numNeuroRemain = numNeuroConnec(randperm(length(numNeuroConnec)));
    
    
    
    for sens = 1:20
        
        indSens = find(strcmp(names,['Sens_',side,num2str(sens)]));
        startP = find(numNeuroRemain > 0);
        
        if isempty(startP)
            
            ind = [];
            ind = find(strcmp(['Neuromast_',side,num2str(15)],names));
            
            adj(ind,indSens) = 1;
            
            indSens = find(strcmp(names,['Sens_',side,num2str(sens-1)]));
            
            adj(ind,indSens) = 0;
            
        else
            startP = startP(1);
        end

        
        if startP == 15 |  numNeuroRemain(startP+1) == 0
            
            ind = [];
            ind = find(strcmp(['Neuromast_',side,num2str(startP)],names));
            
            if adj(ind,indSens) == 1
                disp('ERROR! Tried connecting something that is already connected!')
            else
                adj(ind,indSens) = 1;
            end
            
            numNeuroRemain(startP) = numNeuroRemain(startP) - 1;
            
        elseif ~isempty(startP)
            
            ind = [];
            ind = find(strcmp(['Neuromast_',side,num2str(startP)],names));
            
            if adj(ind,indSens) == 1
                disp('ERROR! Tried connecting something that is already connected!')
            else
                adj(ind,indSens) = 1;
            end
            
            numNeuroRemain(startP) = numNeuroRemain(startP) - 1;
            
            ind = [];
            ind = find(strcmp(['Neuromast_',side,num2str(startP+1)],names));
            
            if adj(ind,indSens) == 1
                disp('ERROR! Tried connecting something that is already connected!')
            else
                adj(ind,indSens) = 1;
            end
            
            numNeuroRemain(startP+1) = numNeuroRemain(startP+1) - 1;
            
        end % End of startP If statement
        
    end % End of sens For Loop
    
    
    
    if sum(numNeuroRemain) == 0
        disp('Everything Connected!')
    else
        disp('Did NOT connect everything!')
        
        indSens = find(strncmpi(names,['Sens_',side],6));
        
        for i = 1:nnz(numNeuroRemain)
            
           neuromast = find(numNeuroRemain > 0);
           neuromast = neuromast(1);
           
           ind = [];
           ind = find(strcmp(['Neuromast_',side,num2str(neuromast)],names));
           
           for cons = 1:numNeuroRemain(neuromast)
               
              tempSensInds = find(adj(ind,:)>0);
              
              if isempty(tempSensInds)
                  
                  tempSensInds = find(strcmp(names,['Sens_',side,num2str(neuromast)]));
                  adj(ind,tempSensInds) = 1;
                  numNeuroRemain(neuromast) = numNeuroRemain(neuromast) - 1;
              
              elseif find(indSens == (tempSensInds(end)+1))
                  
                  adj(ind,tempSensInds(end)+1) = 1;
                  numNeuroRemain(neuromast) = numNeuroRemain(neuromast) - 1;
                  
              elseif find(indSens == (tempSensInds(1)-1))
                  
                  adj(ind,tempSensInds(1)-1) = 1;
                  numNeuroRemain(neuromast) = numNeuroRemain(neuromast) - 1;
                  
              end
              
               
           end
            
        end
        
    end
    
end
        
            




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%  Make RS the right size  %%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


G = adj(71:1196+70,71:1196+70);

connRMV = nnz(G) - 70000;

Tind = find(strncmp(names,'T_',2));


if connRMV > 0
    
    [I J] = find(adj(Tind,Tind) == 1);
    
    Tconns = sub2ind([2616 2616], Tind(I),Tind(J));
    
    Tconns = Tconns(randperm(length(Tconns)));
    
    adj(Tconns(1:connRMV)) = 0;
    
elseif connRMV < 0
    
    T = adj + eye(length(adj));
    
    connRMV = - connRMV;
    
    [I J] = find(T(Tind,Tind) == 0);
    
    Tconns = sub2ind([2616 2616], Tind(I),Tind(J));
    
    Tconns = Tconns(randperm(length(Tconns)));
    
    adj(Tconns(1:connRMV)) = 1;
    
end
        
    
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        

