function [G names] = makeSMP(file,dedOrDIST)
% makeSMP(file,dedOrDistLL,dedOrDistSpine)

% 1 is DEDICATED
% 0 is DISTRIBUTED

if nargin < 2
    dedOrDIST = 1;
end

tic;

fprintf(1,'\n \n Starting Lateral Line ..............')

begFile = fopen('lateral_line_stoch.dot');
lateral_line = fscanf(begFile,'%c');

outputFile = fopen(file,'w');
fprintf(outputFile,lateral_line);
% fclose(outputFile);

fprintf(1,'............. Finished Lateral Line \n \n ')


fprintf(1,'\n \n Starting the RS ..............')

makeRSdot2(file)

fprintf(1,'..............  Finished RS')


fprintf(1,'\n \n Starting LL to RS ..............')

LLtoRS(file,dedOrDIST);

fprintf(1,'............. Finished')


fprintf(1,'\n \n Starting RS to Spine ..............')

RS2Spine3(file,dedOrDIST)

fprintf(1,'............. Finished')


fprintf(1,'\n \n Starting Spine ..............')

makeSpinedot(file)

fprintf(1,'........... Finished')


fprintf(1,'\n \n Making Adjacency Matrix ...........')

myaxes = axes; [G names] = drawDot(file,myaxes,'LayoutEngine','twopi');

fprintf(1,'......... Finished')


fprintf(1,'\n \n Updating Adjacency Matrix.........')

G = updateMatrixDed(G,names);

fprintf(1,'......... Finished \n\n')

toc;

% outputFile = fopen(file,'a');
% endFile = fopen('final_rs.dot');
% RS_section = fscanf(endFile,'%c');


