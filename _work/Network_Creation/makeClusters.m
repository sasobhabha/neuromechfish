function makeClusters(type,side,N)

count = 0;
    

fprintf(1,'\n\n subgraph cluster_%s { \n\n',[cat(2,cat(2,type,'_'),side)])

if nargin<3
    fprintf(1,'"%s";  ',[  cat(2,type,cat(2,'_',side))  ])
else
    for i=1:N
        fprintf(1,'"%s";  ',[  cat(2,cat(2,type,cat(2,'_',side)),num2str(i))  ])
        count = count+1;
        if mod(count,5) == 0
            fprintf(1,'\n')
        end
    end
end

fprintf(1,'\n\n } \n\n\n')

