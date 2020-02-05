% Script for plotting counties
S = shaperead('County.shp');
l = length(S);
% Sort so in the same order as data
[~,index] = sortrows({S.COUNTY}.'); S = S(index); clear index

%% Give the struct its Infecteds timeseries

timestep = 15;

for i = 1:l
    S(i).I = double(Epidemic(i,2,timestep) + Epidemic(i,5,timestep));
end


%%
%Make symbol specification for the polygon geometry
     ColoredConstituencies = makesymbolspec('Polygon',...
         {'I',[0 10],'FaceColor',[0.7 0.7 0.7]},...
         {'I',[11 100],'FaceColor',[1 0 0]*0.1 + [0 0 1]*0.9 },...
         {'I',[101 500],'FaceColor',[1 0 0]*0.25 + [0 0 1]*0.75 },...      
         {'I',[501 1000],'FaceColor',[1 0 0]*0.5 + [0 0 1]*0.5 },...              
         {'I',[1001 5000],'FaceColor',[1 0 0]*0.75 + [0 0 1]*0.25 },...
         {'I',[5001 100000000],'FaceColor',[1 0 0]}...
         );
     
% 
%%
% figure(1)
% clf
% mapshow(S,'LineWidth',0.5,'SymbolSpec',ColoredConstituencies);
% set(gca,'XTickLabels',{},'YTickLabels',{},'Visible','off')


%%

 v = VideoWriter('KenyaEpidemic.avi');
 open(v);
 figure(1)
 
 for t = 1:length(T)
     for i = 1:l
         S(i).I = double(Epidemic(i,2,t) + Epidemic(i,5,t));
     end
     
     fig = figure(1);
     clf
     mapshow(S,'LineWidth',0.5,'SymbolSpec',ColoredConstituencies);
     set(gca,'XTickLabels',{},'YTickLabels',{},'Box','on')
     title(strcat('t = ',num2str(T(t)),'days,  (' , num2str(T(t)./365.25), ' years  )' ),'FontSize', 18);
     Frame = getframe(fig);
     writeVideo(v,Frame);
 end


close(v);
