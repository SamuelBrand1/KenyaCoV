%% Spatial map for movements
load('/Users/sam/Github/KenyaCoV/output/peaktimes_by_county.mat');
S = shaperead('County.shp');
l = length(S);
% Sort so in the same order as data
[~,index] = sortrows({S.COUNTY}.'); S = S(index); clear index
%% Add the time to peak 
for i = 1:47
      S(i).timetopeak = median(peaktimes_by_county(:,i));
end

%% Make symbol specification for the polygon geometry
    x = linspace(0,1,7);
     ColoredConstituencies = makesymbolspec('Polygon',...
         {'timetopeak',[55 60],'FaceColor',[1 0 0]*(1-x(1)) + [0 0 1]*x(1) },...
         {'timetopeak',[60 65],'FaceColor',[1 0 0]*(1-x(2)) + [0 0 1]*x(2) },...
         {'timetopeak',[65 70],'FaceColor',[1 0 0]*(1-x(3))+ [0 0 1]*x(3) },...      
         {'timetopeak',[70 75],'FaceColor',[1 0 0]*(1-x(4)) + [0 0 1]*x(4) },...              
         {'timetopeak',[75 80],'FaceColor',[1 0 0]*(1-x(5)) + [0 0 1]*x(5) },...
         {'timetopeak',[80 85],'FaceColor',[1 0 0]*(1-x(6)) + [0 0 1]*x(6) },...
         {'timetopeak',[85 90],'FaceColor',[1 0 0]*(1-x(7)) + [0 0 1]*x(7) });

     %%
     simplecmap = [];
     for i = 1:7
         simplecmap = [simplecmap;[1 0 0]*(1-x(i)) + [0 0 1]*x(i)];
     end
     

%%
figure(1)
clf
colormap(simplecmap);
h_kenya = axes;
h_kenya.Position =  [0.100 0.100 0.55 0.85];
h_kenya.FontSize = 18;
h_Nairobi = axes;
h_Nairobi.Position =  [0.7 0.600 0.25 0.35];
h_Nairobi.FontSize = 18;
h_mombasa = axes;
h_mombasa.Position =  [0.7 0.100 0.25 0.35];
h_mombasa.FontSize = 18;
axes(h_kenya);
mapshow(S,'LineWidth',0.5,'SymbolSpec',ColoredConstituencies);
title('Kenya - median time to epidemic peak (Nairobi seed)','FontSize',20)
cb = colorbar;
cb.Ticks = x;
cb.TickLabels = [55;60;65;70;75;80;85];
% 
axes(h_Nairobi);
mapshow(S(30),'SymbolSpec',ColoredConstituencies);
title('Nairobi')
% 
axes(h_mombasa);
mapshow(S(28),'SymbolSpec',ColoredConstituencies);
title('Mombasa')




%%
