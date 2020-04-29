clear;
load selected_incidence.mat;
load conversion_matrix.mat
S = shaperead('County.shp');
l = length(S);
% Sort so in the same order as data
[~,index] = sortrows({S.COUNTY}.'); S = S(index); clear index;

%% convert the daily incidence from risk region into county --- choose which scenario 

DI = conversion_matrix*inc_025_CI;
% DI = conversion_matrix*inc_025_nc;
% DI = conversion_matrix*inc_0_nc;


%% Add the median to daily incidence --- choose which day
max_inc = 10001;
for i = 1:47
      S(i).ID = i;
end

%% Make symbol specification for the polygon geometry
map = [0.7 0.7 0.7;jet(99)];



times = 7:7:70;
%%

for t = times
    I = floor(100*(log10(DI(1,t)+1)/log10(max_inc+1))) + 1;
    % I = round(100*((DI(1,t))/(max_inc))) + 1;

    ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',map(I,:)} );
    for i = 2:47
        I = floor(100*(log10(DI(i,t)+1)/log10(max_inc))) + 1;
    %     I = round(100*((DI(i,t))/(max_inc))) + 1;
        ColoredConstituencies.FaceColor{i,1} = 'ID';
        ColoredConstituencies.FaceColor{i,2} = i;
        ColoredConstituencies.FaceColor{i,3} = map(I,:);
    end


    figure(1)
    clf
    % colormap(simplecmap);
    colormap(map);
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
    title(['Daily incidence on day' ' ' num2str(t) ' ' '(median)'],'FontSize',20)
    cb = colorbar;
    cb.Ticks = [log10(1)/log10(max_inc);log10(11)/log10(max_inc);log10(101)/log10(max_inc);log10(1001)/log10(max_inc);log10(10001)/log10(max_inc)];
    cb.TickLabels = [0;10;100;1000;10000];
    cb.Limits = [0 1];
    % 
    axes(h_Nairobi);
    mapshow(S(30),'SymbolSpec',ColoredConstituencies);
    title('Nairobi')
    % 
    axes(h_mombasa);
    mapshow(S(28),'SymbolSpec',ColoredConstituencies);
    title('Mombasa')
    week_num = round(t/7);
    saveas(gcf,['week_' num2str(week_num) '.png'])
end


%% Make a movie

times = 1:150;
v = VideoWriter('daily_incidence.avi');
open(v);
for t = times
    I = floor(100*(log10(DI(1,t)+1)/log10(max_inc+1))) + 1;
    % I = round(100*((DI(1,t))/(max_inc))) + 1;

    ColoredConstituencies = makesymbolspec('Polygon',{ 'ID',1,'FaceColor',map(I,:)} );
    for i = 2:47
        I = floor(100*(log10(DI(i,t)+1)/log10(max_inc))) + 1;
    %     I = round(100*((DI(i,t))/(max_inc))) + 1;
        ColoredConstituencies.FaceColor{i,1} = 'ID';
        ColoredConstituencies.FaceColor{i,2} = i;
        ColoredConstituencies.FaceColor{i,3} = map(I,:);
    end


    figure(1)
    clf
    % colormap(simplecmap);
    colormap(map);
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
    title(['Daily incidence on day' ' ' num2str(t) ' ' '(median)'],'FontSize',20)
    cb = colorbar;
    cb.Ticks = [log10(1)/log10(max_inc);log10(11)/log10(max_inc);log10(101)/log10(max_inc);log10(1001)/log10(max_inc);log10(10001)/log10(max_inc)];
    cb.TickLabels = [0;10;100;1000;10000];
    cb.Limits = [0 1];
    % 
    axes(h_Nairobi);
    mapshow(S(30),'SymbolSpec',ColoredConstituencies);
    title('Nairobi')
    % 
    axes(h_mombasa);
    mapshow(S(28),'SymbolSpec',ColoredConstituencies);
    title('Mombasa')
%     week_num = round(t/7);
%     saveas(gcf,['week_' num2str(week_num) '.png'])
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);


%%



