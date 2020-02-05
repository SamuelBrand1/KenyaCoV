% Script for plotting counties

%Script for plotting Kenyan district data

% S = shaperead('ke_district_boundaries.shp');
 S = shaperead('County.shp');
l = length(S);
%%
% Create random color on grey to red scale
% r = rand(l,1);
% 
% for i = 1:l
% %    S(i).RandColor = (1-r(i))*[0.7 0.7 0.7] + r(i)*[1 0 0]; 
%     name = S(i).CONSTITUEN;
%     F = find(strcmp(KenyaTbl.Constituency,name));
%      S(i).TotalPop = KenyaTbl.Total_ConstitTable(F);
%      S(i).Density = KenyaTbl.DensityPeoplePerSq_Km(F);
% 
% end
%%
%Make symbol specification for the polygon geometry

% ColoredCounties = makesymbolspec('Polygon',{'Default','FaceColor',S(1).RandColor});
% ColoredCounties = makesymbolspec('Polygon',{'RandColor',[0 0.3333],'FaceColor',[0.7 0.7 0.7]},...
%         {'RandColor',[0.3333 0.6667],'FaceColor',0.5*[0.7 0.7 0.7] + 0.5*[1 0 0]},{'RandColor',[0.6667 1],'FaceColor',[1 0 0]});
% 
% 
ColoredConstituencies = makesymbolspec('Polygon',{'Density',[0 500],'FaceColor',[0.7 0.7 0.7]},...
        {'Density',[500 15000],'FaceColor',[1 0.5 0]},{'Density',[15000 26000],'FaceColor',[1 0 0]});
% 
%%
figure(1)
clf
mapshow(S,'LineWidth',0.5);
%%

set(gca,'XTickLabels',{},'YTickLabels',{},'Visible','off')
% hold on
% mapshow(S2)
