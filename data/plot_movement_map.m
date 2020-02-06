%% Spatial map for movements
load('mixing_matrix.mat');
T = readtable('combined_population_estimates.csv');
lats = T.y_location/110.57;
longs = T.x_location/111.32;
%%


S = shaperead('County.shp');
l = length(S);
% Sort so in the same order as data
[~,index] = sortrows({S.COUNTY}.'); S = S(index); clear index

%%

M = movement_mat;
for i = 1:47
    M(i,i) = 0;
end

M = M(:);
% figure(2)
% clf
% histogram(M)
% set(gca,'YScale','log','XScale','log')


M_less = M(M > 0.001);
max_M = max(M_less);
min_M = min(M_less);
%%

figure(1)
clf
mapshow(S,'LineWidth',0.5);
set(gca,'XTickLabels',{},'YTickLabels',{},'Visible','off')
hold on
scatter(longs,lats)
for i = 1:47
    for j = 1:47
        if movement_mat(i,j) > 0.001 & i ~= j
            scale = (log(movement_mat(i,j)) - log(min_M) )/(log(max_M) - log(min_M));
            c = (1-scale)*[0 0 1 0.5] + scale*[1 0 0 0.5];
        quiver(longs(j),lats(j),longs(i)-longs(j),lats(i)-lats(j),'LineWidth',3,'MaxHeadSize',0.5,'Color',c);
        end
    end
end
title('Mixing intensity')


%%
