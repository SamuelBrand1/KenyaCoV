%In scenario 1 with maximal contact capacity of 1e3:       25% detection: 33201.0,     50%: 70726.5,       75%:113773.0,        90%: 142053.5
%In scenario 2 with maximal contact capacity of 5e3:       25% detection: 42891.0,     50%: 92065.5,       75%:149219.5,        90%: 188101.5
%In scenario 3 with maximal contact capacity of 1e4:       25% detection: 85281.0,     50%: 124767.5,    75%: 169125.5,        90%: 198568.5

%data=[ MAX_CAP1e3[0%detection+50tracedPerDetected 0%100] [25%50 25%100]
%... ; MAX_CAP5e3... ]
data=[[0 0] [33201 34201] [70726 75000] [113773 100000] [142053 150000] ; [0 0] [42891 45000] [92065 80000] [149219 160000] [188101 180000]  ; [0 0] [85281 90000] [124767 130000] [169125 200000] [198568 220000]  ];

data2_1=[0 33201 70726  113773 142053 ; 0 42891 92065 149219 188101  ; 0 85281 124767 169125 198568 ];
data2_2=[ 0 34201 75000 100000 150000 ; 0 45000 80000 160000 180000  ; 0 90000 130000 200000 220000 ];


riskregionnames = ["Machakos/Muranga" "Mandera" "Baringo/Nakuru" "Nairobi" "Uasin Gishu" "Marsabit" "Garissa" "Nakuru/Narok" "Turkana" "Taita Taveta" "Kitui/Meru" "Kilifi/Mombasa" "Kericho/Kisumu" "Kilifi/Lamu" "Kakamega/Kisumu" "Wajir" "Kajiado/Kisumu" "Homa bay/Migori" "Samburu/Laikipia" "Kilifi/Kwale" "Total"];
wa_coords=[300 450;515 85; 165 360; 235 465; 115 300; 300 140; 510 380; 180 430; 120 130; 355 615; 340 375; 465 630; 100 380; 495 530; 40 355; 490 255; 155 495; 30 440; 250 290; 400 670];
wa_max_x=max(wa_coords(:,1));wa_max_y=max(wa_coords(:,2));
data_wa=randi(1000,20,1);

data_map=[wa_coords data_wa];

data_matrix=zeros(10,10);
% for i=1:20
%     disp(floor(data_map(i,1)*10/wa_max_x)+"/"+floor(data_map(i,2)*10/wa_max_x));
% end
for i=1:20
    if data_matrix(ceil(data_map(i,1)*10/wa_max_x),ceil(data_map(i,2)*10/wa_max_y))==0
        %disp(i)
    end
%     x=floor(data_map(i,1)*10/wa_max_x);
%     y=floor(data_map(i,2)*10/wa_max_y);
%     if x==0 
%         x=1;
%     end
%     if y==0
%         y=1;
%    data_matrix(x,y)=data_map(i,3);
    data_matrix(ceil(data_map(i,1)*10/wa_max_x),ceil(data_map(i,2)*10/wa_max_y))=data_map(i,3);
end
h=heatmap(data_matrix);
h.YDisplayData = flipud(h.YDisplayData);

wa=[transpose(riskregionnames(1:20)) data_wa];