%Conversion script from risk regions to counties

T= readtable('2019_population_estimates_withriskregions.xls');
conversion_matrix = zeros(47,20);
N = T.Total_Population19(1:47);
N_rr = T.Total_population_2019(1:20);
%%
conversion_matrix(1,3) = (N(1)/N_rr(3));%Baringo

conversion_matrix(2,13) = (N(2)/N_rr(13));%Bomet

conversion_matrix(3,15) = (N(3)/N_rr(15));%Bungoma

conversion_matrix(4,15) = (N(4)/N_rr(15));%Busia

conversion_matrix(5,5) = (N(5)/N_rr(5));%Elgeyo Marakwet

conversion_matrix(6,1) = (N(6)/N_rr(1));%Embu

conversion_matrix(7,7) = (N(7)/N_rr(7));%Garissa

conversion_matrix(8,18) = (N(8)/N_rr(18));%Homa Bay

conversion_matrix(9,11) = (N(9)/N_rr(11));%Isiolo

conversion_matrix(10,4) = (0.5*N(10)/N_rr(4));%Kajiado
conversion_matrix(10,17) = (0.5*N(10)/N_rr(17));%Kajiado

conversion_matrix(11,15) = (N(11)/N_rr(15));%Kakamega
conversion_matrix(12,13) = (N(12)/N_rr(13));%Kericho
conversion_matrix(13,4) = (N(13)/N_rr(4));%Kiambu

conversion_matrix(14,12) = (N(14)/N_rr(12))/3;%Kilifi
conversion_matrix(14,14) = (N(14)/N_rr(14))/3;%Kilifi
conversion_matrix(14,20) = (N(14)/N_rr(20))/3;%Kilifi

conversion_matrix(15,1) = (N(15)/N_rr(1));%Kirinyaga

conversion_matrix(16,18) = (N(16)/N_rr(18));%Kisii

conversion_matrix(17,13) = (0.5*N(17)/N_rr(13));%Kisumu
conversion_matrix(17,15) = (0.5*N(17)/N_rr(15));%Kisumu

conversion_matrix(18,11) = (0.5*N(18)/N_rr(11));%Kitui
conversion_matrix(18,17) = (0.5*N(18)/N_rr(17));%Kitui

conversion_matrix(19,12) = (0.5*N(19)/N_rr(12));%Kwale
conversion_matrix(19,20) = (0.5*N(19)/N_rr(20));%Kwale

conversion_matrix(20,3) = (0.5*N(20)/N_rr(3));%Laikipia
conversion_matrix(20,19) = (0.5*N(20)/N_rr(19));%Laikipia

conversion_matrix(21,14) = (N(21)/N_rr(14));%Lamu

conversion_matrix(22,1) = (0.5*N(22)/N_rr(1));%Machakos
conversion_matrix(22,4) = (0.5*N(22)/N_rr(4));%Machakos

conversion_matrix(23,17) = (N(23)/N_rr(17));%Makueni

conversion_matrix(24,2) = (N(24)/N_rr(2));%Mandera

conversion_matrix(25,6) = (N(25)/N_rr(6));%Marsabit

conversion_matrix(26,11) = (N(26)/N_rr(11));%Meru

conversion_matrix(27,18) = (N(27)/N_rr(18));%Migori

conversion_matrix(28,12) = (N(28)/N_rr(12));%Mombasa

conversion_matrix(29,1) = (N(29)/N_rr(1));%Murang'a

conversion_matrix(30,4) = (N(30)/N_rr(4));%Nairobi

conversion_matrix(31,3) = (0.5*N(31)/N_rr(3));%Nakuru
conversion_matrix(31,8) = (0.5*N(31)/N_rr(8));%Nakuru

conversion_matrix(32,15) = (N(32)/N_rr(15));%Nandi

conversion_matrix(33,8) = (0.5*N(33)/N_rr(8));%Narok
conversion_matrix(33,18) = (0.5*N(33)/N_rr(18));%Narok

conversion_matrix(34,13) = (N(34)/N_rr(13));%Nyamira

conversion_matrix(35,3) = (0.5*N(35)/N_rr(3));%Nyandarua
conversion_matrix(35,8) = (0.5*N(35)/N_rr(8));%Nyandarua

conversion_matrix(36,1) = (N(36)/N_rr(1));%Nyeri

conversion_matrix(37,19) = (N(37)/N_rr(19));%Samburu

conversion_matrix(38,15) = (N(38)/N_rr(15));%Siaya

conversion_matrix(39,10) = (N(39)/N_rr(10));%Taita Taveta

conversion_matrix(40,7) = (0.5*N(40)/N_rr(7));%Tana River
conversion_matrix(40,14) = (0.5*N(40)/N_rr(14));%Tana River

conversion_matrix(41,1) = (N(41)/N_rr(1));%Tharaka Nithi

conversion_matrix(42,5) = (N(42)/N_rr(5));%Trans Nzoia

conversion_matrix(43,9) = (N(43)/N_rr(9));%Turkana

conversion_matrix(44,5) = (N(44)/N_rr(5));%Uasin Gishu

conversion_matrix(45,15) = (N(45)/N_rr(15));%Vihiga

conversion_matrix(46,16) = (N(46)/N_rr(16));%Wajir

conversion_matrix(47,5) = (N(47)/N_rr(5));%West Pokot

%%

conversion_matrix_c_to_rr = zeros(20,47);


conversion_matrix_c_to_rr(3,1) = 1 ;%Baringo

conversion_matrix_c_to_rr(13,2) = 1 ;%Bomet

conversion_matrix_c_to_rr(15,3) = 1;%Bungoma

conversion_matrix_c_to_rr(15,4) = 1;%Busia

conversion_matrix_c_to_rr(5,5) = 1;%Elgeyo Marakwet

conversion_matrix_c_to_rr(1,6) = 1;%Embu

conversion_matrix_c_to_rr(7,7) = 1;%Garissa

conversion_matrix_c_to_rr(18,8) =1;%Homa Bay

conversion_matrix_c_to_rr(11,9) = 1;%Isiolo

conversion_matrix_c_to_rr(4,10) = 0.5;%Kajiado
conversion_matrix_c_to_rr(17,10) = 0.5;%Kajiado

conversion_matrix_c_to_rr(15,11) =1;%Kakamega

conversion_matrix_c_to_rr(13,12) = 1;%Kericho

conversion_matrix_c_to_rr(4,13) = 1;%Kiambu

conversion_matrix_c_to_rr(12,14) = 1/3;%Kilifi
conversion_matrix_c_to_rr(14,14) = 1/3;%Kilifi
conversion_matrix_c_to_rr(20,14) = 1/3;%Kilifi

conversion_matrix_c_to_rr(1,15) = 1;%Kirinyaga

conversion_matrix_c_to_rr(18,16) = 1;%Kisii

conversion_matrix_c_to_rr(13,17) = 0.5;%Kisumu
conversion_matrix_c_to_rr(15,17) = 0.5;%Kisumu

conversion_matrix_c_to_rr(11,18) = 0.5;%Kitui
conversion_matrix_c_to_rr(17,18) = 0.5;%Kitui

conversion_matrix_c_to_rr(12,19) = 0.5;%Kwale
conversion_matrix_c_to_rr(20,19) = 0.5;%Kwale

conversion_matrix_c_to_rr(3,20) = 0.5;%Laikipia
conversion_matrix_c_to_rr(19,20) = 0.5;%Laikipia

conversion_matrix_c_to_rr(14,21) = 1;%Lamu

conversion_matrix_c_to_rr(1,22) = 0.5;%Machakos
conversion_matrix_c_to_rr(4,22) = 0.5;%Machakos

conversion_matrix_c_to_rr(17,23) = 1;%Makueni

conversion_matrix_c_to_rr(2,24) = 1;%Mandera

conversion_matrix_c_to_rr(6,25) = 1;%Marsabit

conversion_matrix_c_to_rr(11,26) = 1;%Meru

conversion_matrix_c_to_rr(18,27) = 1;%Migori

conversion_matrix_c_to_rr(12,28) = 1;%Mombasa

conversion_matrix_c_to_rr(1,29) = 1;%Murang'a

conversion_matrix_c_to_rr(4,30) = 1;%Nairobi

conversion_matrix_c_to_rr(3,31) = 0.5;%Nakuru
conversion_matrix_c_to_rr(8,31) = 0.5;%Nakuru

conversion_matrix_c_to_rr(15,32) = 1;%Nandi

conversion_matrix_c_to_rr(8,33) = 0.5;%Narok
conversion_matrix_c_to_rr(18,33) = 0.5;%Narok

conversion_matrix_c_to_rr(13,34) = 1;%Nyamira

conversion_matrix_c_to_rr(3,35) = 0.5;%Nyandarua
conversion_matrix_c_to_rr(8,35) = 0.5;%Nyandarua

conversion_matrix_c_to_rr(1,36) = 1;%Nyeri

conversion_matrix_c_to_rr(19,37) = 1;%Samburu

conversion_matrix_c_to_rr(15,38) = 1;%Siaya

conversion_matrix_c_to_rr(10,39) = 1;%Taita Taveta

conversion_matrix_c_to_rr(7,40) = 0.5;%Tana River
conversion_matrix_c_to_rr(14,40) =0.5;%Tana River

conversion_matrix_c_to_rr(1,41) = 1;%Tharaka Nithi

conversion_matrix_c_to_rr(5,42) = 1;%Trans Nzoia

conversion_matrix_c_to_rr(9,43) = 1;%Turkana

conversion_matrix_c_to_rr(5,44) = 1;%Uasin Gishu

conversion_matrix_c_to_rr(15,45) = 1;%Vihiga

conversion_matrix_c_to_rr(16,46) = 1;%Wajir

conversion_matrix_c_to_rr(5,47) = 1;%West Pokot

