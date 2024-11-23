clc;
clear;
close all;

%%

fileID1 = fopen("./3b - R - NIR BW = 20nm.txt");
fileID2 = fopen("./3b - T - NIR BW = 20nm.txt");

data_R = [];
data_T = [];

for i=1:18
    fgetl(fileID1);
    fgetl(fileID2);
end

while ~feof(fileID1)
    
    new_R = fgetl(fileID1);
    new_T = fgetl(fileID2);

    
    data_R = [data_R; str2num(new_R)];
    data_T = [data_T; str2num(new_T)];

end

figure
plot( data_R(: , 1),  data_R(: , 2) , data_T(: , 1),  data_T(: , 2) ...
    , 'LineWidth',1.5, 'LineStyle', '-' )

title("UV-VIS measurement")

xlabel("Wavelenght")
ylabel("something")
legend("Reflectance", "Transmittance")

xlim([400 1200])

hold on

% plot(data_R(: , 1),  data_R(: , 2), 'LineWidth',1.5, 'Color', 'blue', 'LineStyle', '-')


