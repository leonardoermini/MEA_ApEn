% FUNZIONE CHE GENERA IL RASTER PLOT DEL NETWORK 

%Paths
root=fileparts(which(mfilename));
pathResults=[root '\Results'];
critico=[root '\01_CRITICO\exp360_2h_10kHz_25div_1\360_01_PeakDetectionMAT_files_2msec\ptrain_360_01_DIV25_Nbasl'];
subcritico=[root '\02_SUBCRITICO\exp360_2h_10kHz_25div_2\360_02_PeakDetectionMAT_files_2ms\ptrain_360_02_DIV25_Nbasl'];
supercritico=[root '\03_SUPERCRITICO\exp376_1h_25kHz\376_PeakDetectionMAT_files_2msec\ptrain_376_01_DIV23_Nbasal'];

%%%Choose the folder%%%%%
cd(critico); Fs=1e4;
%%%%%%%%%%%%%%%%%%%%%%%%%

list=dir('*.mat');
Data=zeros(length(list),60*Fs);
for ii=1:length(list)
    filename=list(ii).name;
    m=matfile(filename);
    Data(ii,:)=m.peak_train(1:60*Fs,1)';
end

indeces=find(Data);
L_peak_train=size(Data,2);
[X,Y]=meshgrid(linspace(0,(L_peak_train-1)/Fs,L_peak_train),1:size(Data,1));
figure
plot(X(indeces),Y(indeces),'.');
xlabel('time, \it s')
ylabel('electrode #')
grid minor
ylim([0 size(Data,1)+1])
xlim([0 size(Data,2)/Fs])


