%% Spike train analysis with an entropy index
home, close all, clear all

%Paths
root=fileparts(which(mfilename));
pathResults=[root '\Results'];
critico=[root '\01_CRITICO\exp360_2h_10kHz_25div_1\360_01_PeakDetectionMAT_files_2msec\ptrain_360_01_DIV25_Nbasl'];
subcritico=[root '\02_SUBCRITICO\exp360_2h_10kHz_25div_2\360_02_PeakDetectionMAT_files_2ms\ptrain_360_02_DIV25_Nbasl'];
supercritico=[root '\03_SUPERCRITICO\exp376_1h_25kHz\376_PeakDetectionMAT_files_2msec\ptrain_376_01_DIV23_Nbasal'];
cd(root)

%Parameters
type=string(["subcritico","critico","supercritico"]);
MinFiringRate=1; %spikes/sec
Lmax=2500;
overlap_percentage=0.5;
m=3;
r=0.001;
total_in_minutes=60;
partial_in_minutes=10;

fprintf('\nStart Processing\n\n')
rstrms=[num2str(r*1000) 'ms'];
fprintf(['\nApEn parameters: m=' num2str(m) ' r=' rstrms '\n\n\n'])

newFolder=['SampEn_Kijon_ISI_subepochsFixed_' num2str(partial_in_minutes) 'min_m' num2str(m) '_r' strrep(rstrms,'.','') '_MFR' num2str(MinFiringRate) '_L' num2str(Lmax)];
cd(pathResults)
mkdir(newFolder)

for jj=1:3
    for ss=1:floor(total_in_minutes/partial_in_minutes)
        
        %Sampling frequency initialization
        cd (eval(type{jj}))
        list=dir('*.mat');
        if strcmp(type{jj},'supercritico')
            fs=2.5*10^4;
            Nsamp=partial_in_minutes*60*fs;
        else
            fs=10^4;
            Nsamp=partial_in_minutes*60*fs;
        end
        actualtype=type{jj};
        inizio=1+(ss-1)*Nsamp;
        fine=ss*Nsamp;
        
        %variables initialization
        MeanFiringRate=zeros(1,length(list));
        EpocsLength=zeros(1,length(list));
        NumberOfEpocs=zeros(1,length(list));
        electrode=cell(1,length(list));
        apen_all=cell(1,length(list));
        PHI_all=cell(1,length(list));
        
        for ii=1:length(list)
            
            %ISI extraction
            names=list(ii).name;
            idx=strfind(names,'_');
            electrode{ii}=names((idx(end)+1):(idx(end)+2));
            load(names);
            peak_train_full=full(peak_train);
            clear peak_train            
            peak_train=peak_train_full(inizio:fine);
            clear peak_train_full
            spike=find(peak_train~=0);
            MeanFiringRate(ii)=length(spike)/(length(peak_train)/fs);
            if MeanFiringRate(ii) <= MinFiringRate
                apen_all{ii}=nan;
                PHI_all{ii}=nan;
                continue
            end
            clear peak_train
            ISI=(diff(spike)./fs)';
            L=length(ISI);
            
            %ISI segmentation into epochs
            if L>(1+overlap_percentage)*Lmax
                epochs=buffer(ISI,Lmax,floor(Lmax*overlap_percentage),'nodelay');
                N=size(epochs,2);
                EpocsLength(ii)=Lmax;
            else
                epochs=ISI';
                N=1;
                EpocsLength(ii)=L;
            end
            
            %Variable initialization for apen analysis
            apen0=zeros(1,N);
            PHI=zeros(N,2);
            
            %Apen analysis
            for j=1:N
                signal=epochs(:,j)';
                signal(signal==0)=[];
                if length(signal)>(3/4)*Lmax || N==1 %protection against the last epochs (few samples)
                    %[apen0(j), PHI(j,:)]=ApEn( m, r, signal);
                    apen0(j)=SampEn( m, r, signal);
                    if apen0(j)==inf %it's used only with SampEn
                        apen0(j)=nan;
                    end
                else
                    apen0(j)=nan;
                    PHI(j,:)=nan;
                end
            end
            
            %Store results
            apen_all{ii}=apen0;
            PHI_all{ii}=PHI;
            NumberOfEpocs(ii)=N;
            
            %display the results in real time in the workspace
            results=cell(2,4);
            results{1}=upper(actualtype);
            results{2}=['Electrode=' electrode{ii}];
            results{3}='n° samples';
            results{4}=L;
            results{5}='n° epochs';
            results{6}=N;
            results{7}='ApEn';
            results{8}=nanmean(apen0);
            disp(results)
        end
        
        % Save results
        cd([pathResults '\' newFolder])
        save([actualtype '_' mfilename '_m' num2str(m) '_r' strrep(rstrms,'.','') '_' num2str(partial_in_minutes*ss) 'min_#' num2str(ss) '.mat'],...
            'apen_all','PHI_all','partial_in_minutes','ss',...
            'r','m','MeanFiringRate','MinFiringRate',...
            'total_in_minutes','fs','actualtype','NumberOfEpocs')
    end
end
clear all


%% Organization of the results
clear all, close all
%Variables initialization
type=string(["subcritico","critico","supercritico"]);
list=dir('*.mat');
L=length(list);
rho=zeros(1,L);
pvalue=zeros(1,L);
names=strings(1,L);
apen_SUB=cell(1,L);
apen_TOT=cell(1,L);

%Order the file list
dates={list.date}';
[~,I]=sort(dates);
listordered=list(I);
load(listordered(1).name,'partial_in_minutes')


%Correlation Boxplot
for i=1:L
    apen_SUB{i}=load(listordered(i).name,'apen_all');
    apen_SUB{i}=cellfun(@nanmean,apen_SUB{i}.apen_all)';
    name=listordered(i).name;
    names{i}=name(1:regexp(name,'_','once')-1);
    actualdirectory=pwd;
    cd('..\..\results\SampEn_Kijon_ISI_iaaft_NS30_m3_r1ms_MFR1_L2500')
    list=dir;
    nameT={list.name}';
    for j=1:length(nameT)
        idx=regexp(nameT{j}, '_', 'once');
        if ~isempty(idx) && strcmp(names{i},nameT{j}(1:idx-1))==1
            break
        end
    end
    apen_TOT{i}=load(nameT{j},'apen_all');
    apen_TOT{i}=cellfun(@nanmean,apen_TOT{i}.apen_all)';
    [rho(i),pvalue(i)]=corr(apen_SUB{i},apen_TOT{i},'Rows','pairwise'); %try also pairwise
    cd(actualdirectory)
end

%partial-total
figure;
boxplot(rho,names);
ylabel('\rho [a.u.]');
ylim([0 1])
title(['Pearson correlation coefficent - ' num2str(partial_in_minutes) 'min vs 60min'])
saveppt('Figures.ppt','','-f')

%partial-partial
for i=1:3
    RHOtemp=corr(cell2mat(apen_SUB(1+(i-1)*(L/3):i*(L/3))),'Rows','complete');
    RHOtemp=triu(RHOtemp);
    RHOtemp(1:1+size(RHOtemp,1):end)=0;
    RHOtemp=RHOtemp(:);
    RHOtemp(RHOtemp==0)=[];
    RHO(:,i)=RHOtemp;
    clear RHOtemp
end
figure;
boxplot(RHO,type);
ylim([0 1])
ylabel('\rho [a.u.]');
title(['Pearson correlation coefficent - ' num2str(partial_in_minutes) 'min vs ' num2str(partial_in_minutes) 'min all combo'])
saveppt('Figures.ppt','','-f')

%Map
for i=1:L/3:L
    for ii=[1 3 6]
        load(listordered(ii-1+i).name)
        disp(listordered(ii-1+i).name)
        y=cellfun(@nanmean,apen_all);
        mApS(y,'',[],[0.8 2.5],1);
        saveppt('Figures.ppt','','-f')
    end
end
mApS(y,'SampEn',[],[0.8 2.5],1);
saveppt('Figures.ppt','','-f')