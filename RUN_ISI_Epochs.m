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
NS=30;

fprintf('\nStart Processing\n\n')
rstrms=[num2str(r*1000) 'ms'];
fprintf(['\nApEn parameters: m=' num2str(m) ' r=' rstrms '\n\n\n'])

cd(pathResults)
newFolder=['ApEn_Kijon_ISI_phaseran_NS' num2str(NS) '_m' num2str(m) '_r' strrep(rstrms,'.','') '_MFR' num2str(MinFiringRate) '_L' num2str(Lmax)];
mkdir(newFolder)

for jj=1:3
    
    %Variables initialization
    cd (eval(type{jj}))
    list=dir('*.mat');
    if strcmp(type{jj},'supercritico')
        fs=2.5*10^4;
        Nsamp=total_in_minutes*60*fs;
    else
        fs=10^4;
        Nsamp=total_in_minutes*60*fs;
    end
    inizio=1;
    fine=Nsamp;
    actualtype=type{jj};
    MeanFiringRate=zeros(1,length(list));
    EpocsLength=zeros(1,length(list));
    NumberOfEpocs=zeros(1,length(list));
    TtestSuccessRate=zeros(1,length(list));
    electrode=cell(1,length(list));
    apen_all=cell(1,length(list));
    PHI_all=cell(1,length(list));
    apenSH_all=cell(1,length(list));
    TtestResult=cell(1,length(list));
    EpochsFiringRate=cell(1,length(list));
    
    
    for ii=1:length(list)
        
        %ISI extraction
        name=list(ii).name;
        idx=strfind(name,'_');
        electrode{ii}=name((idx(end)+1):(idx(end)+2));
        load(name);
        peak_train_full=full(peak_train);
        clear peak_train
        peak_train=peak_train_full(inizio:fine);
        clear peak_train_full
        spike=find(peak_train~=0);
        MeanFiringRate(ii)=length(spike)/(length(peak_train)/fs);
        if MeanFiringRate(ii) <= MinFiringRate
            apen_all{ii}=NaN;
            PHI_all{ii}=NaN;
            apenSH_all{ii}=NaN;
            TtestResult{ii}=NaN;
            EpochsFiringRate{ii}=NaN;
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
        
        %Epochs firing rate
        EpochsFiringRate{ii}=1./mean(epochs);
        
        %Variable initialization for apen analysis
        apen0=zeros(1,N);
        hp=zeros(N,2);
        PHI=zeros(N,2);
        apen0sh=zeros(N,NS);
        
        %Apen analysis and surrogate data testing
        for j=1:N
            signal=epochs(:,j)';
            signal(signal==0)=[];
            if length(signal)>(3/4)*Lmax || N==1 %protection against the last epochs (few samples)
                [apen0(j), PHI(j,:)]=ApEn( m, r, signal);
                %apen0(j)=SampEn( m, r, signal);
                %surrogate data
                surrogates=phaseran(signal,NS);
                %[surrogates,iter]=IAAFT(signal,NS,1000);
                for s=1:NS
                    [apen0sh(j,s), ~]=ApEn( m, r, surrogates(:,s)');
                    %apen0sh(j,s)=SampEn( m, r, surrogates(:,s)');
                end
                %T-test: h=1-->surrogates greater than original :)
                [hp(j,1),hp(j,2)]=ttest(apen0sh(j,:),apen0(j),'Tail','right');
            else
                apen0(j)=nan;
                hp(j,:)=nan;
                PHI(j,:)=nan;
                apen0sh(j,:)=nan;
                EpochsFiringRate{ii}(j)=nan;
            end
        end
        
        %Store results
        apen_all{ii}=apen0;
        apenSH_all{ii}=apen0sh;
        PHI_all{ii}=PHI;
        TtestResult{ii}=hp;
        TtestSuccessRate(ii)=sum(hp(:,1)==1)/N;
        NumberOfEpocs(ii)=N;
        
        %display the results in real time in the workspace
        results=cell(2,8);
        results{1}=upper(actualtype);
        results{2}=['Electrode=' num2str(ii)];
        results{3}='n° samples';
        results{4}=L;
        results{5}='n° epochs';
        results{6}=N;
        results{7}='ApEn or';
        results{8}=nanmean(apen0);
        results{9}='PHI_+';
        results{10}=nanmean(PHI(:,1));
        results{11}='PHI_-';
        results{12}=nanmean(PHI(:,2));
        results{13}='ApEn surr';
        results{14}=nanmean(apen0sh(:));
        results{15}='T-test rate';
        results{16}=TtestSuccessRate(ii);
        disp(results)
    end
    
    % Save results
    cd([pathResults '\' newFolder])
    save([actualtype '_' mfilename '_m' num2str(m) '_r' strrep(rstrms,'.','') '_' num2str(total_in_minutes) 'min.mat'],...
        'apen_all','apenSH_all','PHI_all','TtestResult','EpochsFiringRate',...
        'TtestSuccessRate','r','m','MeanFiringRate','MinFiringRate',...
        'total_in_minutes','fs','actualtype','NumberOfEpocs')
    
end



%% Organization of the results
clear all, close all
%Variables initialization
type=string(["subcritico","critico","supercritico"]);
list=dir('*.mat');
idx=zeros(1,3);

%Order the file list
for i=1:3
    name=list(i).name;
    name=name(1:regexp(name,'_','once')-1);
    for ii=1:3
        if strcmp(name,type(ii))
            idx(i)=ii;
        end
    end
    clear name
end

% Good epochs identification based on the surrogates right tailed T-test
for i=1:3
    load(list(idx(i)).name)
    apen_good=zeros(1,60);
    NumberGoodEpochs=zeros(1,60);
    %Select only the epochs with positive T-test
    for j=1:60
        if ~isempty(TtestResult{1,j})
            h=TtestResult{j}(:,1);
            apen_good(j)=mean(apen_all{j}(h==1));
            NumberGoodEpochs(j)=length(apen_all{j}(h==1));
        else
            apen_good(j)=NaN;
            NumberGoodEpochs(j)=0;
        end
    end
    save(list(idx(i)).name,'apen_good','NumberGoodEpochs','-append')
end

% Plot all the results

%Linear plot
for i=1:3
    load(list(idx(i)).name)
    figure; hold on;
    for ii=1:60
        y=apen_all{ii};
        x=ii+(0:length(y)-1)./10;
        plot(x,y,'ok','MarkerFaceColor','none')
        ysh=apenSH_all{ii};
        plot(x,ysh,'.r')
    end
    title(actualtype)
    ylabel('ApEn (red = surrogates)')
    xlabel('Electrodes')
    grid on; grid minor;
end

%Boxplot
GROUP="";
APEN=[];
for i=1:3
    load(list(idx(i)).name)
    y=cellfun(@nanmean,apen_all);
    y(isnan(y))=[];
    APEN=[APEN y];
    GROUP=[GROUP; string(repmat(actualtype,length(y),1))];
end
GROUP(1)=[];
figure;
boxplot(APEN,GROUP);
ylabel('ApEn [a.u.]');
%title(['ALL n° electrodes sub=' num2str(sum(strcmp(GROUP,"subcritico"))) ' cr=' num2str(sum(strcmp(GROUP,"critico"))) ' sup=' num2str(sum(strcmp(GROUP,"supercritico")))])
saveppt('Figures.ppt','','-f')

%Map
lim=round(get(gca,'YLim'),3);
for i=1:3
    load(list(idx(i)).name)
    y=cellfun(@nanmean,apen_all); 
    mApS(y,'','',lim,1);
    %saveppt('Figures.ppt','','-f')
end
mApS(y,'Colorbar','',lim,1)
%saveppt('Figures.ppt','','-f')

% Plot only good results

%Boxplot
GROUP="";
APEN=[];
for i=1:3
    load(list(idx(i)).name)
    y=apen_good(~isnan(apen_good));
    APEN=[APEN y];
    GROUP=[GROUP; string(repmat(actualtype,length(y),1))];
end
GROUP(1)=[];
figure; 
boxplot(APEN,GROUP);
ylabel=('ApEn [a.u.]');
%title(['GOOD n° electrodes sub=' num2str(sum(strcmp(GROUP,"subcritico"))) ' cr=' num2str(sum(strcmp(GROUP,"critico"))) ' sup=' num2str(sum(strcmp(GROUP,"supercritico")))])
saveppt('Figures.ppt','','-f')

% Map
lim=round(get(gca,'YLim'),3);
for i=1:3
    load(list(idx(i)).name)
    mApS(apen_good,'','',lim,1)
    saveppt('Figures.ppt','','-f')
end
mApS(apen_good,'Colorbar','',lim,1)
saveppt('Figures.ppt','','-f')