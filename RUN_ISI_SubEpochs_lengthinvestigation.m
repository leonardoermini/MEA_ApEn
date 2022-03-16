%% Spike train analysis with an entropy index
home, close all, clear all

%Paths
root=fileparts(which(mfilename));
pathResults=[root '\results'];
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

fprintf('\nStart Processing\n\n')
rstrms=[num2str(r*1000) 'ms'];
fprintf(['\nApEn parameters: m=' num2str(m) ' r=' rstrms '\n\n\n'])

newFolder=['SampEn_Kijon_ISI_SubepochsInvestigation_m' num2str(m) '_r' strrep(rstrms,'.','') '_MFR' num2str(MinFiringRate) '_L' num2str(Lmax)];
cd(pathResults)
mkdir(newFolder)

for jj=3:3
    
    for tt=1:total_in_minutes
        
        %culture differences: sampling freq
        cd (eval(type{jj}))
        list=dir('*.mat');
        if strcmp(type{jj},'supercritico')
            fs=2.5*10^4;
        else
            fs=10^4;
        end
        actualtype=type{jj};
        Nsamp=tt*60*fs;
        inizio=1;
        fine=Nsamp;
        
        %variables initialization
        MeanFiringRate=zeros(1,length(list));
        EpocsLength=zeros(1,length(list));
        NumberOfEpocs=zeros(1,length(list));
        electrode=cell(1,length(list));
        apen_all=cell(1,length(list));
        PHI_all=cell(1,length(list));
        
        for ii=1:length(list)
            
            %sparse matrix reconstruction
            names=list(ii).name;
            idx=strfind(names,'_');
            electrode{ii}=names((idx(end)+1):(idx(end)+2));
            load(names);
            peak_train_full=full(peak_train);
            clear peak_train
            
            %period of interest
            
            
            %ISI extraction
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
            
            %variable initialization for apen analysis
            apen0=zeros(1,N);
            PHI=zeros(N,2);
            
            %apen analysis
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
            
            %results storage
            apen_all{ii}=apen0;
            PHI_all{ii}=PHI;
            NumberOfEpocs(ii)=N;
            
            %real time results visualization
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
        
        %results saving
        cd([pathResults '\' newFolder])
        save([actualtype '_' num2str(tt) 'min_' mfilename '_m' num2str(m) '_r' strrep(rstrms,'.','') '.mat'],...
            'apen_all','PHI_all','tt',...
            'r','m','MeanFiringRate','MinFiringRate',...
            'total_in_minutes','fs','actualtype','NumberOfEpocs')
    end
end


%% Result files initialization
clear all, close all
%Variables initialization
type=string(["subcritico","critico","supercritico"]);
list=dir('*.mat');
L=length(list);

%Order the file list
dates={list.date}';
[~,I]=sort(dates);
listordered=list(I);

% Correlation among partial and total
rho=zeros(L/3,3);
apen_SUB=cell(1,L);
apen_TOT=cell(1,L);
for i=1:L
    apen_SUB{i}=load(listordered(i).name,'apen_all');
    apen_SUB{i}=cellfun(@nanmean,apen_SUB{i}.apen_all)';
    name=listordered(i).name;
    name=name(1:regexp(name,'_','once')-1);
    actualdirectory=pwd;
    cd('..\..\results\SampEn_Kijon_ISI_iaaft_NS30_m3_r1ms_MFR1_L2500')
    list=dir;
    nameT={list.name}';
    for j=1:length(nameT)
        idx=regexp(nameT{j}, '_', 'once');
        if ~isempty(idx) && strcmp(name,nameT{j}(1:idx-1))==1
            break
        end
    end
    apen_TOT{i}=load(nameT{j},'apen_all');
    apen_TOT{i}=cellfun(@nanmean,apen_TOT{i}.apen_all)';
    if strcmp(name,'subcritico')
        rho(i,1)=corr(apen_SUB{i},apen_TOT{i},'Rows','complete');
    elseif strcmp(name,'critico')
        rho(i-L/3,2)=corr(apen_SUB{i},apen_TOT{i},'Rows','complete');
    else
        rho(i-2*(L/3),3)=corr(apen_SUB{i},apen_TOT{i},'Rows','complete');
    end
    cd(actualdirectory)
end
clear apen_SUB apen_TOT nameT list I idx dates name i j
%plot pearson correlation coefficent
figure
plot(1:L/3,rho(:,1),'k--');
hold on
plot(1:L/3,rho(:,2),'k');
plot(1:L/3,rho(:,3),'k-.');
hold off
%title('Pearson correlation coefficent respect to the record length')
ylabel('\rho [a.u.]')
xlabel('Record length [min]')
legend({'subcritical' 'critical' 'supercritical'})
saveppt('Figures.ppt','','-f')

%Normality and Homoscedasticity test
pN=zeros(3,L/3);
hN=zeros(3,L/3);
pV=zeros(1,L/3);
for i=1:L/3
    apen_sub=load(listordered(i).name,'apen_all');
    apen_sub=cellfun(@nanmean,apen_sub.apen_all)';
    apen_cr=load(listordered(i+L/3).name,'apen_all');
    apen_cr=cellfun(@nanmean,apen_cr.apen_all)';
    apen_super=load(listordered(i+2*(L/3)).name,'apen_all');
    apen_super=cellfun(@nanmean,apen_super.apen_all)';
    y=[apen_sub apen_cr apen_super];
    for j=1:3
        x=y(:,j);
        x=x(~isnan(x));
        x=x(~isinf(x));
        x=(x-mean(x))./std(x);
        [hN(j,i), pN(j,i), ~] = swtest(x,0.05);
    end
    pV(i) = vartestn(y,'TestType','BrownForsythe','Display','off');
end
figure
semilogy(1:L/3,pN(1,:),'k');
hold on
semilogy(1:L/3,pN(2,:),'k--');
semilogy(1:L/3,pN(3,:),'k-.');
plot(1:L/3,repmat(0.05,L/3),'k--','LineWidth',2)
ylim([0.001 1])
ylabel('p-value (nomrality test) [a.u.]')
yyaxis right
ylabel('p-value (homoscedasticity test) [a.u.]')
plot(1:L/3,pV);
ylim([0 1])
xlabel('Record length [min]')
legend({'subcritical' 'critical' 'supercritical','confidence level'})
title({'Normality test (Shapiro-Wilk)', 'Homoscedasticity test (Brown-Forsythe)'}')
saveppt('Figures.ppt','','-f')

%  Mann-Whitney U-test + Bonferroni correction
p_sub_cr=zeros(1,L/3);
p_sub_super=zeros(1,L/3);
p_cr_super=zeros(1,L/3);
for i=1:L/3
    apen_sub=load(listordered(i).name,'apen_all');
    apen_sub=cellfun(@nanmean,apen_sub.apen_all)';
    apen_cr=load(listordered(i+L/3).name,'apen_all');
    apen_cr=cellfun(@nanmean,apen_cr.apen_all)';
    apen_super=load(listordered(i+2*(L/3)).name,'apen_all');
    apen_super=cellfun(@nanmean,apen_super.apen_all)';
    p_sub_cr(i)=ranksum(apen_sub,apen_cr)*3;
    p_sub_super(i)=ranksum(apen_sub,apen_super)*3;
    p_cr_super(i)=ranksum(apen_cr,apen_super)*3;
end
clear apen_sub apen_cr apen_super
figure
semilogy(1:L/3,p_sub_cr,'k');
hold on
semilogy(1:L/3,p_sub_super,'k--');
semilogy(1:L/3,p_cr_super,'k-.');
plot(1:L/3,repmat(0.05,L/3),'k--','LineWidth',2)
hold off
ylim([0.0001 1])
ylabel('p-value [a.u.]')
xlabel('Record length [min]')
legend({'subcritical vs. critical' 'subcritical vs. supercritical' 'critical vs. supercritical'})
%title('Mann-Whitney U-test + Bonferroni correction')
saveppt('Figures.ppt','','-f')


% % Dunn
% p_sub_cr=zeros(1,L/3);
% p_sub_super=zeros(1,L/3);
% p_cr_super=zeros(1,L/3);
% for i=1:L/3
%     apen_sub=load(listordered(i).name,'apen_all');
%     apen_sub=cellfun(@nanmean,apen_sub.apen_all)';
%     apen_cr=load(listordered(i+L/3).name,'apen_all');
%     apen_cr=cellfun(@nanmean,apen_cr.apen_all)';
%     apen_super=load(listordered(i+2*(L/3)).name,'apen_all');
%     apen_super=cellfun(@nanmean,apen_super.apen_all)';
%     y=[apen_sub apen_cr apen_super];
%     S=size(y,2);
%     couples=combnk(1:3,2);
%     group=repmat(1:size(y,2),size(y,1),1);
%     Y=y(:);
%     group=group(:);
%     [q, q_crt]=dunn(Y(~isnan(Y))',group(~isnan(Y))');
%     q=[q{:,2}];
%     p_sub_cr(i)=q(1);
%     p_sub_super(i)=q(3);
%     p_cr_super(i)=q(2);
% end
% clear apen_sub apen_cr apen_super Y group
% figure 
% plot(1:L/3,p_sub_cr,'k');
% hold on
% plot(1:L/3,p_sub_super,'k--');
% plot(1:L/3,p_cr_super,'k-.');
% plot(1:L/3,repmat(q_crt,L/3),'k--','LineWidth',2)
% hold off
% set(gca,'Ydir','reverse')
% ylabel('q-value [a.u.]')
% xlabel('Record length [min]')
% legend({'subcritical vs. critical' 'subcritical vs. supercritical' 'critical vs. supercritical'})
% title('stepdown Dunn test for multiple comparison')
% saveppt('Figures.ppt','','-f')