clear

FlagSingle=0;
rowEnd=1000;
%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=2;
DotHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
NumOfGyouInSubPlot=2;
NumOfRetsuInSubPlot=7;

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
%a=[37 29 21 13 5];%%gyou5 retsu4
a=[22 8];%%gyou2 retsu7
SSSS=1;
for BossSet=SSSS:-1:1
    PhaseDifferenceJudge=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
    
    SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
    SubPlotHistogram=figure('Name','Histogram','Position',[0 0 5000 5000]);
    
    DisSize=5000;
    SubPlotWaveForm=figure('Name','Wave forms','Position',[0 0 DisSize DisSize]);
    SinglePlotWaveForm=figure('Name','Wave forms','Resize','on','Position',[0 00 DisSize DisSize]);
    folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/DrugInstall/normalized/ParameterModified/ALIFE追試/MikamiWithNoise/Working/Results1';
    %folderName0=strcat(folderName0,'/',num2str(BossSet));
    minDist=2;
    maxDist=12;
    Sets=12;
    
    for ParameterSet=Sets:-1:0
        ThrePeak=1.0;
        PeakPro=100.0;
        ParameterSet;
        q=fix(ParameterSet/NumOfRetsuInSubPlot);%q:商  r:余り
        IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-ParameterSet);
        
        folderName=strcat(folderName0,'/ParameterSet');
        folderName=strcat(folderName,num2str(ParameterSet));
        cd(folderName);
        fileName=strcat(folderName,'/SumV');
        %fileName=strcat(fileName,num2str(ParameterSet));
        fileName=strcat(fileName,'.csv')
        
        fileData=[];
        fileData=csvread(fileName);
        
        %%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
        [row,col]=size(fileData);
        start=1;
        dend=row-1;
        b=cast((col-2)/2,'int64');
        DataTop=[];
        DataLeft=[];
        DataRight=[];
        DataTop=fileData(start:dend,2);
        DataLeft=fileData(start:dend,3);
        XData=(fileData(start:dend,1)/10000);%Fullpaper
        dDataLeft=diff(DataLeft);
        dDataTop=diff(DataTop);
        
        [row,col]=size(dDataTop);
        
        figure(SinglePlotWaveForm);
        %極大値の値とそのインデックスを取得
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        %yyaxis left
        plot(locsLeft,pks1,'b*','LineWidth',DotHaba);
        hold on;
        %yyaxis left
        plot(dDataLeft,'-.','color','#0072BD','LineWidth',LineHaba);
        hold on;
        %yyaxis left
        plot(dDataTop,'-','color','#A2142F','LineWidth',LineHaba)
        hold on;
        %yyaxis left
        plot(locsTop,pks3,'r*','LineWidth',DotHaba);
        hold on;
        %yyaxis right
        %plot(Distance,'LineWidth',LineHaba);
        %グラフの体
        xlim([1 row-1]);
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%軸の太さ
        %グラフ保存
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'.png');
        saveas(gcf,filenameBlend);
        hold off;
        %clf
        
        
        figure(SubPlotWaveForm);
        subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
        plot(dDataTop,'LineWidth',LineHaba)
        hold on;
        plot(dDataLeft,'LineWidth',LineHaba)
        titlename=strcat('i=',num2str(ParameterSet));
        title(titlename);
        hold on;
        
        %%%%バイオフィルム間の位相差を求める
        Top=[];
        Left=[];
        Right=[];
        mTop=1;
        mRight=1;
        mLeft=1;%足をついた回数
        Top=[];
        Right=[];
        Left=[];
        if isempty(locsTop)==0&&isempty(locsLeft)==0
            for i=1:row-1%時間
                %%%%%%%%%Top%%%%%%%%%
                if  mTop == 1 % mTop番目のPeak
                    if i<locsTop(mTop)% 1個目の極大値よりも前の時
                        Top(i)=NaN;
                    else
                        i;
                        mTop=mTop+1;
                        TopPre=locsTop(mTop-1);
                        TopNext=locsTop(mTop);
                        Bunsi=cast((i-TopPre),'double');
                        Bumbo=cast((TopNext-TopPre),'double');
                        Top(i)=(Bunsi/Bumbo)*2.0*pi;
                    end
                else
                    if i<locsTop(length(locsTop)) %ステップ数が最後のpeakがあるステップよりも小さい時
                        if (locsTop(mTop)<i)
                            mTop=mTop+1;
                        end
                        TopPre=locsTop(mTop-1);
                        TopNext=locsTop(mTop);
                        Bunsi=cast((i-TopPre),'double');
                        Bumbo=cast((TopNext-TopPre),'double');
                        Top(i)=(Bunsi/Bumbo)*2.0*pi;
                        
                    else
                        Top(i)=NaN;
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%Top 終わり%%%%%%%%%%
                %%%%%%%%%Left%%%%%%%%%
                if  mLeft == 1
                    if i<locsLeft(mLeft)% 1個目の極大値よりも前の時
                        Left(i)=NaN;
                    else
                        mLeft=mLeft+1;
                        LeftPre=locsLeft(mLeft-1);
                        LeftNext=locsLeft(mLeft);
                        Bunsi=cast((i-LeftPre),'double');
                        Bumbo=cast((LeftNext-LeftPre),'double');
                        Left(i)=(Bunsi/Bumbo)*2.0*pi;
                    end
                else
                    if i<locsLeft(length(locsLeft)) %ステップ数が最後のpeakがあるステップよりも小さい時
                        if (locsLeft(mLeft)<i)
                            mLeft=mLeft+1;
                        end
                        LeftPre=locsLeft(mLeft-1);
                        LeftNext=locsLeft(mLeft);
                        Bunsi=cast((i-LeftPre),'double');
                        Bumbo=cast((LeftNext-LeftPre),'double');
                        Left(i)=(Bunsi/Bumbo)*2.0*pi;
                        
                    else
                        Left(i)=NaN;
                    end
                end
            end
            Top=transpose(Top);
            Left=transpose(Left);
            TopAndLeft=[];
            TopAndLeft=Left-Top;
            for i=1:length(TopAndLeft)
                if abs(TopAndLeft(i))>2.0*pi
                    TopAndLeft(i)=mod(TopAndLeft(i),2.0*pi);
                end
                if TopAndLeft(i)<0.0
                    TopAndLeft(i)=TopAndLeft(i)+2.0*pi;
                end
            end
            [Num,col]=size(TopAndLeft);
            
            if any(TopAndLeft(:)>0.0)~=0
                BossSet
                ParameterSet
            end
            
            figure(SubPlotWaveForm);
            %yyaxis left
            plot(TopAndLeft,'o')
            %xlim([XLlimMin XLlimMax]);
            ylim([0 2.0*pi]);
            xlim([1 row-1]);

            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            filenameBlend=strcat(folderName,'/PhaseDiference',num2str(ParameterSet),'.png');
            saveas(gcf,filenameBlend);
            %clf
            
            
            figure(SubPlotPhaseDiffernce);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            %yyaxis left
            plot(TopAndLeft,'o')
            ylim([0 2.0*pi]);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            
           
            figure(SubPlotHistogram);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            histogram(TopAndLeft);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            
        end
    end
    
    
    figure(SubPlotPhaseDiffernce);
    fileName1=strcat(folderName0,'/SubplotPhaseDifferences.png')
    saveas(gcf,fileName1);
    %pause(1)
    close
    
    figure(SubPlotHistogram);
    fileName2=strcat(folderName0,'/Histogram.png');
    saveas(gcf,fileName2);
    %pause(1)
    close
    
    figure(SubPlotWaveForm);
    close
    figure(SinglePlotWaveForm);
    close
    
    
    
end

