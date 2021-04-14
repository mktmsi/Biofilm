clear

FlagSingle=0;
flagInhibited=0;
rowEnd=1000;
%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=2;
XlimStart=2.0;
XlimEnd=3.1;

NumOfGyouInSubPlot=4;
NumOfRetsuInSubPlot=6;
IndexOfSubPlot=0;
%a=[36 22 8];
%a=[91 71 51 31 11];
a=[86 76 66 56 46 36 26 16 6];
a=[22 16 10 4];
a=[43 31 19 7]

for NN=1:1:7
    FieldSize=[51,61,71,81,91,101,111];
    for BossSet=1:1:1
        Results=[];
        if flagInhibited==1
            f1Name=strcat('Signal_PhaseDifferenceFieldSize',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f2Name=strcat('Signal_Histogram',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f3Name=strcat('Signal_WaveForms',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f5Name=strcat('Signal_Heatmap',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f1=figure('Name',f1Name,'Position',[0 0 5000 5000]);
            f2=figure('Name',f2Name,'Position',[0 0 5000 5000]);
            f3=figure('Name',f3Name,'Position',[0 0 5000 5000]);
            f5=figure('Name',f5Name,'Resize','on');
        else
            f1Name=strcat('Cell_PhaseDifferenceFieldSize',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f2Name=strcat('Cell_Histogram',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f3Name=strcat('Cell_WaveForms',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f5Name=strcat('Cell_Heatmap',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet));
            f1=figure('Name',f1Name,'Position',[0 0 5000 5000]);
            f2=figure('Name',f2Name,'Position',[0 0 5000 5000]);
            f3=figure('Name',f3Name,'Position',[0 0 5000 5000]);
            f5=figure('Name',f5Name,'Resize','on');
        end
        
        folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/2D/Hex/Paper/Now/0918/WithNoise/N=';
        folderName0=strcat(folderName0,num2str(FieldSize(NN)));
        folderName0=strcat(folderName0,'/',num2str(BossSet));
        Results=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
        for ParameterSet=0:1:23
            ParameterSet
            if flagInhibited==1
                f6Name=strcat('Signal_WaveForm',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet),'_Paraset',num2str(ParameterSet));
            else
                f6Name=strcat('Cell_WaveForm',num2str(FieldSize(NN)),'_BossSet',num2str(BossSet),'_Paraset',num2str(ParameterSet));
            end
            
            f6=figure('Name',f6Name,'Resize','on','Position',[0 0 5000 5000]);
            q=fix(ParameterSet/NumOfRetsuInSubPlot);%q:商  r:余り
            r=rem(ParameterSet,NumOfRetsuInSubPlot);
            IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-ParameterSet);
            
            folderName=strcat(folderName0,'/Results/ParameterSet');
            folderName=strcat(folderName,num2str(ParameterSet))
            cd(folderName);
            
            if flagInhibited==1
                fileName=strcat(folderName,'/InhibitedCellcounts');
            else
                fileName=strcat(folderName,'/Cellcounts');
            end
            
            fileName=strcat(fileName,num2str(ParameterSet));
            fileName=strcat(fileName,'.csv');
            fileData=[]
            fileData=csvread(fileName);
            
            %%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
            %fileData=fileData(1:2401,:);
            [row,col]=size(fileData);
            start=2;
            dend=row;
            b=cast((col-2)/2,'int64');
            DataTop=[];
            DataLeft=[];
            DataRight=[];
            DataTop=fileData(start:dend,2);
            DataLeft=fileData(start:dend,3);
            %DataRight=fileData(start:dend,4);
            %XData=(fileData(start:dend,1)/1000000);%ALIFE2020
            XData=(fileData(start:dend,1)/10000);%Fullpaper
            [row,col]=size(DataTop);
            if flagInhibited==0
                dDataLeft=diff(DataLeft);
                dDataTop=diff(DataTop);
                [row,col]=size(dDataTop);
                
                if min(dDataTop)<0.0
                    dDataTop=dDataTop+abs(min(dDataTop));
                end
                if min(dDataLeft)<0.0
                    dDataLeft=dDataLeft+abs(min(dDataLeft));
                end
            end
            
            if flagInhibited==1
                PlotLeft=DataLeft(start:row-1);
                PlotTop=DataTop(start:row-1);
            else
                PlotLeft=dDataLeft(start:row-1);
                PlotTop=dDataTop(start:row-1);
            end
            
            threshPeak=0.15*min([abs(max(PlotTop)) abs(max(PlotLeft))])
            figure(f6);
            %極大値の値とそのインデックスを取得
            [pks1,locsLeft] = findpeaks(PlotLeft,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            %[pks1,TF1] = findpeaks(FilterdDataLeft,'MinPeakProminence',0.05e+0);
            % [pks2,locsRight] = findpeaks(DataRight(start:row-1),'MinPeakHeight',15);
            [pks3,locsTop] = findpeaks(PlotTop,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            % [pks2,TF2] = findpeaks(FilterdDataRight,'MinPeakProminence',0.05e+0);
            plot(PlotTop,'LineWidth',LineHaba)
            hold on;
            findpeaks(PlotTop,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            hold on;
            plot(PlotLeft,'LineWidth',LineHaba)
            hold on;
            findpeaks(PlotLeft,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            %グラフの体
            set(gca,'FontSize',MemoriSize);
            set(gca,'LineWidth',LineHaba);%軸の太さ
            %グラフ保存
            filenameBlend=strcat(folderName,'/',f6Name,num2str(ParameterSet),'.png');
            saveas(gcf,filenameBlend);
            hold off;
            close
            
            figure(f3);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            plot(PlotTop,'LineWidth',LineHaba)
            hold on;
            plot(PlotLeft,'LineWidth',LineHaba)
            hold on;
            findpeaks(PlotTop,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            hold on;
            findpeaks(PlotLeft,'MinPeakHeight',threshPeak,'MinPeakProminence',threshPeak);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            
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
            if isempty(locsTop)==1 | isempty(locsLeft)==1  |length(locsTop)==1|length(locsLeft)==1
                
            else
                for i=1:row-1%時間
                    %%%%%%%%%Top%%%%%%%%%
                    if  mTop == 1
                        if i<locsTop(mTop)% 1個目の極大値よりも前の時
                            Top(i)=NaN;
                        else
                            mTop=mTop+1;
                            TopPre=locsTop(mTop-1);
                            TopNext=locsTop(mTop);
                            Bunsi=cast((i-TopPre),'double');
                            Bumbo=cast((TopNext-TopPre),'double');
                            Top(i)=(Bunsi/Bumbo)*2.0*pi;
                        end
                    else
                        if mTop<length(locsTop)
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
                        if mLeft<length(locsLeft)
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
            figure(f1);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            plot(TopAndLeft,'o')
            ylim([0.0 2.0*pi]);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            
            figure(f2);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            Hist=histogram(TopAndLeft)
            Hist.BinLimits=([0.0 2.0*pi]);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            
            [M,I]=max(Hist.Values);%ヒストグラムの最大値のインデックスIを求める
            %Results(IndexOfSubPlot)=Hist.BinEdges(I)+Hist.BinWidth/2.0;
            Results(NumOfGyouInSubPlot-q,r+1)=Hist.BinEdges(I)+Hist.BinWidth/2.0;
        end
        
        
        if flagInhibited==1
            fileName1=strcat(folderName0,'/Signal_PhaseDifferences.png')
            fileName2=strcat(folderName0,'/Signal_Histogram.png');
            fileName3=strcat(folderName0,'/Signal_WaveForms.png');
            fileName5=strcat(folderName0,'/Signal_Heatmap',num2str(BossSet),'.png');
        else
            fileName1=strcat(folderName0,'/Cells_PhaseDifferences.png')
            fileName2=strcat(folderName0,'/Cells_Histogram.png');
            fileName3=strcat(folderName0,'/Cells_WaveForms.png');
            fileName5=strcat(folderName0,'/Cells_Heatmap',num2str(BossSet),'.png');            
        end
        
        j1=figure(f1);
        %fileName1=strcat(folderName0,'/PhaseDifferences.png')
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        set(j1,'Units','Inches');
        pos = get(j1,'Position');
        set(j1,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        saveas(gcf,fileName1);
        close
        
        j2= figure(f2);
        %fileName2=strcat(folderName0,'/Histogram.png');
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        set(j2,'Units','Inches');
        pos = get(j2,'Position');
        set(j2,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        saveas(gcf,fileName2); 
        close
        
        j3=figure(f3);
        %fileName3=strcat(folderName0,'/WaveForms.png');
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        set(j3,'Units','Inches');
        pos = get(j3,'Position');
        set(j3,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
        saveas(gcf,fileName3);
        close
        
        figure(f5);
        mymap = [1 1 1
            0.75 0.75 0.75
            0.5 0.5 0.5
            0.25 0.25 0.25
            0 0 0
            0.25 0.25 0.25
            0.5 0.5 0.5
            0.75 0.75 0.75
            1 1 1];
        %h=heatmap(reshape(Results,[NumOfGyouInSubPlot,NumOfRetsuInSubPlot]))
        h=heatmap(Results);
        colormap(mymap);
        %fileName5=strcat(folderName0,'/Heatmap',num2str(BossSet),'.png');
        saveas(gcf,fileName5);
        close
    end
end




