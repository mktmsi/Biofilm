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
NumOfGyouInSubPlot=7;
NumOfRetsuInSubPlot=10;

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
%a=[37 29 21 13 5];%%gyou5 retsu4
for BossSet=1:1:1
    PhaseDifferenceJudge=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
    f1=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
    f2=figure('Name','Histogram','Position',[0 0 5000 5000]);
    f3=figure('Name','Wave forms','Position',[0 0 500 500]);
    f4=figure('Name','Wave forms','Resize','on','Position',[0 00 5000 500]);
    %folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/2D/Hex/ChangeK2AndAlpha';
    folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/2D/Hex/ChangeK2AndAlpha_Revised/DzConst0_0125Take2';
    %folderName0=strcat(folderName0,'/',num2str(BossSet));
    %ThrePeaks=[200 200 200 200 200 200 200 200 200 230 230 230 230 230 230 230 230 230 230 230 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260 260];
    %PeakPros=[200 200 200 200 200 200 200 200 200 230 230 230 230 230 230 230 230 230 230 230 260 260 260 260 200 200 200 250 250 250 250 250 250 250 250 250 250 250 250 250];
                minDist=2;
            maxDist=12;
    Sets=69;
    for ParameterSet=Sets:-1:0
        if ParameterSet>=60
            ThrePeak=500.0;%ThrePeaks(ParameterSet+1);
            PeakPro=500.0;%PeakPros(ParameterSet+1);
        elseif ParameterSet>=50
            ThrePeak=400.0;
            PeakPro=400.0; 
        elseif ParameterSet>=40
            ThrePeak=300.0;
            PeakPro=300.0; 
        elseif ParameterSet>=30
            ThrePeak=200.0;
            PeakPro=200.0;
        elseif ParameterSet>=20
            ThrePeak=250.0;
            PeakPro=250.0; 
        elseif ParameterSet>=10
            ThrePeak=200.0;
            PeakPro=200.0; 
        else 
            ThrePeak=200.0;
            PeakPro=200.0; 
        end
            ParameterSet
        %ThrePeak=ThrePeaks(ParameterSet+1);
        %PeakPro=PeakPros(ParameterSet+1);
        q=fix(ParameterSet/NumOfRetsuInSubPlot);%q:商  r:余り
        IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-ParameterSet);
        
        folderName=strcat(folderName0,'/Results/ParameterSet')
        folderName=strcat(folderName,num2str(ParameterSet));
        cd(folderName);
        fileName=strcat(folderName,'/Cellcounts');
        fileName=strcat(fileName,num2str(ParameterSet));
        fileName=strcat(fileName,'.csv');
        fileData=[];
        fileData=csvread(fileName);
        
        %%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
        %fileData=fileData(1:2401,:);
        [row,col]=size(fileData);
        start=2;
        dend=row-1;
        b=cast((col-2)/2,'int64');
        DataTop=[];
        DataLeft=[];
        DataRight=[];
        DataTop=fileData(start:dend,2);
        DataLeft=fileData(start:dend,3);
        Distance=fileData(start:dend,4);
        %DataRight=fileData(start:dend,4);
        %XData=(fileData(start:dend,1)/1000000);%ALIFE2020
        XData=(fileData(start:dend,1)/10000);%Fullpaper
        dDataLeft=diff(DataLeft);
        dDataTop=diff(DataTop);
        
        [row,col]=size(dDataTop);
        
        figure(f4);
        %極大値の値とそのインデックスを取得
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        %[pks1,TF1] = findpeaks(FilterdDataLeft,'MinPeakProminence',0.05e+0);
        % [pks2,locsRight] = findpeaks(DataRight(start:row-1),'MinPeakHeight',15);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        %[pks2,TF2] = findpeaks(FilterdDataRight,'MinPeakProminence',0.05e+0);
        yyaxis left
        plot(locsLeft,pks1,'b*','LineWidth',DotHaba);
        hold on;
        yyaxis left
        plot(dDataLeft,'-.','color','#0072BD','LineWidth',LineHaba);
        hold on;
        yyaxis left
        plot(dDataTop,'-','color','#A2142F','LineWidth',LineHaba)
        hold on;
        yyaxis left
        plot(locsTop,pks3,'r*','LineWidth',DotHaba);
        hold on;
        yyaxis right
        plot(Distance,'LineWidth',LineHaba);
        %グラフの体
        xlim([0 row-1]);
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%軸の太さ
        %グラフ保存
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'.png');
        saveas(gcf,filenameBlend);
        hold off;
        clf
        %{
        figure(f3);
        subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
        plot(dDataTop,'LineWidth',LineHaba)
        hold on;
        plot(dDataLeft,'LineWidth',LineHaba)
        titlename=strcat('i=',num2str(ParameterSet));
        title(titlename);
        hold on;
        %}
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
                    %{
                    if mTop<=length(locsTop)% mTopの値がPeak数を超えていなければ
                        if (locsTop(mTop)<i)&&mTop<length(locsTop)
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
                    %}
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
                    %{
                    if mLeft<=length(locsLeft)
                        if (locsLeft(mLeft)<i  && mLeft<length(locsLeft))
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
                    %}
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
            %f1=figure('Name','Phase Differnce');
           %{

            Distance2=Distance;
            Distance2(Distance2<=minDist)=NaN;
            Distance2(Distance2>=maxDist)=NaN;
            Distance2(Distance2>=minDist & Distance2<=maxDist)=1;
            XData2=XData.*Distance2;
            %}
            Distance=Distance(1:Num,1);
            XLlimMin=min(find(Distance<=maxDist));
            XLlimMax=max(find(Distance>=minDist));
            figure(f3);
            yyaxis left
            plot(TopAndLeft,'o')
            xlim([XLlimMin XLlimMax]);
            ylim([0 2.0*pi]);
            yyaxis right
            plot(Distance,'LineWidth',LineHaba);
            xlim([XLlimMin XLlimMax]);
            %xlim([0 Num]);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            filenameBlend=strcat(folderName,'/PhaseDiference',num2str(ParameterSet),'.png');
            saveas(gcf,filenameBlend);
            clf
            %{
            figure(f1);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            yyaxis left
            plot(TopAndLeft,'o')
            ylim([0 2.0*pi]);
            xlim([XLlimMin XLlimMax]);
            hold on;
            yyaxis right
            plot(Distance,'LineWidth',LineHaba);
            xlim([XLlimMin XLlimMax]);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            
            figure(f2);
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            histogram(TopAndLeft);
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            hold on;
            %}
            Gyou=NumOfGyouInSubPlot-fix(ParameterSet/NumOfRetsuInSubPlot);
            Retsu=rem(ParameterSet,NumOfRetsuInSubPlot)+1;
            if (TopAndLeft(XLlimMax)<= 1.5*pi) && (TopAndLeft(XLlimMax)>= 0.5*pi)
                PhaseDifferenceJudge(Gyou,Retsu)=1.0;
            else
                PhaseDifferenceJudge(Gyou,Retsu)=0.0;
            end
            
        end
    end
    figure(f1);
    fileName1=strcat(folderName0,'/PhaseDifferences.png')
    saveas(gcf,fileName1);
    %pause(1)
    close
    
    figure(f2);
    fileName2=strcat(folderName0,'/Histogram.png');
    saveas(gcf,fileName2);
    %pause(1)
    close
    %{
    figure();
    %PhaseDifferenceJudge=reshape(PhaseDifferenceJudge,[NumOfGyouInSubPlot NumOfRetsuInSubPlot] );
    heatmap(PhaseDifferenceJudge);
    fileName2=strcat(folderName0,'/DifferencesHeatmap.png');
    saveas(gcf,fileName2);
    close
    %}
    %{
    figure(f3);
    fileName3=strcat(folderName0,'/WaveForms.png');
    saveas(gcf,fileName3);
    pause(1)
    close
    %}
    
    
end

