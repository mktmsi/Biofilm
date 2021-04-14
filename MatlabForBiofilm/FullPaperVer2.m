clear

FlagSingle=0;
rowEnd=1000;
%%ã°ã©ãã?®å¶å¾¡
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=2;
DotHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
NumOfGyouInSubPlot=3;
NumOfRetsuInSubPlot=391;

ModeXAxis=0;%0:æ¨ªè»¸æéã?1:æ¨ªè»¸ããã?

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
%a=[37 29 21 13 5];%%gyou5 retsu4
%a=[22 8];%%gyou2 retsu7
a=[16 10 4];%%gyou3 retsu3

DisSize=5000;
SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
SubPlotHistogram=figure('Name','Histogram','Position',[0 0 5000 5000]);
SubPlotWaveForm=figure('Name','Wave forms','Position',[0 0 DisSize DisSize]);
SinglePlotWaveForm1=figure('Name','Wave forms All','Resize','on','Position',[0 00 DisSize DisSize]);
SinglePlotWaveForm2=figure('Name','Wave forms 3000','Resize','on','Position',[0 00 DisSize DisSize]);

SSSS=1;
for BossSet=SSSS:-1:1
    PhaseDifferenceJudge=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
    
    folderName00='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Working/Division';
    folderName0=strcat(folderName00,'/',num2str(BossSet),'/Results');
    minDist=2;
    maxDist=12;
    Sets=1;
    for ParameterSet=Sets:-1:0
        ThrePeak=1.0;
        PeakPro= 3.0;%interval100 ã®ã¨ã? 3.0   å¾å 25.0 paraset=20-8ã®æ?:10.0 3-8ã®æ?:5.0
        %interval100ã®æï¼Bossset=1-6ã¾ã§ã¯ãã£ã¨3.0ã§ok
        %{
        PeakPro= 100.0;
        if BossSet >8 && BossSet <16%interval1000ã®æ?
            PeakPro= 500.0;
        else
            PeakPro= 500.0;
        end
        %}
        %{
        if BossSet <6
            PeakPro= 3.0;
        else
            if ParameterSet>=20
                PeakPro=25.0;
            elseif ParameterSet>=8
                PeakPro=10.0;
            else
                PeakPro=5.0;
            end
        end
        %}
        if BossSet >=23
            PeakPro=10.0;
            if BossSet==23&&ParameterSet>=32
                PeakPro=15.0;
            end
        elseif BossSet==11
            if ParameterSet>19
                PeakPro=50.0;
            elseif ParameterSet>10
                PeakPro=25.0;
            else
                PeakPro=10.0;
            end
        end
        
        if ParameterSet>=20
            tempParaSet=ParameterSet-20;
        else
            tempParaSet=ParameterSet;
        end
        q=fix(tempParaSet/NumOfRetsuInSubPlot);%q:å?  r:ä½ã
        IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-tempParaSet);
        
        folderName=strcat(folderName0,'/ParameterSet');
        folderName=strcat(folderName,num2str(ParameterSet));
        folderName='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Working/Division/base/Results/ForPaperize/2/ParameterSet0';
        cd(folderName);
        
        fileName=strcat(folderName,'/SumV');
        %fileName=strcat(fileName,num2str(ParameterSet));
        fileName=strcat(fileName,'.csv')
        
        fileData=[];
        fileData=csvread(fileName);
        
        %%%2ã¤ã®ãã¤ãªãã£ã«ã?ã®ãµã¤ãºã®æéæ¨ç§»ãè¡¨ç¤º%%%%%%%%%%%
        [row,col]=size(fileData);
        start=1;
        dend=row-1;
        b=cast((col-2)/2,'int64');
        DataTop=[];
        DataLeft=[];
        DataRight=[];
        DataTop=fileData(start:dend,2);
        DataLeft=fileData(start:dend,3);
        Distance=fileData(start:dend,4);
        %
        if ModeXAxis==0
            XData=(fileData(start:dend,1)/10000);%Fullpaper
        else
            XData=(fileData(start:dend,4));
        end
        
        dDataLeft=diff(DataLeft);%å¾®å?
        dDataTop=diff(DataTop);
        
        [row,col]=size(dDataTop);
        
        figure(SinglePlotWaveForm1);
        colororder({'#0072BD','#238C2A'})
        %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
        %æ¥µå¤§å¤ã®å¤ã¨ãã?®ã¤ã³ã?ã?ã¯ã¹ãåå¾?
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
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
        plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
        %ã°ã©ãã?®ä½?
        xlim([1 row-1]);
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%è»¸ã®å¤ªã?
        %ã°ã©ãä¿å­?
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'.png');
        saveas(gcf,filenameBlend);
        hold off;
        clf
        
        
        
        figure(SinglePlotWaveForm2);
        colororder({'#0072BD','#238C2A'})
        %æ¥µå¤§å¤ã®å¤ã¨ãã?®ã¤ã³ã?ã?ã¯ã¹ãåå¾?
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
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
        plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
        %ã°ã©ãã?®ä½?
        if locsTop(length(locsTop))>8000
        xlim([8000 row-1]);
        end
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%è»¸ã®å¤ªã?
        %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
        %ã°ã©ãä¿å­?
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'3000-.png');
        saveas(gcf,filenameBlend);
        hold off;
        clf
        
        %{
        figure(SubPlotWaveForm);
        subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
        colororder({'#0072BD','#238C2A'})
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
        plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
        %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
        Du=(BossSet)*0.2;
        %}
        %{
        if ParameterSet>=20
            aaa=' uin=2100-4000';
        else
            aaa=' uin=100-2000';
        end
        
        sgtitleName=strcat('Du=',num2str(Du),aaa);
        
        sgtitle(sgtitleName)
        %}
        titlename=strcat('i=',num2str(ParameterSet));
        title(titlename);
        hold on;
        
        
        %%%%ãã¤ãªãã£ã«ã?éã?®ä½ç¸å·®ãæ±ãã?
        Top=[];
        Left=[];
        Right=[];
        mTop=1;
        mRight=1;
        mLeft=1;%è¶³ãã¤ã?ãåæ°
        Top=[];
        Right=[];
        Left=[];
        if isempty(locsTop)==0&&isempty(locsLeft)==0&&length(locsTop)>2&&length(locsLeft)>2
            for i=1:row-1%æé
                %%%%%%%%%Top%%%%%%%%%
                if  mTop == 1 % mTopçªç®ã®Peak
                    if i<locsTop(mTop)% 1åç®ã®æ¥µå¤§å¤ãããåã®æ?
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
                    if i<locsTop(length(locsTop)) %ã¹ã?ã?ãæ°ãæå¾ã?®peakãããã¹ã?ã?ããããå°ãã?æ?
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
                %%%%%%%%%%%%%%%%%%%%%%%%Top çµãã?%%%%%%%%%%
                %%%%%%%%%Left%%%%%%%%%
                if  mLeft == 1
                    if i<locsLeft(mLeft)% 1åç®ã®æ¥µå¤§å¤ãããåã®æ?
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
                    if i<locsLeft(length(locsLeft)) %ã¹ã?ã?ãæ°ãæå¾ã?®peakãããã¹ã?ã?ããããå°ãã?æ?
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
                if i<locsTop(3)||i<locsLeft(3)% 1åç®ã®æ¥µå¤§å¤ãããåã®æ?
                    TopAndLeft(i)=NaN;
                end
            end
            [Num,col]=size(TopAndLeft);
            
            if any(TopAndLeft(:)>0.0)~=0
                BossSet
                ParameterSet
            end
            
            figure();
            if ModeXAxis==0
                colororder({'#0072BD','#238C2A'})
                yyaxis left
                plot(TopAndLeft,'o')
                ylim([0 2.0*pi]);
                hold on;
                yyaxis right
                plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
                %xlim([XLlimMin XLlimMax]);
                xlim([1 row-1]);
                %set(gca,'defaultAxesColorOrder',[[0.0 0.0 0.0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
            else
                colororder({'#0072BD','#238C2A'})
                XData=XData(1:length(TopAndLeft));
                plot(XData,TopAndLeft,'o','MarkerSize',3)
                ax=gca
                ylim([0 2.0*pi]);
                ax.XDir='reverse';
                ax.YDir='reverse';
                xlim([min(XData) max(XData)]);
            end
            titlename=strcat('i=',num2str(ParameterSet));
            title(titlename);
            if ModeXAxis==0
                filenameBlend=strcat(folderName,'/PhaseDiference',num2str(ParameterSet),'.png');
            else
                filenameBlend=strcat(folderName,'/PhaseDiferenceDistance',num2str(ParameterSet),'.png');
            end
            saveas(gcf,filenameBlend);
            close
            
            
            figure(SubPlotPhaseDiffernce);
            colororder({'#0072BD','#238C2A'})
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            if ModeXAxis==0
                
                yyaxis left
                plot(TopAndLeft,'o')
                ylim([0 2.0*pi]);
                hold on;
                yyaxis right
                plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
                xlim([0 row-1]);
                
            else
                XData=XData(1:length(TopAndLeft));
                plot(TopAndLeft,XData,'o','MarkerSize',3)
                ax=gca;
                xlim([0 2.0*pi]);
                ylim([min(XData) max(XData)]);
                ax.YDir='reverse';
                ax.XDir='reverse';
                view(90,90);
            end
            Du=(BossSet)*0.2;
            %{
            if ParameterSet>=20
                aaa=' uin=2100-4000';
            else
                aaa=' uin=100-2000';
            end
            
            sgtitleName=strcat('Du=',num2str(Du),aaa);
            sgtitle(sgtitleName)
            %}
            %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
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
        
        if ParameterSet==20||ParameterSet==0
            if ModeXAxis==0
                if ParameterSet==20
                    %aaa='later';
                    aaa=strcat('later',num2str(BossSet))
                else
                    %aaa='faster';
                    aaa=strcat('faster',num2str(BossSet))
                end
            else
                if ParameterSet==20
                    %aaa='later';
                    aaa=strcat('laterDistance',num2str(BossSet))
                else
                    %aaa='faster';
                    aaa=strcat('fasterDistance',num2str(BossSet))
                end
            end
            figure(SubPlotPhaseDiffernce);
            fileName1=strcat(folderName0,'/SubplotPhaseDifferences',aaa,'.png')
            saveas(gcf,fileName1);
            %pause(1)
            clf
            
            figure(SubPlotHistogram);
            fileName2=strcat(folderName0,'/Histogram',aaa,'.png');
            saveas(gcf,fileName2);
            %pause(1)
            clf
            
            figure(SubPlotWaveForm);
            fileName2=strcat(folderName0,'/WaveForm',aaa,'.png');
            saveas(gcf,fileName2);
            %pause(1)
            clf
        end
        
        if ParameterSet==0
            figure(SubPlotHistogram);
            close
            figure(SubPlotPhaseDiffernce);
            close
            figure(SinglePlotWaveForm1);
            close
            figure(SinglePlotWaveForm2);
            close
            figure(SubPlotWaveForm);
            close
            SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
            SubPlotHistogram=figure('Name','Histogram','Position',[0 0 5000 5000]);
            SubPlotWaveForm=figure('Name','Wave forms','Position',[0 0 DisSize DisSize]);
            SinglePlotWaveForm1=figure('Name','Wave forms All','Resize','on','Position',[0 00 DisSize DisSize]);
            SinglePlotWaveForm2=figure('Name','Wave forms 3000','Resize','on','Position',[0 00 DisSize DisSize]);
        end
    end
end

