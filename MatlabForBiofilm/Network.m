clear

FlagSingle=0;
rowEnd=1000;
%%グラフ�?�制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=2;
DotHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
NumOfGyouInSubPlot=5;
NumOfRetsuInSubPlot=4;

ModeXAxis=0;%0:�������ԁ@1:��������

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
a=[37 29 21 13 5];%%gyou5 retsu4
%a=[22 8];%%gyou2 retsu7

DisSize=5000;
SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
SubPlotHistogram=figure('Name','Histogram','Position',[0 0 5000 5000]);
SubPlotWaveForm=figure('Name','Wave forms','Position',[0 0 DisSize DisSize]);
SinglePlotWaveForm1=figure('Name','Wave forms All','Resize','on','Position',[0 00 DisSize DisSize]);
SinglePlotWaveForm2=figure('Name','Wave forms 3000','Resize','on','Position',[0 00 DisSize DisSize]);


SSSS=1;
for BossSet=SSSS:-1:1
    PhaseDifferenceJudge=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
    
    folderName00='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Network';
    folderName0=strcat(folderName00,'/',num2str(BossSet),'/Results');
    minDist=2;
    maxDist=12;
    Sets=35;
    for ParameterSet=0:1:0
        
        ThrePeak=1.0;
        PeakPro= 010.0;
        %{
        PeakPro= 100.0;
        if BossSet >8 && BossSet <16%interval1000の�?
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
        q=fix(tempParaSet/NumOfRetsuInSubPlot);
        IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-tempParaSet);
        
        folderName=strcat(folderName0,'/ParameterSet');
        folderName=strcat(folderName,num2str(ParameterSet));
        cd(folderName);
        fileName=strcat(folderName,'/V_para');
        %fileName=strcat(fileName,num2str(ParameterSet));
        fileName=strcat(fileName,'.csv')
        
        fileData=[];
        fileData=csvread(fileName);
        
        [row,col]=size(fileData);
        start=1;
        dend=row-1;
        b=cast((col-2)/2,'int64');
        DataTop=[];
        DataLeft=[];
        DataRight=[];
        DataTop=fileData(start:dend,2);
        DataLeft=fileData(start:dend,4);
        %
        if ModeXAxis==0
            XData=(fileData(start:dend,1)/10000);%Fullpaper
        else
            XData=(fileData(start:dend,4));
        end
        
        dDataLeft=diff(DataLeft);%微�?
        dDataTop=diff(DataTop);
        
        [row,col]=size(dDataTop);
        
        figure(SinglePlotWaveForm1);
        colororder({'#0072BD','#238C2A'})
        %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
        
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        plot(locsLeft,pks1,'b*','LineWidth',DotHaba);
        hold on;
        plot(dDataLeft,'-.','color','#0072BD','LineWidth',LineHaba);
        hold on;
        plot(dDataTop,'-','color','#A2142F','LineWidth',LineHaba)
        hold on;
        plot(locsTop,pks3,'r*','LineWidth',DotHaba);
        hold on;       
        %グラフ�?��?
        xlim([1 row-1]);
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%軸の太�?
        %グラフ保�?
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'.png');
        saveas(gcf,filenameBlend);
        hold off;
        clf
        
        
        figure(SinglePlotWaveForm2);
        colororder({'#0072BD','#238C2A'})
        %極大値の値とそ�?�イン�?�?クスを取�?
        [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
        plot(locsLeft,pks1,'b*','LineWidth',DotHaba);
        hold on;
        plot(dDataLeft,'-.','color','#0072BD','LineWidth',LineHaba);
        hold on;
        plot(dDataTop,'-','color','#A2142F','LineWidth',LineHaba)
        hold on;
        plot(locsTop,pks3,'r*','LineWidth',DotHaba);
        hold on;
        %グラフ�?��?
        xlim([30000 row-1]);
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%軸の太�?
        %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
        %グラフ保�?
        filenameBlend=strcat(folderName,'/dExpansion',num2str(ParameterSet),'3000-.png');
        saveas(gcf,filenameBlend);
        hold off;
        clf
        
        figure(SubPlotWaveForm);
        subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
        colororder({'#0072BD','#238C2A'})
        plot(locsLeft,pks1,'b*','LineWidth',DotHaba);
        hold on;
        plot(dDataLeft,'-.','color','#0072BD','LineWidth',LineHaba);
        hold on;
        plot(dDataTop,'-','color','#A2142F','LineWidth',LineHaba)
        hold on;
        plot(locsTop,pks3,'r*','LineWidth',DotHaba);
        Du=(BossSet)*0.2;
        if ParameterSet>=20
            aaa=' uin=2100-4000';
        else
            aaa=' uin=100-2000';
        end
        sgtitleName=strcat('Du=',num2str(Du),aaa);
        sgtitle(sgtitleName)
        titlename=strcat('i=',num2str(ParameterSet));
        title(titlename);
        hold on;
        
        
        %%%%バイオフィル�?間�?�位相差を求め�?
        Top=[];
        Left=[];
        Right=[];
        mTop=1;
        mRight=1;
        mLeft=1;%足をつ�?た回数
        Top=[];
        Right=[];
        Left=[];
        if isempty(locsTop)==0&&isempty(locsLeft)==0&&length(locsTop)>2&&length(locsLeft)>2
            for i=1:row-1%時間
                %%%%%%%%%Top%%%%%%%%%
                if  mTop == 1 % mTop番目のPeak
                    if i<locsTop(mTop)% 1個目の極大値よりも前の�?
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
                    if i<locsTop(length(locsTop)) %ス�?�?プ数が最後�?�peakがあるス�?�?プよりも小さ�?�?
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
                %%%%%%%%%%%%%%%%%%%%%%%%Top 終わ�?%%%%%%%%%%
                %%%%%%%%%Left%%%%%%%%%
                if  mLeft == 1
                    if i<locsLeft(mLeft)% 1個目の極大値よりも前の�?
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
                    if i<locsLeft(length(locsLeft)) %ス�?�?プ数が最後�?�peakがあるス�?�?プよりも小さ�?�?
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
                if i<locsTop(3)||i<locsLeft(3)% 1個目の極大値よりも前の�?
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
            end
            saveas(gcf,filenameBlend);
            close
            
            
            figure(SubPlotPhaseDiffernce);
            colororder({'#0072BD','#238C2A'})
            subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
            if ModeXAxis==0
                %{
                yyaxis left
                plot(TopAndLeft,'o')
                ylim([0 2.0*pi]);
                hold on;
                yyaxis right
                plot(Distance,'color','#238C2A','LineWidth',3.0*LineHaba);
                xlim([0 row-1]);
                %}
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
            if ParameterSet>=20
                aaa=' uin=2100-4000';
            else
                aaa=' uin=100-2000';
            end
            sgtitleName=strcat('Du=',num2str(Du),aaa);
            sgtitle(sgtitleName)
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

