clear

%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=2;
DotHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
NumOfGyouInSubPlot=1;
NumOfRetsuInSubPlot=5;
SubPlotTitleFontSize=25;

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
a=[37 29 21 13 5];%%gyou5 retsu4
%a=[22 8];%%gyou2 retsu7

IndexOfParaset=[2 9 20 30 39];
IndexOfU0Value=[200 1000 2100 3100 4000];

DisSizeX=5000;
DisSizeY=250;
SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 DisSizeX DisSizeY]);
folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Now/Paper/Ok/Distance/Du/CompareNormWithDifferentTrial/Interval1000/Done';
OutputFolder='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Now/Paper/Ok/Distance/Du/CompareNormWithDifferentTrial/Interval1000';

SSSS=24;
for BossSet=SSSS:-1:1
    FolderPatterns=["/Norm","/Norm2","/Norm3"];
    
    minDist=2;
    maxDist=12;
    Sets=0;
    for ParameterSet=0:1:4
        for PatternNum=1:1:3%比較用のフォルダ名
            ThrePeak=1.0;
            PeakPro= 2.0;%interval100 のとき 3.0   後半 25.0 paraset=20-8の時:10.0 3-8の時:5.0
            %interval100の時，Bossset=1-6まではずっと3.0でok
            %{
            if BossSet <6
                PeakPro= 3.0;
            else
                if ParameterSet>=20
                    PeakPro=25.0;
                elseif ParameterSet>=8
                    PeakPro=10.0
                else
                    PeakPro=5.0
                end
            end
            %}
            IndexOfParaset(ParameterSet+1)>=20
            threshh=15
            if 8 < BossSet  && BossSet <=threshh%interval 1000の時
                PeakPro= 500.0;
                %{
                if ParameterSet==0
                    PeakPro=200.0;
                end
                %}
                if IndexOfParaset(ParameterSet+1)<20
                    PeakPro=100.0;
                end
            else
                PeakPro=100.0;
            end
            %{
            if ParameterSet>=20
                tempParaSet=ParameterSet-20;
            else
                tempParaSet=ParameterSet;
            end
            %}
            
            IndexOfSubPlot=ParameterSet+1;
            
            folderName=strcat(folderName0,FolderPatterns(1,PatternNum),'/',num2str(BossSet),'/Results','/ParameterSet',num2str(IndexOfParaset(ParameterSet+1)));
            cd(folderName);
            %/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Now/Paper/Ok/Distance/Du/CompareNormWithDifferentTrial/Interval1000/Done/Norm/24/Results/ParameterSet2
            fileName=strcat(folderName,'/SumV');
            %fileName=strcat(fileName,num2str(ParameterSet));
            fileName=strcat(fileName,'.csv')
            fileID=fopen(fileName);
            % if fileID~=-1
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
            %Distance=fileData(start:dend,4);
            XData=(fileData(start:dend,1)/10000);%Fullpaper
            dDataLeft=diff(DataLeft);%微分
            dDataTop=diff(DataTop);
            
            [row,col]=size(dDataTop);
            
            [pks1,locsLeft] = findpeaks(dDataLeft(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
            [pks3,locsTop] = findpeaks(dDataTop(start:row-1),'MinPeakHeight',ThrePeak,'MinPeakProminence',PeakPro);
            
            %%%%バイオフィルム間の位相差を求める
            Top=[];
            Left=[];
            Right=[];
            mTop=1;
            mRight=1;
            mLeft=1;%足をついた回数
            if isempty(locsTop)==0&&isempty(locsLeft)==0&&length(locsTop)>2&&length(locsLeft)>2
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
                    if i<locsTop(3)||i<locsLeft(3)% 1個目の極大値よりも前の時
                        TopAndLeft(i)=NaN;
                    end
                end
                
                [Num,col]=size(TopAndLeft);
                
                if any(TopAndLeft(:)>0.0)~=0
                    BossSet;
                    ParameterSet;
                end
                
                figure(SubPlotPhaseDiffernce);
                subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
                plot(TopAndLeft,'o','MarkerSize',3)
                ylim([0 2.0*pi]);
                hold on;
                xlim([0 row-1]);
                Du=(BossSet)*0.2;
                sgtitleName=strcat('Du=',num2str(Du));
                sgt=sgtitle(sgtitleName)
                sgt.FontSize=SubPlotTitleFontSize;
                titlename=strcat('u_{in}=',num2str(IndexOfU0Value(ParameterSet+1)));
                Tit=title(titlename);
                Tit.FontSize=SubPlotTitleFontSize-5;
                hold on;
                
            end
            %end
        end
        
        if ParameterSet==4
            figure(SubPlotPhaseDiffernce);
            fileName1=strcat(OutputFolder,'/',num2str(BossSet),'.png');
            saveas(gcf,fileName1);
            clf
        end
        
        
    end
    
end

