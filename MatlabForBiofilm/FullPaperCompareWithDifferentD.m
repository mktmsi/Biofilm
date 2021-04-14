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
NumOfGyouInSubPlot=5;
NumOfRetsuInSubPlot=4;

IndexOfSubPlot=0;
%a=[91 71 51 31 11];
a=[131 111 91 71 51 31 11];%gyou7 retsu10
%a=[71 51 31 11];%gyou4 retsu10
a=[37 29 21 13 5];%%gyou5 retsu4
%a=[22 8];%%gyou2 retsu7

DisSize=5000;
SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
folderName0='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Now/Paper/Ok/Distance/Du';
OutputFolder='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Now/Paper/Ok/Distance/Du/CompareNormalWithDifferentDu';

minDist=2;
maxDist=12;

Sets=39;
for ParameterSet=Sets:-1:0
    
    SSSS=22;
    
    for BossSet=1:1:SSSS
        
        PhaseDifferenceJudge=zeros(NumOfGyouInSubPlot,NumOfRetsuInSubPlot);
        FolderPatterns=["/Norm","/NormTake2","/NormTake3"];
        
        for PatternNum=1:1:1%比較用のフォルダ名
            
            ThrePeak=1.0;
            PeakPro= 100.0;
            if BossSet >8 && BossSet <16
                PeakPro= 500.0;
            end
            
            if ParameterSet>=20
                tempParaSet=ParameterSet-20;
            else
                tempParaSet=ParameterSet;
            end
            q=fix(tempParaSet/NumOfRetsuInSubPlot);%q:商  r:余り
            IndexOfSubPlot=a(1,q+1)-(NumOfRetsuInSubPlot*NumOfGyouInSubPlot-tempParaSet);
            
            folderName=strcat(folderName0,FolderPatterns(1,PatternNum),'/',num2str(BossSet),'/Results','/ParameterSet',num2str(ParameterSet));
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
            Top=[];
            Right=[];
            Left=[];
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
                
                Bluee=BossSet/25;
                hsv=[0.3 1.0 Bluee];
                rgb=hsv2rgb(hsv);
                figure(SubPlotPhaseDiffernce);
                subplot(NumOfGyouInSubPlot,NumOfRetsuInSubPlot,IndexOfSubPlot);
                plot(TopAndLeft,'o','MarkerSize',1.0,'MarkerEdgeColor',rgb)
                ylim([0 2.0*pi]);
                hold on;
                xlim([0 row-1]);
                %set(gca,'defaultAxesColorOrder',[[0 0 0]; [cast(35/255,'double') cast(140/255,'double') cast(42/255,'double')]]);
                titlename=strcat('i=',num2str(ParameterSet));
                title(titlename);
                hold on;
                
            end
        end
        
        if ParameterSet==20||ParameterSet==0
            if ParameterSet==20
                aaa='later2';
            else
                aaa='faster';
            end
            
            figure(SubPlotPhaseDiffernce);
            fileName1=strcat(OutputFolder,'/',num2str(BossSet),aaa,'.png');
            %fileName1=strcat(folderName0,'/SubplotPhaseDifferences',aaa,'.png')
            saveas(gcf,fileName1);
            %pause(1)
            clf
        end
        
        if ParameterSet==0
            figure(SubPlotPhaseDiffernce);
            close
            SubPlotPhaseDiffernce=figure('Name','Phase Differnce','Position',[0 0 5000 5000]);
        end
    end
end

