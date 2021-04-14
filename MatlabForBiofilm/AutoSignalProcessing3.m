clear

FlagSingle=0;
rowEnd=1000;
%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
for j=0:1:1
    
    folderName='/Users/mikamitaishi/Dropbox/Programs/Biofilm/MembranePotentialModel/Oscillation2D/Hex/FullPaper/Now/vzBig/1/Results/ParameterSet';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ALIFE2020CameraReady/AntiPhase';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/OldForALIFE2020/SmallField/InPhase';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ParameterSet3';
    folderName=strcat(folderName,num2str(j));
    cd(folderName);
    fileName=strcat(folderName,'/Cellcounts');
    fileName=strcat(fileName,num2str(j));
    fileName=strcat(fileName,'.csv');
    %fileName=strcat(folderName,'/V_para.csv');
    %fileName=strcat(folderName,'/SumV.csv');
    fileData=csvread(fileName);
    [row,col]=size(fileData);
    %%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
    if FlagSingle==0
        start=2;
        dend=row;
        b=cast((col-2)/2,'int64');
        DataLeft=[];
        DataRight=[];
        DataLeft=fileData(start:dend,2);
        DataRight=fileData(start:dend,3);
        %XData=(fileData(start:dend,1)/1000000);%ALIFE2020
        XData=(fileData(start:dend,1)/10000);%Fullpaper
        figure;
        plot(XData,DataLeft,'LineWidth',LineHaba);
        hold on;
        plot(XData,DataRight,'LineWidth',LineHaba);
        c=legend('Biofilm1','Biofilm2');%凡例を制御
        c.FontSize=HanreiSize;
        c.Location='northwest';
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%軸の太さ
        XlimStart=min(XData);
        XlimEnd=max(XData);
        xlim([XlimStart XlimEnd]);
        % xlim([min(XData) inf]);
        box off%余分な右のy軸と上のx軸を削除
        %軸ラベルの追加
        
        ax=gca;
        %ax.XTick=[0 200 400 600 800 1000];%強制的に軸の刻みを200ずつにする
        %axis([0 1000 0 inf])
        %xticklabels({'0','20','40','60','80','100'})
        xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
        ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold')
        axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%メモリのない上側のx軸と右側のy軸を追加
        %グラフ保存
        filenameBlend=strcat(folderName,'/Expansion',num2str(j),'.png');
        saveas(gcf,filenameBlend);
        close
        
        %グラフの微分値を取得
        dDataLeft=diff(DataLeft);
        dDataRight=diff(DataRight);
        [row2,col2]=size(dDataLeft);
        %plot(dDataLeft)
        
        figure;
        XdData=XData(1:row2);
        XlimStart=min(XdData);
        XlimEnd=max(XdData);
        ax=plot(XdData,dDataLeft,'LineWidth',LineHaba);
        %findpeaks(FilterdDataLeft,XDataFiltered)
        hold on;
        %%右側のバイオフィルム をプロット
        plot(XdData,dDataRight,'-.','LineWidth',LineHaba);

        %XlimStartとXlimEndの範囲のみ計算するために，データ範囲を調整．
        interval= XdData(2)- XdData(1);
        StepStart=cast((XlimStart+interval-min(XdData))/interval,'int64');
        StepEnd=cast((XlimEnd+interval-min(XdData))/interval,'int64');
        XdData=XdData(StepStart:StepEnd);
        
        %極小値の値とそのインデックスを取得
        TF1 = islocalmin(dDataLeft(StepStart:StepEnd), 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        TF2 = islocalmin(dDataRight(StepStart:StepEnd), 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        plot(XdData(TF1),dDataLeft(TF1),'r*','LineWidth',LineHaba);%極小値をグラフに追加
        plot(XdData(TF2),dDataRight(TF2),'g*','LineWidth',LineHaba);%極小値をグラフに追
        A=find(dDataLeft(TF1) <= 0.0)%極小値があるインデックスのうち最小のものを見つける
        B=find(dDataRight(TF2) <= 0.0)
        min(A)
        if min(A)>min(B)
            area=min(A):StepEnd
            XdData=XdData();
            dDataLeft(StepStart:StepEnd)
        else
            XdData=XdData(min(B):StepEnd);
        end
        
        
        %極大値の値とそのインデックスを取得
        [pks1,locs1] = findpeaks(dDataLeft(StepStart:StepEnd),'MinPeakProminence',0.05);
        %[pks1,TF1] = findpeaks(FilterdDataLeft,'MinPeakProminence',0.05e+0);
        [pks2,locs2] = findpeaks(dDataRight(StepStart:StepEnd),'MinPeakProminence',0.05);
        % [pks2,TF2] = findpeaks(FilterdDataRight,'MinPeakProminence',0.05e+0);
        %{
        if TF1(1)>TF2(2)
            dDataLeft=dDataLeft(TF1(1)+1:StepEnd);
            dDataRight=dDataRight(TF1(1)+1:StepEnd);
        else
            dDataLeft=dDataLeft(TF2(1)+1:StepEnd);
            dDataRight=dDataRight(TF2(1)+1:StepEnd);   
        end
        %}
        %plot(XdData(locs2),pks2,'g*','LineWidth',LineHaba);%極大値をグラフに追
        %plot(XdData(locs1),pks1,'r*','LineWidth',LineHaba);%極大値をグラフに追加
        
        
        %グラフの体裁
        xlim([XlimStart XlimEnd]);
        c=legend('Biofilm1','Biofilm2');
        c.FontSize=HanreiSize;
        c.Location='northwest';
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);
        xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize);
        ylabel('Expanding rate [/step]','FontSize',LabelFontSize);
        ax=gca;
        axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);
        filenameBlend=strcat(folderName,'/dExpansion',num2str(j),'.png');
        %%グラフを保存
        saveas(gcf,filenameBlend);
        close;
        
        figure;
        tempnorm = dDataRight - mean(dDataRight);
        fs = length(dDataRight);
        t = (0:length(dDataRight) - 1)/fs;
        [pxx,f] = periodogram(tempnorm,[],[],fs);
        plot(f,pxx)
        pxx=find(pxx == max(pxx))
        f=f(pxx)%バイオフィルムが振動している周波数（のはず）
        L=fs/f;%バイオフィルムの振動周期（のはず）
        
        %{
        figure;
        %FilterdDataLeft=lowpass(dDataLeft,0.3,'ImpulseResponse','iir','Steepness',0.96);
        FilterdDataRight=lowpass(dDataRight,f,fs,'ImpulseResponse','iir','Steepness',0.99);
        plot(FilterdDataRight,'-.','LineWidth',LineHaba);
        filenameBlend=strcat(folderName,'/Filtered',num2str(j),'.png');
        %%グラフを保存
        saveas(gcf,filenameBlend);
        %}
        



        
        %%%%%%ローパスフィルタに通した拡大率を表示%%%%%%%
        %ローパスフィルタ の設計
        %{
        fc =01.0;
        Wn = 0.05;%(2.0/row2)*fc;%カットオフ
        % t = linspace(0,1,row2);
        %b = fir1(200,Wn,'low');%2D
        b = fir1(30000,Wn,'low');%1D
        %フィルタ を適用
        FilterdDataLeft=[];
        FilterdDataRight=[];
        FilterdDataLeft=filter(b,1,dDataLeft);
        FilterdDataRight=filter(b,1,dDataRight);
        [rowFilterd,colFiltered]=size(FilterdDataLeft);
        XDataFiltered=XData(1:rowFilterd);
        %}
        
        %%%%グラフ表示
        %{
        figure;
        %%左側のバイオフィルム
        ax=plot(XDataFiltered,FilterdDataLeft,'LineWidth',LineHaba);
        %findpeaks(FilterdDataLeft,XDataFiltered)
        hold on;
        %%右側のバイオフィルム をプロット
        plot(XDataFiltered,FilterdDataRight,'-.','LineWidth',LineHaba);
        
        %findpeaks(FilterdDataRight,XDataFiltered)
        
        %%各グラフの極小点を検出し，グラフに追加
        TF1 = islocalmin(FilterdDataLeft, 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        plot(XDataFiltered(TF1),FilterdDataLeft(TF1),'r*','LineWidth',LineHaba);%極小値をグラフに追加
        TF2 = islocalmin(FilterdDataRight, 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        plot(XDataFiltered(TF2),FilterdDataRight(TF2),'g*','LineWidth',LineHaba);%極小値をグラフに追
       
        %%各グラフの極大点を検出し，グラフに追加
        XlimStart=min(XDataFiltered);
        XlimEnd=max(XDataFiltered);
         
        %XlimStartとXlimEndの範囲のみ計算するために，データ範囲を調整．
        interval=XDataFiltered(2)-XDataFiltered(1);
        StepStart=cast((XlimStart+interval-min(XDataFiltered))/interval,'int64');
        StepEnd=cast((XlimEnd+interval-min(XDataFiltered))/interval,'int64');
        XDataFiltered=XDataFiltered(StepStart:StepEnd);
        [pks1,TF1] = findpeaks(FilterdDataLeft(StepStart:StepEnd),'MinPeakProminence',0.05);
        %[pks1,TF1] = findpeaks(FilterdDataLeft,'MinPeakProminence',0.05e+0);
        plot(XDataFiltered(TF1),pks1,'r*','LineWidth',LineHaba);%極小値をグラフに追加
        [pks2,TF2] = findpeaks(FilterdDataRight(StepStart:StepEnd),'MinPeakProminence',0.05);
        % [pks2,TF2] = findpeaks(FilterdDataRight,'MinPeakProminence',0.05e+0);
        plot(XDataFiltered(TF2),pks2,'g*','LineWidth',LineHaba);%極小値をグラフに追
        xlim([XlimStart XlimEnd]);
        
        %xlim([min(XData) inf]);
        %xl = xticklabels
        
        box off
        %xticklabels({'0','20','40','60','80','100'})
        %xticklabels({'70','75','80','85','90','95','100'})
        %%グラフの体裁を整える
        
        c=legend('Biofilm1','Biofilm2');
        c.FontSize=HanreiSize;
        c.Location='northwest';
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);
        xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize);
        ylabel('Expanding rate [/step]','FontSize',LabelFontSize);
        ax=gca;
        axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);
        filenameBlend=strcat(folderName,'/Filtered',num2str(j),'.png');
        %%グラフを保存
        saveas(gcf,filenameBlend);
        
        close
        %}
        flagLocalMin=0;
        %%%%バイオフィルム間の位相差を求める
        k1={};
        k2={};
        if flagLocalMin%極小値に基づいて位相差を決める時
            k1=find(TF1)%極小点のインデックスを検出
            k2=find(TF2)%極小点のインデックスを検出
        else
            k1=TF1;
            k2=TF2;
        end
        SizeK1=size(k1)
        SizeK2=size(k2)
        Num={};
        if SizeK1(1)>SizeK2(1)%極値の数が少ない方をNumに代入
            Num=SizeK2(1);
        else
            Num=SizeK1(1);
        end
        
        if Num<2%極値がないとき
            PhiDiff=[]
        else
            if k1(1)<k2(1)
                later=k2;
                faster=k1;
            else
                later=k1;
                faster=k2;
            end
            NumLater=1;
            NumFaster=1;%各グラフのローカルインデックス
            %飛ばしがないか確認（グラフの最初）
            if later(1)>faster(2)%飛ばしがあるとき
                for i=3:Num-1
                    if later(1)<faster(i)%飛ばしでなくなればbreak
                        tmp=later;
                        later=faster;
                        faster=later;
                        NumLater=i;
                        NumFaster=1;
                        break;
                    end
                end
            end
            for i=1:Num-1%このIはグローバルインデックス
                if NumLater>Num-1 || NumFaster>Num-1%このパラメータセットの終了条件
                    break;
                end
                %飛ばしがないか確認（グラフのとちゅう）
                if later(NumLater)>faster(NumFaster+1)%飛ばしがあるとき
                    for k=2:Num-1
                        if NumLater+k>Num-1 || NumFaster+k>Num-1%このパラメータセットの終了条件
                            break;
                        end
                        if later(NumLater)<faster(NumFaster+k)%飛ばしでなくなったら
                            tmp=later;
                            later=faster;
                            faster=later;
                            tmpNum=NumLater;
                            NumLater=NumFaster+k;
                            NumFaster=NumLater;
                            break;
                        end
                    end
                end
                %分母を計算
                BunboFaster={};
                BunboLater={};
                BunboFaster=faster(NumFaster+1)-faster(NumFaster);
                BunboLater=later(NumLater+1)-later(NumLater);
                %%位相差計算
                PhiDiff(1,i)=2*pi*(((faster(NumFaster+1)-faster(NumFaster))/BunboFaster)-((faster(NumLater+1)-later(NumLater))/BunboLater));
                if abs(PhiDiff(i))>pi
                    PhiDiff(1,i)=2.0*pi-abs(PhiDiff(1,i));
                end
                NumFaster=NumFaster+1;
                NumLater=NumLater+1;
            end
        end
    end
    %{
        if flagLocalMin%極小値に基づいて位相差を決める時
            k1=find(TF1)%極小点のインデックスを検出
            k2=find(TF2)%極小点のインデックスを検出
        else
            k1=TF1;
            k2=TF2;
        end
        SizeK1=size(k1)
        SizeK2=size(k2)
        
        if SizeK1(1)>SizeK2(1)%小さい方を代入
            Num=SizeK2(1);
        else
            Num=SizeK1(1);
        end
        if Num==0
            PhiDiff=[]
        else
            for i=1:Num-1
                Bunbo1(i)=k1(i+1)-k1(i);
                Bunbo2(i)=k2(i+1)-k2(i);
                if k1(i)<k2(i)
                    k=k1;
                else
                    k=k2;
                end
                PhiDiff(1,i)=2*pi*(((k(i+1)-k1(i))/Bunbo1(i))-((k(i+1)-k2(i))/Bunbo2(i)));
                %PhiDiff(1,i)=2*pi*(((k(i+1)-k1(i))/Bunbo1(i))-((k(i+1)-k2(i))/Bunbo2(i)));
                if abs(PhiDiff(i))>pi
                    PhiDiff(1,i)=2.0*pi-abs(PhiDiff(i));
                end
            end
        end
    %}
    j
    PhiDiff/pi
    abc=mean(PhiDiff)/pi
    Results(1,j+1)=abc;
    PhiDiff=[];
    %{
        if abs(Results(j+1))>1.0
            Results(j+1)=2.0-abs(Results(j+1));
        end
    %}
    
end



figure;
Results=fliplr(reshape(abs(Results),[91 96]))
%xvalues={'20','30','40','50','60'};
%yvalues={'1000','100','10','1','0.1','0.01'};
h=heatmap((Results))
h.xlabel('Food Source u_{in}');
h.ylabel('Signalling Capability k_{vz}');

