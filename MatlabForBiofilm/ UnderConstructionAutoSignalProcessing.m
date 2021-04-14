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
%Results={};
for j=0:0
    
    folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/FullPaper/AutoModes/Results0609/ParameterSet';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ALIFE2020CameraReady/AntiPhase';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/OldForALIFE2020/SmallField/InPhase';
    %folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ParameterSet3';
    folderName=strcat(folderName,num2str(j));
    cd(folderName);
    fileName=strcat(folderName,'/Vsum.csv');
    %fileName=strcat(folderName,'/V_para.csv');
    %fileName=strcat(folderName,'/SumV.csv');
    fileData=csvread(fileName);
    [row,col]=size(fileData);
    %%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
    if FlagSingle==0
        start=2;
        dend=row;
        b=cast((col-2)/2,'int64');
        %DataLeft=sum(fileData(start:dend,2:b),2);
        % DataRight=sum(fileData(start:dend,b:col),2);
        % DataLeft=DataLeft.';
        % DataRight=DataRight.';
        DataLeft=[];
        DataRight=[];
        DataLeft=fileData(start:dend,2);
        DataRight=fileData(start:dend,3);
        %XData=(fileData(start:dend,1)/1000000);%ALIFE2020
        XData=(fileData(start:dend,1)/10000);%Fullpaper
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
        filenameBlend=strcat(folderName,'/Expansion.png');
        saveas(gcf,filenameBlend);
        close
        %グラフの微分値を取得
        dDataLeft=diff(DataLeft);
        dDataRight=diff(DataRight);
        [row2,col2]=size(dDataLeft);
        %plot(dDataLeft)
        
        %%%%%%ローパスフィルタに通した拡大率を表示%%%%%%%
        %ローパスフィルタ の設計
        fc =0.40;
        Wn = (2.0/row2)*fc;%カットオフ
        % t = linspace(0,1,row2);
        %b = fir1(200,Wn,'low');%2D
        b = fir1(5.0,Wn,'low');%1D
        %フィルタ を適用
        FilterdDataLeft=[];
        FilterdDataRight=[];
        FilterdDataLeft=filter(b,1,dDataLeft);
        FilterdDataRight=filter(b,1,dDataRight);
        [rowFilterd,colFiltered]=size(FilterdDataLeft);
        XDataFiltered=XData(1:rowFilterd);
        
        %%%%グラフ表示
        figure;
        %%左側のバイオフィルム
        ax=plot(XDataFiltered,FilterdDataLeft,'LineWidth',LineHaba);
        %findpeaks(FilterdDataLeft,XDataFiltered)
        hold on;
        %%右側のバイオフィルム をプロット
        plot(XDataFiltered,FilterdDataRight,'-.','LineWidth',LineHaba);
        %findpeaks(FilterdDataRight,XDataFiltered)
        %{
        %%各グラフの極小点を検出し，グラフに追加
        TF1 = islocalmin(FilterdDataLeft, 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        plot(XDataFiltered(TF1),FilterdDataLeft(TF1),'r*','LineWidth',LineHaba);%極小値をグラフに追加
        TF2 = islocalmin(FilterdDataRight, 'FlatSelection', 'center','MinProminence',0.05);%極小値検出．起伏が0.05以上の点のみ抽出
        plot(XDataFiltered(TF2),FilterdDataRight(TF2),'g*','LineWidth',LineHaba);%極小値をグラフに追
        %}
    %%各グラフの極大点を検出し，グラフに追加
    XlimStart=4.0;
    XlimEnd=6.0;
    %XlimStartとXlimEndの範囲のみ計算するために，データ範囲を調整．
    interval=XDataFiltered(2)-XDataFiltered(1);
    StepStart=(XlimStart+interval-min(XDataFiltered))/interval;
    StepEnd=(XlimEnd+interval-min(XDataFiltered))/interval;
    XDataFiltered=XDataFiltered(StepStart:StepEnd);
    [pks1,TF1] = findpeaks(FilterdDataLeft(StepStart:StepEnd),'MinPeakProminence',0.05e+0);
    %[pks1,TF1] = findpeaks(FilterdDataLeft,'MinPeakProminence',0.05e+0);
    plot(XDataFiltered(TF1),pks1,'r*','LineWidth',LineHaba);%極小値をグラフに追加
    [pks2,TF2] = findpeaks(FilterdDataRight(StepStart:StepEnd),'MinPeakProminence',0.05e+0);
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
        filenameBlend=strcat(folderName,'/dExpansion.png');
        %%グラフを保存
        saveas(gcf,filenameBlend);
        close
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
        %Phi1=int64.empty(Num,0);
        %Phi1=int64.empty(Num,0);
        %sizePhi=size(Phi1);
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
        j
        PhiDiff/pi;
        abc=mean(PhiDiff)/pi
        Results(1,j+1)=abc;
        PhiDiff=[];
        %{
        if abs(Results(j+1))>1.0
            Results(j+1)=2.0-abs(Results(j+1));
        end
        %}
        
    end
    
end
figure;
Results=fliplr(reshape(abs(Results),[91 65]))
xvalues={'20','30','40','50','60'};
yvalues={'1000','100','10','1','0.1','0.01'};
h=heatmap(xvalues,yvalues,(Results))
h.xlabel('Food Source u_{in}');
h.ylabel('Signalling Capability k_{vz}');

%%%単バイオフィルム%%%%
if FlagSingle==1
    Data=fileData(1:rowEnd,2);
    XData=(fileData(1:rowEnd,1)/1000);
    f=figure;
    plot(XData,Data,'LineWidth',LineHaba);
    xlim([min(XData) max(XData)]);
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%軸の太さ
    box off%余分な右のy軸と上のx軸を削除
    %軸ラベルの制御
    ax=gca;
    %ax.XTick=[0 100 200 300 400 500 600 700]%強制的に軸の刻みを100ずつにする
    %axis([0 700 0 500])
    %xticklabels({'0','10','20','30','40','50','60','70'})%目盛りの名前付け直し
    %xticklabels(strXData);
    xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold');
    ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold');
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);%目盛りのない上側のx軸と右側のy軸を追加
    %グラフ保存
    filenameBlend=strcat(folderName,'/Expansion.png');
    saveas(gcf,filenameBlend);
    %グラフの微分値を取得
    dData=diff(Data);
    [row2,col2]=size(dData);
    %plot(dDataLeft)
    %%%%%%ローパスフィルタに通した拡大率を表示%%%%%%%
    %ローパスフィルタ の設計
    fc =10000.0;
    Wn = (2.0/row2)*fc;
    t = linspace(0,0.0001,row2);
    b = fir1(200,Wn,'low');
    %フィルタ を適用
    %dData=dData(450:925);
    FilterdDataLeft=filter(b,1,dData);
    size(FilterdDataLeft);
    %グラフ表示
    [row3,col3]=size(dData);
    %XData=XData(450:925);
    max(XData);
    figure;
    ax=plot(XData,FilterdDataLeft,'LineWidth',LineHaba);
    xlim([min(XData) max(XData)]);
    %xticklabels({'0','20','40','60','80','100'})
    %xticklabels({'70','75','80','85','90','95','100'})
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);
    xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold');
    ylabel('Expanding rate [/step]','FontSize',LabelFontSize,'FontWeight','bold');
    ax=gca;
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba);
    filenameBlend=strcat(folderName,'/dExpansion.png');
    saveas(gcf,filenameBlend);
end

%{
%周波数パワースペクトルを調べる
[row2,col2]=size(dDataLeft)
y = fft(dDataLeft,251);    % Compute DFT of x
n = 251       % number of samples
fs=20;
f = (0:n-1)*(fs/n);
Pyy = y.*conj(y)/row2;
%plot(f(1:50),Pyy(1:50))
%plot(f,Pyy(1:128))
title('Power spectral density')
xlabel('Frequency (Hz)')
%}