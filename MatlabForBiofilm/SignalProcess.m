clear

FlagSingle=0;
rowEnd=1000;
%%グラフの制御
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
%XlimStart=2.0;
%XlimEnd=3.1;

folderName='/Users/mikamihiroshi/Desktop/Single/Results/ForPoster';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ALIFE2020CameraReady/AntiPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/OldForALIFE2020/SmallField/InPhase';
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/oscillation_1D/DrugInstall/normalized/ParameterModifying/Coupling/Results/ParameterSet3';
%folderName=strcat(folderName,num2str(79));
%cd(folderName);
%fileName=strcat(folderName,'/Vsum.csv');
%fileName=strcat(folderName,'/V_para.csv');
fileName=strcat(folderName,'/SumV.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
%%%2つのバイオフィルムのサイズの時間推移を表示%%%%%%%%%%%
if FlagSingle==0
    start=1;
    dend=row;
    figure;
    b=cast((col-2)/2,'int64');
    DataLeft=fileData(start:dend,2);
    XData=(fileData(start:dend,1)/1000000);%ALIFE2020
    %XData=(fileData(start:dend,1)/10000);%Fullpaper
    plot(XData,DataLeft,'LineWidth',LineHaba);
    hold on;
    c=legend('Biofilm1','Biofilm2');%凡例を制御
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%軸の太さ
    XlimStart=4.0;
    XlimEnd=4.20;
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
    %グラフの微分値を取得
    dDataLeft={};
    dDataLeft=diff(DataLeft);
    [row2,col2]=size(dDataLeft);
    figure;
    %%左側のバイオフィルム
    XData=XData(1:row2)
    ax=plot(XData,dDataLeft,'LineWidth',LineHaba);
    %findpeaks(FilterdDataLeft,XDataFiltered)
    c=legend('Biofilm1');%凡例を制御
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%軸の太さ

    xlim([XlimStart XlimEnd]);
    % xlim([min(XData) inf]);
    box off%余分な右のy軸と上のx軸を削除
    %軸ラベルの追加
    ax=gca;
    xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Expansion rate ','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%メモリのない上側のx軸と右側のy軸を追加
   
    filenameBlend=strcat(folderName,'/dExpansion.png');
    %%グラフを保存
    saveas(ax,filenameBlend);

end


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

