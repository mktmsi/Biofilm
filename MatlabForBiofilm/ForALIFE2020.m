clear
FlagSingle=0;
rowEnd=1000;
%%�O���t�̐���
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
%XlimStart=2.0;
%XlimEnd=3.1;
folderName='/Users/mikamitaishi/Dropbox/Programs/Biofilm/Oscillation/1D/Working/Single/1/Results/ParameterSet0';
cd(folderName);
fileName=strcat(folderName,'/SumV.csv');
%fileName=strcat(folderName,'/V_para.csv');
%fileName=strcat(folderName,'/SumV.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
if row==282
    fileData=fileData(1:row-1,:);
end
[row,col]=size(fileData);
%%%2�̃o�C�I�t�B�����̃T�C�Y�̎��Ԑ��ڂ�\��%%%%%%%%%%%
if FlagSingle==0
    start=1;
    dend=row;
    figure;
    b=cast((col-2)/2,'int64');
    %DataLeft=sum(fileData(start:dend,2:b),2);
    % DataRight=sum(fileData(start:dend,b:col),2);
    % DataLeft=DataLeft.';
    % DataRight=DataRight.';
    DataLeft=fileData(start:dend,5);
    %DataRight=fileData(start:dend,3);
     XData=(fileData(start:dend,1)/1000000);%ALIFE2020
   % XData=(fileData(start:dend,1)/1000);%ALIFE2020
    plot(XData,DataLeft,'LineWidth',LineHaba);
    hold on;
    %plot(XData,DataRight,'LineWidth',LineHaba);
    %c=legend('Biofilm1','Biofilm2');%�}��𐧌�
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���
        XlimStart=4.90;%for alife
   % XlimEnd=5.00;%for alife
   % xlim([XlimStart XlimEnd]);% for alife
    % xlim([min(XData) inf]);
    box off%�]���ȉE��y���Ə��x�����폜
    %�����x���̒ǉ�
    ax=gca;
    xlabel('Time [\times 10^6 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%�������̂Ȃ��㑤��x���ƉE����y����ǉ�

    %�O���t�ۑ�

    filenameBlend=strcat(folderName,'/Expansion.png');
    saveas(gcf,filenameBlend);
    %�O���t�̔����l���擾
    dDataLeft=diff(DataLeft);
    %dDataRight=diff(DataRight);
    [row2,col2]=size(dDataLeft);  
    XDataFiltered=XData(1:row2)
    figure;
    %%�����̃o�C�I�t�B����
    ax=plot(XDataFiltered,dDataLeft,'LineWidth',LineHaba);
    %findpeaks(FilterdDataLeft,XDataFiltered)
    hold on;
    %%�E���̃o�C�I�t�B���� ���v���b�g
    %plot(XDataFiltered,dDataRight,'-.','LineWidth',LineHaba);
    c=legend('Biofilm1','Biofilm2');%�}��𐧌�
    c.FontSize=HanreiSize;
    c.Location='northwest';
    set(gca,'FontSize',MemoriSize);
    set(gca,'LineWidth',LineHaba);%���̑���

    %xlim([-5 XlimEnd]);
    % xlim([min(XData) inf]);
    box off%�]���ȉE��y���Ə��x�����폜
    %�����x���̒ǉ�
    ax=gca;
    XlimStart=1.5;
    XlimEnd=3.0;
    xlim([XlimStart XlimEnd]);% for alife
    %xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
    ylabel('Expansion rate ','FontSize',LabelFontSize,'FontWeight','bold')
    axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%�������̂Ȃ��㑤��x���ƉE����y����ǉ�
    
    filenameBlend=strcat(folderName,'/dExpansionNotFiltered.png');
    saveas(ax,filenameBlend);
end