clear

FlagSingle=0;
rowEnd=1000;
%%�O���t�̐���
LabelFontSize=30;
MemoriSize=25;
HanreiSize=25;
LineHaba=5;
XlimStart=2.0;
XlimEnd=3.1;
%Results={};
for j=5272:-1:5254
    j
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
    %%%2�̃o�C�I�t�B�����̃T�C�Y�̎��Ԑ��ڂ�\��%%%%%%%%%%%
    if FlagSingle==0
        start=2;
        dend=row;
        figure;
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
        c=legend('Biofilm1','Biofilm2');%�}��𐧌�
        c.FontSize=HanreiSize;
        c.Location='northwest';
        set(gca,'FontSize',MemoriSize);
        set(gca,'LineWidth',LineHaba);%���̑���
        XlimStart=min(XData);
        XlimEnd=max(XData)/3;
        xlim([XlimStart XlimEnd]);
        % xlim([min(XData) inf]);
        box off%�]���ȉE��y���Ə��x�����폜
        %�����x���̒ǉ�
        ax=gca;
        %ax.XTick=[0 200 400 600 800 1000];%�����I�Ɏ��̍��݂�200���ɂ���
        %axis([0 1000 0 inf])
        %xticklabels({'0','20','40','60','80','100'})
        xlabel('Time [\times 10^3 steps]','FontSize',LabelFontSize,'FontWeight','bold')
        ylabel('Biofilm size','FontSize',LabelFontSize,'FontWeight','bold')
        axes('position',ax.Position,'box','on','ytick',[],'xtick',[],'color','none','LineWidth',LineHaba)%�������̂Ȃ��㑤��x���ƉE����y����ǉ�
        %�O���t�ۑ�
        filenameBlend=strcat(folderName,'/ExpansionZenhan2.png');
        saveas(gcf,filenameBlend);
        close

        
       
    end
    
end
