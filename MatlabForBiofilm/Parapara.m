clear
count1=0;
count2=0;

folderName='/Users/mikamihiroshi/Dropbox/Programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/FullPaper/Revised/0623/AutoModes0623_revised1/1/Results/ParameterSet2';
MovieFolder='/Users/mikamihiroshi/Dropbox/Programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/FullPaper/Revised/0623/AutoModes0623_revised1/1/Results/Numbers';
cd(folderName);
fileName=strcat(folderName,'/Vsum.csv');
fileData=csvread(fileName);
[row,col]=size(fileData);
start=2;
dend=row;
DataLeft=fileData(start:dend,2);
XData=(fileData(start:dend,1)/1000000);%ALIFE2020
for j=0:1:8735;
    
    

    if j==0
        j
        count1
        count2
    else
        j
        if mod(count1 , 90) == 0 && count1 ~= 0
            count1 = 0
            count2=count2+1
        else
            count1=count1+1
            count2
        end
    end
    % //パラメータの更新
    % hiku = 50.0 / (pow(10.0, (double)(count1 / 18)));
    10.0^fix(cast(count1 / 18,'double'))
    hiku = 50.0 / (10.0^(fix(cast(count1 / 18,'double'))))
    %k_vz = (1000.0 / (pow(10.0, (count1 / 18)))) - hiku * ((double)(count1 % 18));
    k_vz = (1000.0 / ((10.0^ fix(count1 / 18))) - hiku * cast(mod(count1 , 18),'double'))
    u0 = (1000.0 - 10.0 * cast(count2,'double')) %20-60くらい？10刻み
    
    figure
    XData=(fileData(start:dend,1)/10000);%Fullpaper
    plot(XData,DataLeft);
    title(['Setcount=',num2str(j),' u_0=',num2str(u0),' k_{vz}=',num2str(k_vz)],'FontSize',20)
    filenameBlend=strcat(MovieFolder,'/dExpansion',num2str(j),'.png');
    saveas(gcf,filenameBlend);
    close
    
    
    
end



