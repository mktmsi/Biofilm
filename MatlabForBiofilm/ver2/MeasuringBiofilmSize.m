clear
%folderName='/Users/mikamihiroshi/Document/Research/programs/Biofilm/MembranePotentialModel/OscillationModel/Oscillation_2D/Hex/Syuron/Single/1/Results/ParameterSet0/ArrayV';

cd(folderName);
%dir('*.csv')
MyFolderInfo = dir('*.csv');
NumFiles=length(MyFolderInfo)

BiofilmSize=zeros(1,700)

for i=1:700%0:NumFiles-1
    fileName=strcat(folderName,'/',num2str(i),'.csv');
    i;
    fileData=csvread(fileName);
    [row,col]=size(fileData);
    count=0;
    for j=1:row
        if fileData(50,j)>0.0
            count=count+1;
        end
    end
    BiofilmSize(i)=count;
end

plot(BiofilmSize)