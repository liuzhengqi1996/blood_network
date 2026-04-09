function [faceMx,ptCoordMx,dia,BC,np,nf,nt]=caseReaderMJ(filename)
myfileName = strcat(filename,'.fMx');
fileInfo = dir(myfileName);
if fileInfo.bytes == 0
    faceMx=[];
else
    faceMx = load(myfileName);
end
myfileName = strcat(filename,'.pMx');
fileInfo = dir(myfileName);
if fileInfo.bytes == 0
    ptCoordMx=[];
else
    ptCoordMx = load(myfileName);
end
myfileName = strcat(filename,'.dia');
fileInfo = dir(myfileName);
if fileInfo.bytes == 0
    fprintf('%s\n','dia file is EMPTY');
    dia=[];
else    
    dia = load(myfileName);
    fprintf('%s\n','dia file is loaded');
end
myfileName = strcat(filename,'.BC');
fileInfo = dir(myfileName); 
if fileInfo.bytes == 0
    BC=[];
else   
    BC = load(myfileName);
end
np= length(ptCoordMx); nf= length(faceMx(:,2)); nt= np+nf;
end
