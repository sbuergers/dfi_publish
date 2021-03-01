function [extracted]=ExtractAllSaccades(infiles,skip)
%In command window, call f=dir('wildcards.asc')
%which returns a structure. The 'name' field of the structure contains the
%filename. If you pass the entire structure to this function it will go
%through each file in sequence to return one large 'extracted' matrix. The
%'skip' parameter controls whether trials without saccades between display
%onset and trial end are included in the output. Set to 1 if you do not
%want to include such trials.

extracted=[]; % matrix to be filled

for findex=1:length(infiles)
    
    fin=fopen(infiles(findex).name,'r');
    filestr=fscanf(fin,'%c'); % read entire file into 1 string

    % search the string for the critical messages
    sTrlIdx = strfind(filestr,'TRIALID');
    SacIdx  = strfind(filestr,'ESACC');
    eTrlIdx = strfind(filestr,'END TRIAL');
    trigIdx.iti  = strfind(filestr,'TRIGGER ');
    trigIdx.resp = strfind(filestr,'TRIGGER 10');
    strfind(filestr, 'MSG	*');
    
    AllSaccadesCounter=1;
    for j=1:length(DispIndex)
        
        %get time of display onset
        disptime='';
        index=DispIndex(j)-2;
        while(uint8(filestr(index))~=9)
            index=index-1;
        end
        index=index+1;
        while(uint8(filestr(index))~=32)
            disptime=[disptime,filestr(index)];
            index=index+1;
        end

        %Once the 'DISPLAY ON' message has been identified, work backwards
        %to find the relevant TRIALID
        IDind=length(IDIndex);
        while(IDIndex(IDind)>DispIndex(j))
            IDind=IDind-1;
        end
        %and then extract the display parameters (read until the end of line after the
        %'TRIALID' message)
        ID=[];
        index=IDIndex(IDind)+8;
        temp=1;
        while(uint8(filestr(index))~=13) %until the end of the line
            tempstr='';
            while((uint8(filestr(index))~=32)&(uint8(filestr(index))~=13))
                tempstr=[tempstr,filestr(index)];
                index=index+1;
            end
            ID(temp)=str2num(tempstr);
            if(uint8(filestr(index))~=13)
                index=index+1;
            end
            temp=temp+1;
        end        
        
        %find first saccade after display onset
        while((AllSaccadesCounter<length(AllSacIndex))&&(AllSacIndex(AllSaccadesCounter)<IDIndex(j)))%DispIndex(j)))
            AllSaccadesCounter=AllSaccadesCounter+1;
        end
        %find last saccade of that trial
        LastSaccadeofTrial=AllSaccadesCounter;
        while((LastSaccadeofTrial<=length(AllSacIndex))&&(AllSacIndex(LastSaccadeofTrial)<EndIndex(j)))
            LastSaccadeofTrial=LastSaccadeofTrial+1;
        end

        %first read the response code (number after 'TRIAL_RESULT')
        index=EndIndex(j)+13; %go to the end of the message and skip the blank
        tempstr='';
        while(uint8(filestr(index)~=13)) %end of line
            tempstr=[tempstr filestr(index)];
            index=index+1;
        end
        result=str2num(tempstr);
        
        %check whether saccades were made in between display onset and
        %trial end; if so, proceed to extract the data from the relevant
        %saccades
        if((AllSacIndex(AllSaccadesCounter)>IDIndex(j))&&(AllSacIndex(AllSaccadesCounter)<EndIndex(j)))
            
            SacNr=1;
            for i=AllSaccadesCounter:LastSaccadeofTrial-1
                index=AllSacIndex(i)+9;
                for c=1:8
                    tempstr='';
                    while((uint8(filestr(index))~=32)&(uint8(filestr(index))~=9))
                        tempstr=[tempstr,filestr(index)];
                        index=index+1;
                    end
                    if (isempty(str2num(tempstr))==0)
                        sacdata(c)=str2num(tempstr);
                    else
                        sacdata(c)=-999;
                    end
                    while((uint8(filestr(index))==32)|(uint8(filestr(index))==9))
                        index=index+1;
                    end
                end
                %write the data to the output matrix
                extracted=[extracted;findex ID SacNr sacdata(1)-str2num(disptime) sacdata(3:8) result];
                SacNr=SacNr+1;
            end
        else %if no saccades were made, write the response code to the output
            if(skip==0)
                extracted=[extracted;findex ID -999 -999 -999 -999 -999 -999 -999 -999 result];
            end
        end
        if(LastSaccadeofTrial<=length(AllSacIndex))
            AllSaccadesCounter=LastSaccadeofTrial;
        end
    end
    fclose(fin);
end