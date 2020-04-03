function [tstamps, xpos, ypos] = Read_NVT_File (the_file,ts_start,ts_stop,usebounds)

    %fid = fopen(the_file,'rb','b');   %open file
    fid = fopen(the_file,'rb','a');   %open file
    
    linein = '';
    nvtdir = dir(the_file);
    filelen = nvtdir.bytes;
    
    
    %while isempty(findstr('-Resolution',linein))==1
    %    linein = fgets(fid); % read in single character
    %    fprintf('%s',linein);
    %end
    
    %fprintf('\n\n');
    HEADER_SIZE = 16384;
    header = fread(fid,HEADER_SIZE,'ubit8');      %get rid of header
    
    VT_RECORD_SIZE = 1828;
    NLX_VTREC_NUM_POINTS = 400;
    NLX_VTREC_NUM_TARGETS = 50;
    
    start_ind_set = 0;
    start_ind = 0;
    xpos_tmp = [];
    ypos_tmp = [];
    tstamps_tmp = [];
    tstamps = [];
    xpos = [];
    ypos = [];
    dataidx = 1;
    waithandle= waitbar(0.0,'Extracting data from VT file ..... ');pause(0.01);
    while ftell(fid) < filelen
        waitbar(ftell(fid)/filelen, waithandle,sprintf('Extracting data from VT file ..... '));pause(0.01);
        %fprintf('Processing VT File ..... (position %d/%d) | dataidx: %d\n',ftell(fid),filelen,dataidx);
        for i=1:1000   
            if(ftell(fid) >= filelen)
                break
            end
            swstx = fread(fid,1,'*uint16');                             % 2 bytes
            swid = fread(fid,1,'*uint16');                              % 2 bytes
            swdata_size = fread(fid,1,'*uint16');                       % 2 bytes
            qwTimeStamp = double(fread(fid,1,'uint64'));                % 8 bytes
            dwPoints = fread(fid,NLX_VTREC_NUM_POINTS,'*uint32');       % 1600 bytes
            sncrc = fread(fid,1,'*uint16');                             % 2 bytes
            dnextracted_x = double(fread(fid,1,'*uint32'));             % 4 bytes
            dnextracted_y = double(fread(fid,1,'*uint32'));             % 4 bytes
            dnextracted_angle = fread(fid,1,'*uint32');                 % 4 bytes
            dntargets = fread(fid,NLX_VTREC_NUM_TARGETS,'*uint32');     % 200 bytes
            if isempty(qwTimeStamp) == 1
                break;
            end
            if(usebounds)
                if( (qwTimeStamp >= ts_start) && (qwTimeStamp <= ts_stop) )
                    %if( (i >= ts_start) && (start_ind_set == 0))
                    %    start_ind = dataidx;
                    %    start_ind_set = 1;
                    %end
                    if(dataidx == 1)
                        dataidx = 2;
                    end
                    tstamps(dataidx) = qwTimeStamp;
                    xpos(dataidx) = dnextracted_x;
                    ypos(dataidx) = dnextracted_y;
                    dataidx = dataidx + 1;
                elseif(qwTimeStamp < ts_start)
                    tstamps(dataidx) = qwTimeStamp;
                    xpos(dataidx) = dnextracted_x;
                    ypos(dataidx) = dnextracted_y;
                elseif(qwTimeStamp > ts_stop)
                    tstamps(dataidx) = qwTimeStamp;
                    xpos(dataidx) = dnextracted_x;
                    ypos(dataidx) = dnextracted_y;
                    dataidx = dataidx + 1;
                    %fprintf('Found last point\n');
                    break
                end
            else
                tstamps(dataidx) = qwTimeStamp;
                xpos(dataidx) = dnextracted_x;
                ypos(dataidx) = dnextracted_y;
                dataidx = dataidx + 1;
            end
        
            %fprintf('swstx: %d\n',double(swstx));
            %fprintf('Timestamp: %d\n',double(qwTimeStamp));
        end    %end for loop
        if((usebounds) && (qwTimeStamp > ts_stop))
            break
        end
    end
    fclose(fid);
    %waitbar(ftell(fid)/filelen, waithandle,'Preparing VT Data ..... ');pause(0.01);
    %if(usebounds)
    %    while ( (start_ind > 1) && (xpos(start_ind) == 0))
    %        start_ind = start_ind - 1;
    %    end
   %     tstamps = tstamps_tmp(start_ind:(dataidx-1));
    %    xpos = xpos_tmp(start_ind:(dataidx-1));
    %    ypos = ypos_tmp(start_ind:(dataidx-1));
    %    fprintf('\n  -> Found %d position entries\n\n',dataidx - start_ind);
    %else
        %fprintf('\n  -> Found %d position entries\n\n',dataidx - 1);
        %end
    close(waithandle);
    
    
    
    