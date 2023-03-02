function [data_IMG, scales_YX, FrameRate] = buildBMode_matOnly(filenameRAW)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initiate the parameters for reading in the Bmode data
    % Set up file names
    fname = [filenameRAW '.bmode'];
    fnameXml = [filenameRAW '.xml'];

    % Parse the XML parameter file - DO NOT CHANGE
    param = VsiParseXml(fnameXml, '.bmode');
    BmodeNumSamples = param.BmodeNumSamples;
    BmodeNumLines = param.BmodeNumLines;
    BmodeDepthOffset = param.BmodeDepthOffset; %mm
    BmodeDepth =  param.BmodeDepth; %mm
    BmodeWidth =  param.BmodeWidth; %mm
    scales_YX = [(BmodeDepth-BmodeDepthOffset)/BmodeNumSamples, BmodeWidth/BmodeNumLines];

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This is to strip the header data in the files - DO NOT CHANGE
    size_header = 1;   %  4 bytes (# of bytes per position)
    file_header = 40;  % 40 bytes
    line_header = 0;   %  4 bytes
    frame_header = 56; % 56 bytes

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Iterate through each of the files to build the 2D DataSet
    fid = fopen(fname,'r');
    file_header_info = fread(fid, 10, 'uint');
    num_frames = file_header_info(2);
    data_IMG = zeros(BmodeNumSamples, BmodeNumLines, num_frames, 'uint8');
    timeStamps_ticks = zeros(num_frames,1);
    timeStamps_ms    = zeros(num_frames,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % For each frame in the image
    w1 = waitbar(0,'Compiling Data');
    for j = 1:num_frames
        fid = fopen(fname,'r');

        header =  file_header  + ... % since starting from new, add in the 40 byte file header
                 (frame_header + ... % size of frame header
                  size_header*BmodeNumSamples*BmodeNumLines + ... % size of image
                  BmodeNumLines*line_header) * (j-1); % size of line headers (to n-1 frames)
        
        fseek(fid, header, -1); % move to start of current frame (pre-header)
        timeStamps_ticks(j) = fread(fid,1,'uint');
        timeStamps_ms(j)    = fread(fid,1,'float64');
        fseek(fid,44,0);
        
        for i = 1:BmodeNumLines
            fseek(fid,line_header,0);
            data_IMG(:,BmodeNumLines-(i-1),j) = fread(fid, BmodeNumSamples, 'uint8');
        end

        fclose(fid);
        waitbar(j/num_frames,w1);
    end
    close(w1);
    
    FrameRate = 1/((timeStamps_ms(2)-timeStamps_ms(1))/1000);
end