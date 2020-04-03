%CleanVTData.m
% reads in raw Neuralynx VT files between a specified starting timestamp and
% an ending timestamp.  Get these timestamps after saving cluster data on familiar
% novel runs.  The program then runs Piotr's fixArtifacts.m program to 
% clean up noisy VT data and prints out timestamp, x and y positions to ascii
% file
%
% 25 April 2006
% Victoria Booth

global tstamps
global xpos
global ypos 
global rawVTdata
global cleanVTdata

working_dir = pwd;
DATADIR = 'C:\Sleepscorer04092008\';

%read in VT file
[nvt_file, nvt_path] = uigetfile({'*.Nvt',...
            'Neuralynx Video File (*.Nvt)'},'Select a Corresponding Video Position File');
vt_file = fullfile(nvt_path, nvt_file);

    
% set up variable for starting and ending timestamps
prompt = {'Enter starting timestamp (sec plus 6 zeros):'};
dlg_title = 'Starting timestamp';
prompt2 = {'Enter ending timestamp (sec plus 6 zeros):'};
dlg_title2 = 'Ending timestamp';
num_lines = 1;
def = {'1'};

%prompt for starting and ending timestamps
tempunitnum = inputdlg(prompt,dlg_title,num_lines,def);
first_vt_tstamp=str2double(tempunitnum{1});

tempunitnum2 = inputdlg(prompt2,dlg_title2,num_lines,def);
last_vt_tstamp=str2double(tempunitnum2{1});

if last_vt_tstamp <= first_vt_tstamp, error('Bad time stamps'); end

% extract VT data
cd(working_dir);
%first_vt_tstamp = 0;
%last_vt_tstamp = 0;
usebounds = 1; % = 0 reads whole file, = 1 reads between first_vt_tstamp and last_vt_tstamp
[tstamps xpos ypos] = Read_NVT_File(vt_file,first_vt_tstamp,last_vt_tstamp,usebounds);

% call Piotr's fixArtifacts function
rawVTdata(:,1) = tstamps'; 
rawVTdata(:,2) = xpos'; 
rawVTdata(:,3) = ypos';

%parameters for interpolation
% a=speed tolerance (10 by default)
% b=minimum length of correct data portions in steps (35 by default)
a=5;
b=2;

cleanVTdata = fixArtifacts(rawVTdata,a,b);

% write clean VT data to file
[filename, pathname] = uiputfile({'*.csv',...
      'Ascii VT file (*.csv)'},sprintf('Save VT data as ascii'));

cd(working_dir);
output_file_name= fullfile(pathname, filename);
        
OUTPUT_FILE = fopen(output_file_name,'w');

for i=1:length(cleanVTdata)
    fprintf(OUTPUT_FILE,'%12.7f, %12.7f, %12.7f\n',cleanVTdata(i,1),cleanVTdata(i,2),cleanVTdata(i,3));
end

fclose(OUTPUT_FILE);