% Y = fixArtifacts(R, A, B)
%
% resamples data with 16.715msec sampling intervals and fix artifacts
% within recording R.
%
% R is a 3 column matrix ; 1'st column - TIME value ; 2'nd - X coordinates ; 
% 3'rd - Y coordinates.
%
% Y is 3 column matrix with the fixed data. Format is identical as in R.
%
% Plot of original and fixed data will be shown.
%
% A is a parameter defining speed tolerance (10 by default)
% B defines minimum length of the correct data portions in steps (35 by default)
%
% Exaple:
% Fam1_Fixed = fixArtifacts(Fam1_VT)
%
% by Piotr Jablonski (c) 2006 ( piotr@umich.edu )
function m = fixArtifacts(data, a,b)

%if nargin==1, a=10; b=15; end
if nargin==1, a=10; b=35; end     % original values
if nargin==2 | nargin>3, error('You should enter first one or all 3 arguments.'); end

working_dir = pwd;
DATADIR = 'C:\';

t=data(:,1);%czas krok 16.715 msec
x=data(:,2);%pozycje x
y=data(:,3);%pozycje y

figure(1)
plot(t,x,'r.')%,'markersize',1)
hold on
plot(t,y,'g.')%,'markersize',1)
figure(2)
plot(x,y,'.k')
hold on

%resampling z stalym krokiem
x = interp1(t,x,t(1):16715:t(end));
y = interp1(t,y,t(1):16715:t(end));
t = t(1):16715:t(end);

x(x==640)=0;
x(y==640)=0;
y(y==640)=0;
y(x==640)=0;

x(y==0)=0;
y(x==0)=0;

%szybkosc 2D w czasie
lag = sqrt( (x(2:end)-x(1:end-1)).^2 + (y(2:end)-y(1:end-1)).^2 );

% hist(lag,10000)
tresh=a;%powyzej tej predkosci - szum

lag(lag>tresh)=0;
szum = lag==0;

szum(2:end+1)=szum;
szum(1)=0;
% plot(t,szum*320,'.')

x(szum)=0;
y(szum)=0;
l=zeros(size(x));

%tylko dlugie zostaja
p=0;
for i=1:length(x)
    if x(i)>0
        p=p+1;
        l(i)=p;
    else
        p=0;
    end
end
p=0;
for i=length(x):-1:1
    if x(i)>0
        p=p+1;
        l(i)=l(i)+p;
    else
        p=0;
    end
end
l=l-1;
%dluzsze niz 35 krokow zostaja
x(l<b)=0;
y(l<b)=0;

%znajdowanie odcinkow z zerami
p=0;%numeruje artefakty
on=false; %wlaczone gdy jestesmy wewnatrz artefaktu
%beg i fin pamietaja poczatki i konce artefaktow
for i=1:length(x)
    if x(i)==0 & ~on
        p=p+1;
        on=true;
        beg(p)=i-1;
    elseif x(i)~=0 & on
        fin(p)=i;
        on=false;
    end
end
%jakby sie konczylo artefaktem
l=length(beg)-length(fin);
if l>0
    fin(end+1)=length(x);
end
if beg(1)==0
    beg(1)=1;
end

%interpolacja
for i = 1:length(beg)
    tB=beg(i);
    tF=fin(i);
    vBx=x(tB);
    vFx=x(tF);
    vBy=y(tB);
    vFy=y(tF);
    x(tB+1:tF-1) = interp1( [tB tF], [vBx vFx], (tB+1):(tF-1) );
    y(tB+1:tF-1) = interp1( [tB tF], [vBy vFy], (tB+1):(tF-1) );
end
figure(1)
plot(t,x,'.k')
plot(t,y,'.b')
xlabel('time[microsec]')
ylabel('position')
legend('original X','original Y', 'fixed X', 'fixed Y')

figure(2)
plot(x,y,'.r')
legend('original', 'fixed')

% write output to .csv files
%cd(DATADIR);
%[filename, pathname] = uiputfile({'*.csv',...
%      'Cleaned VT data (*.csv)'},sprintf('Save Cleaned VT data'));
%
%cd(working_dir);
%output_file_name= fullfile(pathname, filename);
%        
%OUTPUT_FILE = fopen(output_file_name,'w');
%
%for i=1:length(t)
%    fprintf(OUTPUT_FILE,'%12.7f, %12.7f, %12.7f\n',t(i),x(i),y(i));
%end
%fclose(OUTPUT_FILE);

m(:,1)=t';
m(:,2)=x';
m(:,3)=y';