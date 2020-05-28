function []= ColorMapping(dem_name,a,varargin)

%Test DEM:
%dem_name='ap14_dtm_tmerc.cub';

dem=cuberead(dem_name);

if a < 0
    fprintf('\nERROR: a = alpha value. It must be between 0 (exclusive) and 90 (exclusive)\n');
    return;
elseif a == 0
    a=0.000001;
elseif a > 90
    fprintf('\nERROR: a = alpha value. It must be betweeen 0 and 90\n');
    return;
elseif a == 90
    a=89.99999;
end

%alpha values:
%a=89.9999999999; %angle between 0 and 90.
%a=45
%a=1
%a=45; %input alpha value

% Import Color Map
if sum(ismember(varargin,{'lut','LUT'}))>0
    fprintf('\nImporting LUT file\n');
    %TODO: Read in LUT file.
elseif sum(ismember(varargin,{'m','M'}))>0
    fprintf('\nImporting matlab colormap\n');
    %TODO: assign map from colormap file.
elseif sum(ismember(varargin,{'viridis'}))>0
    % Viridis (17)
    norm_map=viridis(17);
    %norm_map=[0.267004, 0.004874, 0.329415;0.282327, 0.094955, 0.417331;0.278826, 0.175490, 0.483397;0.258965, 0.251537, 0.524736;0.229739, 0.322361, 0.545706;0.199430, 0.387607, 0.554642;0.172719, 0.448791, 0.557885;0.149039, 0.508051, 0.557250;0.127568, 0.566949, 0.550556;0.120638, 0.625828, 0.533488;0.157851, 0.683765, 0.501686;0.246070, 0.738910, 0.452024;0.369214, 0.788888, 0.382914;0.515992, 0.831158, 0.294279;0.678489, 0.863742, 0.189503;0.845561, 0.887322, 0.099702;0.993248, 0.906157, 0.143936];
    map=round(norm_map*255);
elseif sum(ismember(varargin,{'inferno'}))>0
    %Inferno (17)
    norm_map=[0.0305,0.0221,0.1171;0.1020,0.0467,0.2492;0.1961,0.0386,0.3665;0.2915,0.0458,0.4188;0.3814,0.0772,0.4328;0.4703,0.1098,0.4286;0.5593,0.1412,0.4102;0.6476,0.1745,0.3777;0.7328,0.2143,0.3321;0.8111,0.2659,0.2756;0.8782,0.3324,0.2120;0.9307,0.4135,0.1437;0.9666,0.5064,0.0707;0.9852,0.6076,0.0241;0.9861,0.7144,0.1053;0.9690,0.8242,0.2437;0.9465,0.9289,0.4379];
    map=round(norm_map*255);
elseif sum(ismember(varargin,{'turbo'}))>0
    %Turbo (17)
    norm_map=turbomap(17);
    %norm_map=[0.2432,0.2225,0.5755;0.2729,0.3755,0.8341;0.2740,0.5191,0.9780;0.2131,0.6599,0.9790;0.1195,0.7922,0.8618;0.1021,0.8941,0.7180;0.2338,0.9604,0.5569;0.4466,0.9958,0.3708;0.6384,0.9910,0.2365;0.7835,0.9370,0.2033;0.9064,0.8429,0.2220;0.9810,0.7271,0.2208;0.9955,0.5800,0.1660;0.9623,0.4122,0.0936;0.8920,0.2684,0.0395;0.7874,0.1607,0.0133;0.6453,0.0752,0.0040];
    map=round(norm_map*255);
elseif sum(ismember(varargin,{'plasma'}))>0
    norm_map=plasma(17);
    map=round(norm_map*255);
else
    map=[30,83,151;69,117,180;97,144,195;125,172,209;154,197,223;185,217,234;215,238,245;233,247,233;255,255,191;255,244,174;255,233,157;254,215,138;253,183,116;253,150,95;243,116,76;229,82,57;215,48,39];
    norm_map=map/255;
end



%% Get Dimensions and Sort Values
%Ordered values.
[m,n]=size(dem);
[D,I]=sort(dem(:));

clear dem;

%%Find max index of non-null values%%
D_nonull=D';
D_nonull= D_nonull(~isnan(D'));
max_elev = max(D_nonull);
min_elev = min(D_nonull);
max_idx=length(D_nonull);

clear D_nonull;

%% EQUATION %%
% Let v_diag = v_n - v_1 and p_i = < v_i, v_diag>/|v_diag|^2
% Data Driven Color Mapping by Martin Eisemann, Georgia Albuquerque, Marcus
% Magnor
% https://pdfs.semanticscholar.org/1145/6c8c21cd2b740561b6b82516fd61ce125051.pdf

% x is position in sorted array
% y is data value

% To change angle of diagonal:
% vdiag = vn - v1 = [xn,yn] - [x1,y1] = [xn-x1,yn-y1] = [dx,dy]
% dx = (yn-y1)/tan(a)

%% Data Projection
 
% Calculate v_diag & V:
%dx = max_idx-1;
dx = (D(max_idx)-D(1))/tand(a);
dy = (D(max_idx)-D(1));
v_diag=[dx, dy];

denominator = norm(v_diag)*norm(v_diag); %Square of magnitude of v_diag.

v_idx=(1:1:length(D))/length(D)*dx; %vector of positions (indices) (why??)
%v_idx=( 1:1:length(D) ) ; %vector of positions (indices)
v_i=[v_idx(:),D(:)]; % 2 by n matrix v_i

V_diag=repmat(v_diag,length(D),1);
P=dot(v_i,V_diag,2)/denominator;
%min(P)
%max(P)

clear v_i v_diag V_diag denominator v_min v_max v_idx;

%% Put Values back into Matrix
C(I)=P-min(P);
color_dem=reshape(C,[m,n]);
 
%% Colormaps 

% % Starry Night
% sn_map=[0,0,102;60,83,163;94,192,238;174,245,255;255,255,255;253,245,124;251,217,40;240,165,0;200,90,0;162,40,10;88,1,1];
% sn_norm_map=sn_map/255;
% 
% % Eyeball
% eb_map=[109,61,103;218,122,206;190,110,158;170,113,217;134,79,182;81,30,145;32,16,200;16,64,200;16,110,220;32,150,240;32,190,255;32,224,255;90,225,245;60,190,220;120,190,220;170,190,220;210,220,210;150,220,150;66,200,75;90,175,75;120,150,85;100,100,60;130,80,20;160,100,25;190,130,30;230,170,20;255,239,10;250,140,0;230,70,0;255,0,0;250,95,95;250,180,155;250,250,250];
% eb_norm_map=eb_map/255;

% ColorBrewer
% map=[30,83,151;69,117,180;97,144,195;125,172,209;154,197,223;185,217,234;215,238,245;233,247,233;255,255,191;255,244,174;255,233,157;254,215,138;253,183,116;253,150,95;243,116,76;229,82,57;215,48,39];
% norm_map=map/255;
% 
% %Normal
% n_map=[255,120,255;120,120,255;120,255,255;120,255,120;255,255,120;255,120,120;255,255,255];
% n_norm_map=n_map/255;

% Viridis (17)
% norm_map=[0.267004, 0.004874, 0.329415;0.282327, 0.094955, 0.417331;0.278826, 0.175490, 0.483397;0.258965, 0.251537, 0.524736;0.229739, 0.322361, 0.545706;0.199430, 0.387607, 0.554642;0.172719, 0.448791, 0.557885;0.149039, 0.508051, 0.557250;0.127568, 0.566949, 0.550556;0.120638, 0.625828, 0.533488;0.157851, 0.683765, 0.501686;0.246070, 0.738910, 0.452024;0.369214, 0.788888, 0.382914;0.515992, 0.831158, 0.294279;0.678489, 0.863742, 0.189503;0.845561, 0.887322, 0.099702;0.993248, 0.906157, 0.143936];
% map=norm_map*255;

% %Magma (18)0.001462, 0.000466, 0.013866;
% norm_map=[0.035520, 0.028397, 0.125209;0.102815, 0.063010, 0.257854;0.191460, 0.064818, 0.396152;0.291366, 0.064553, 0.475462;0.384299, 0.097855, 0.501002;0.475780, 0.134577, 0.507921;0.569172, 0.167454, 0.504105;0.664915, 0.198075, 0.488836;0.761077, 0.231214, 0.460162;0.852126, 0.276106, 0.418573;0.925937, 0.346844, 0.374959;0.969680, 0.446936, 0.360311;0.989363, 0.557873, 0.391671;0.996580, 0.668256, 0.456192;0.996727, 0.776795, 0.541039;0.992440, 0.884330, 0.640099;0.987053, 0.991438, 0.749504];
% map=norm_map*255;


%%  Display Images:

%TODO:
% 1) Move below final
% 2) Add display functionality for even bins
if sum(ismember(varargin,{'display'}))>0
    figure,image(color_dem,'CDataMapping','scaled')
    colormap(norm_map)
end

%% Initialize Elevation/Bin Related Variables:

elev_range = max_elev - min_elev;
step = (max(P)-min(P))/(size(map,1)); 
elevs = zeros(1,size(map,1)+1);

%% FIGURE OUT PRELIMINARY ELEVATION VALUES:
%add min & max elev and pcts
elevs(1)=min_elev;
%fill in the rest of the elevations

for l = 0:size(map,1) %17
    [~,closestIndex] = min(abs((P-min(P))-(step*l))); %Remember 0<=P<=1
    elevs(l+1)=D(closestIndex);
end
if sum(ismember(varargin,{'debug'}))>0
    fprintf('Original Elevations: \n');
    elevs/1000
end
%% Figure out Elevation Differences, and Min Bin Size:
elev_diffs=diff(elevs);
min_diff=min(elev_diffs); %TODO: Probably compare this to min bin size for final rounding.
%max_diff=max(elev_diffs); 
%avg_diff=avg(elev_diffs);
%ROUNDING OPTIONS:
round_opts=[1,5,10,20,25,50,100]; %potential ways to round. Can be edited.
%try and find out which minimum bin size is closest to the size of the
%smallest difference:
round_diffs=round_opts(:)-(min_diff);
round_diffs(round_diffs <= -1) = inf; 
[~,round_I]=min(abs(round_diffs)); %rounding option that's > or ~= to minimum

%make sure min_bin_size isn't larger than average and max bin size??
while elev_range/size(map,1) < round_opts(round_I) 
    round_I=round_I-1;
end

%make sure legend can handle minimum bin size:
min_percent=0.02*elev_range; %min percentage of range so legend is readable.
if min_percent > round_opts(round_I)
    min_percent_diffs=round_opts(:)-min_percent;
    min_percent_diffs(min_percent_diffs <= -1) = inf; 
     [~,round_I]=min(min_percent_diffs); %rounding option that's > or ~= to minimu
    min_value=min_percent;
else
    min_value=min_diff;
end

if round_opts(round_I) > elev_range/size(map,1)
    fprintf('ERROR: available minimum bin sizes are incompatible with minimum bin for legend, and uneven bins. Use a colormap with less colors.');
    return;
end

if round_I == 1
    %NEAREST 1:
    min_bin_size=round(min_value); %Do I need this??
    if min_bin_size == 0
        min_bin_size=1;
    end
elseif round_I == 2
    %NEAREST 5:
    min_bin_size=round(min_value/5)*5; %Do I need this??
    if min_bin_size == 0
        min_bin_size=5;
    end
elseif round_I == 3
    %NEAREST 10:
    min_bin_size=round(min_value,-1); %Do I need this??
    if min_bin_size == 0
        min_bin_size=10;
    end
elseif round_I == 4
    %NEAREST 20:
    min_bin_size=round(min_value/20,-1)*20; %Do I need this??
    if min_bin_size == 0
        min_bin_size=20;
    end
elseif round_I == 5
    %NEAREST 25:
    min_bin_size=round(min_value/25,-1)*25; %Do I need this??
    if min_bin_size == 0
        min_bin_size=25;
    end
elseif round_I == 6
    %NEAREST 50:
    min_bin_size=round(min_value/50,-1)*50; %Do I need this??
    if min_bin_size == 0
        min_bin_size=50;
    end
elseif round_I == 7
    %NEAREST 100:
    min_bin_size=round(min_value,-2);
    if min_bin_size == 0
        min_bin_size=100;
    end
else
    fprintf('Oh no! Something went wrong with determining how to round\n');
    return;
end

clear round_diffs round_opts;

%% Which Bins are Too Small - Initial Pass
S = Stack(size(elev_diffs,2)); %stack for center bin numbers
Dir = Stack(size(elev_diffs,2)); %stack for direction to shift bins

%TODO: Maybe too_small is not necessary?
too_small=elev_diffs<min_bin_size;%logical array of bins that are too small: 1=true.

k=1;
while k < size(elev_diffs,2)+1
    if too_small(k) == 1 %find an instance of the first in a streak of too small bins.
        j=k; %elev_diffs index
        b=0; %Number of too small bins.
        while j < size(elev_diffs,2)+1 %TODO make this & if statement more efficient
            if too_small(j) == 0 %find the last instance of the streak of too small bins
                %Create a stack (FILO) of the center bins and the
                %directions that they need to be adjusted:
                %   -1: Adjust the Bin Minimum
                %    0: Adjust both Bin Min and Max
                %    1: Adjust the Bin Maximum
                if mod(b,2) == 0 %even number of too small bins
                    center_bin_big=j-b + (b/2);
                    center_bin_sm=center_bin_big-1;
                    
                    S.Push(center_bin_big);
                    Dir.Push(1);
                    S.Push(center_bin_sm);
                    Dir.Push(-1);
                else %odd number of too small bins
                    center_bin_num=j-b + floor(b/2);
                    S.Push(center_bin_num);
                    Dir.Push(0);
                end
                %b=0;%reset b
                k=j;
                break;%to exit while loop.
            end
            b = b+1;
            j = j+1;
        end
        
    end
    k = k + 1;
end


%% Adjust Bin Sizes
elevs_adj=elevs; %Initialize array for adjusted elevs
%elevs_adj=round_elevs; %TRY THIS FIRST

while ~S.IsEmpty()
    bin_num=S.Pop();
    expandDir=Dir.Pop();
    
    if expandDir == 0
        %TODO - add check if max bin, so following if doesn't break??
        if (elevs_adj(bin_num+1)-elevs_adj(bin_num)) < min_bin_size
            if bin_num == 1 %Min Edge Case
                S.Push(bin_num);
                Dir.Push(1);
            elseif bin_num == size(elevs_adj,2)-1 %Max Edge Case
                S.Push(bin_num);
                Dir.Push(-1);
            else
                binCenter=((elevs_adj(bin_num+1)+elevs_adj(bin_num))/2);
                
                new_min=binCenter-(min_bin_size/2);
                if new_min <= min_elev
                    elevs_adj(bin_num)=min_elev;
                    S.Push(bin_num);
                    Dir.Push(1);
                else
                    elevs_adj(bin_num)=new_min;
                    S.Push(bin_num-1);
                    Dir.Push(-1);
                end
                
                
                new_max=binCenter+(min_bin_size/2);
                if new_max >= max_elev
                    elevs_adj(bin_num+1)=max_elev;
                    S.Push(bin_num);
                    Dir.Push(-1);
                else
                    elevs_adj(bin_num+1)=new_max;
                    S.Push(bin_num+1);
                    Dir.Push(1);
                end
                
                
            end
        end
    elseif expandDir == 1
        if (elevs_adj(bin_num+1)-elevs_adj(bin_num)) < min_bin_size
            if bin_num == size(elevs_adj,2)-1
                S.Push(bin_num);
                Dir.Push(-1);
            else
                new_max=elevs_adj(bin_num)+min_bin_size;
                if new_max >= max_elev
                    elevs_adj(bin_num+1)=max_elev;
                    S.Push(bin_num+1);
                    Dir.Push(-1);
                else
                    elevs_adj(bin_num+1)=new_max;
                    S.Push(bin_num+1);
                    Dir.Push(1);
                end
            end
        end
    elseif expandDir == -1
        if (elevs_adj(bin_num+1)-elevs_adj(bin_num)) < min_bin_size
            if bin_num == 1 %Min Edge Case
                S.Push(bin_num);
                Dir.Push(1);
            else
                new_min=elevs_adj(bin_num+1)-min_bin_size;
                if new_min <= min_elev
                    elevs_adj(bin_num)= min_elev;
                    S.Push(bin_num-1);
                    Dir.Push(1);
                else
                    elevs_adj(bin_num)= new_min;
                    S.Push(bin_num-1);
                    Dir.Push(-1);
                end
            end
        end
    else
        fprintf('Oh god. Something bad has happened. You might want to fix that.\n');
        return;
    end
end
if sum(ismember(varargin,{'debug'}))>0
    fprintf('Adjusted Elevations: \n')
    elevs_adj/1000
end

%% Round Results

final_elevs=elevs_adj; %initial array
max_round_bin=size(elevs_adj,2)-1;
if round_I == 1
    %NEAREST 1:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin),0);
elseif round_I == 2
    %NEAREST 5:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin)/5)*5;
elseif round_I == 3
    %NEAREST 10:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin),-1);
elseif round_I == 4
    %NEAREST 20:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin)/20)*20;
elseif round_I == 5
    %NEAREST 25:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin)/25)*25;
elseif round_I == 6
    %NEAREST 50:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin)/50)*50;
elseif round_I == 7
    %NEAREST 100:
    final_elevs(2:max_round_bin)=round(elevs_adj(2:max_round_bin),-2);
else
    fprintf('Oh no! Something went wrong with rouding\n');
    return;
end
if sum(ismember(varargin,{'debug'}))>0
    fprintf('Final Rounded Elevations: \n')
    final_elevs/1000
    fprintf('Differences: \n')
    diff(final_elevs)
end
%% Calculate Percentages & Make Output File:
pcts = (final_elevs-min_elev)/elev_range * 100;

[~,fbase,~]=fileparts(dem_name);
fname = strcat(fbase,'.lut');

fileID = fopen(fname,'w');
fprintf(fileID,'nv\t0\t0\t0 //noData to black\n'); %NO VALUE LINE
fprintf(fileID,'0%%\t%i\t%i\t%i\n',map(1,1),map(1,2),map(1,3)); %0 Percent Line
for e = 2:size(pcts,2)-1
    fprintf(fileID, '%.4f%%\t%i\t%i\t%i\n',pcts(e),map(e-1,1),map(e-1,2),map(e-1,3));
    fprintf(fileID, '%.4f0001%% %.0f %.0f %.0f \n', pcts(e), map(e,1),map(e,2),map(e,3));
end
fprintf(fileID,'100%%\t%i\t%i\t%i\n',map(size(map,1),1),map(size(map,1),2),map(size(map,1),3)); %100 Percent

fclose(fileID);
%% Print Min/Max Elevs to a Log File:
lname = strcat(fbase,'.log');
logID = fopen(lname,'w');
fprintf(logID,'Min:%.4f%\n',min_elev);
fprintf(logID,'Max:%.4f%\n',max_elev);
fclose(logID);
