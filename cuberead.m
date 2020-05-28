function [ image_data  ] = cuberead(image_name,varargin)
%% Reads in an ISIS cube file, tiled or band sequential.
% [ image_data  ] = cuberead(image_name);
% or 
% [ cropped_image_data  ] = cuberead(image_name, beginsample, endsample, beginline, endline);


switch nargin
    case 1;
        % input the image name only
        beginline=  -inf;
        beginsample=-inf;
        endline=     inf;
        endsample=   inf;
        beginband=  -inf;
        endband=     inf;
    case 7;
        beginsample =varargin{1};
        endsample   =varargin{2};
        beginline   =varargin{3};
        endline     =varargin{4};
        beginband   =varargin{5};
        endband     =varargin{6};
    otherwise
        ERR=MException( 'UserInput:IncorrectArgumentNumber','Please input the correct number of arguments\n1 or 5');
        throw(ERR)
end
%% gather cube information
try
    A=fopen(image_name);
catch ERR
    if ~exist(image_name, 'file')
    ERR=MException( 'IOError:FileNotFound',[ 'file: ' image_name ' not found' ]);
    throw(ERR)
    else
        rethrow(ERR)
    end
end
        

[ exitcode , format ] = system( ['grep -aim1 '' format '' ''' image_name  ''' | sed ''s#^.*=[[:blank:]]*\([[:alpha:]]*\).*#\1#''' ] );
if exitcode~=0
    ERR=MException( 'SystemCall:ExitNotZero','could not read file format');
    throw(ERR)
end
format=format(1:(end-1));

[ exitcode , samples ] = system( ['grep -aim1 '' samples '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
if exitcode~=0
    ERR = MException( 'SystemCall:ExitNotZero','could not read samples');
    throw(ERR)
end
    samples=str2double(samples);

[ exitcode , lines ] = system( ['grep -aim1 '' lines '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
if exitcode~=0
    ERR = MException( 'SystemCall:ExitNotZero','could not read Lines');
    throw(ERR)
end
    lines=str2double(lines);

[ exitcode , bands ] = system( ['grep -aim1 '' bands '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
if exitcode~=0
   ERR = MException( 'SystemCall:ExitNotZero', 'could not read Bands');
    throw(ERR)
end
    bands=str2double(bands);
    
[ exitcode , startbyte ] = system( ['grep -aim1 '' StartByte '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
if exitcode~=0
    ERR = MException( 'SystemCall:ExitNotZero','could not read header length');
    throw(ERR)
end
    startbyte=str2double(startbyte) - 1 ;
    

if strcmpi(format, 'Tile' )
    [ exitcode , tilesamples ] = system( ['grep -aim1 '' tilesamples '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
    if exitcode~=0
        ERR = MException( 'SystemCall:ExitNotZero','could not read TileSamples');
        throw(ERR)
    end
    tilesamples=str2double(tilesamples);
    
    [ exitcode , tilelines ] = system( ['grep -aim1 '' tilelines '' ''' image_name ''' | grep -Eo "[0-9]+" ' ] );
    if exitcode~=0
        ERR = MException( 'SystemCall:ExitNotZero','could not read TileLines');
        throw(ERR)
    end
    tilelines=str2double(tilelines);
    
elseif strcmpi(format, 'BandSequential' )
    tilesamples=1;
    tilelines=1;
else
    ERR = MException( 'SystemCall:ExitNotZero','could not read file format');
    throw(ERR)
end

[ exitcode , bittype ] = system( ['grep -aim1 '' Type '' ''' image_name ''' | sed ''s#^.*=[[:blank:]]*\([[:alpha:]]*\).*#\1#''' ] );
    if exitcode~=0
        ERR = MException( 'SystemCall:ExitNotZero','could not read TileLines');
        throw(ERR)
    end
    bittype=bittype(1:(end-1));
    
switch bittype
    case 'UnsignedByte'
        bytesperpix=1;
        bittype='uint8';
    case 'SignedWord'
        bytesperpix=2;
        bittype='int16';
    case 'Real'
        bytesperpix=4;
        bittype='single';
end

[ exitcode , base ] = system( ['grep -aim1 '' Base '' ''' image_name ''' | grep -Eo "[0-9.]+" ' ] );
if exitcode~=0
   ERR = MException( 'SystemCall:ExitNotZero', 'could not read Base');
    throw(ERR)
end
base=str2double(base);

[ exitcode , multiplier ] = system( ['grep -aim1 '' Multiplier '' ''' image_name ''' | grep -Eo "[0-9.]+" ' ] );
if exitcode~=0
   ERR = MException( 'SystemCall:ExitNotZero', 'could not read Multiplier');
    throw(ERR)
end
multiplier=str2double(multiplier);

%% test user inputs
% to be continued


%% read image data
imagelines=ceil(lines/tilelines)*tilelines;
imagesamples=ceil(samples/tilesamples)*tilesamples;
imagetiles=[ ceil(lines/tilelines) ceil(samples/tilesamples) ];

beginband=max([1 beginband]);
endband  =min([bands endband]);

firsttile=[floor(max([1 beginline])/ tilelines)+1 floor(max([1 beginsample])/ tilesamples)+1 ];
lasttile=[ceil(min([lines endline])/tilelines) ceil(min([samples endsample])/tilesamples) ];

fseekbob= (firsttile(1)-1) * imagesamples * tilelines;
fseekeob= ( imagetiles(1) - lasttile(1)) * imagesamples * tilelines;
fseekbol= (firsttile(2)-1) * tilesamples * tilelines;
fseekeol= (imagetiles(2) - lasttile(2)) * tilesamples * tilelines;

fseek(A,startbyte,-1);

V=nan( (lasttile(1) - firsttile(1)+1)*tilelines , (lasttile(2) - firsttile(2) + 1 ) * tilesamples,  endband-beginband+ 1 );
fseek(A, imagesamples*imagelines * (1-beginband) * bytesperpix ,0);
for i=1:(endband-beginband+ 1)
    fseek(A, fseekbob * bytesperpix ,0); % seek to current line
%     if strcmpi(format, 'Tile' )
        for j=1:(lasttile(1) - firsttile(1)+1)
            fseek(A,fseekbol * bytesperpix  ,0); %seek to begin of samples
            for k=1:(lasttile(2)-firsttile(2)+1)
                V( ((j-1)*tilelines + 1) : j*tilelines , ((k-1)*tilesamples + 1) : k*tilesamples , i)=fread(A,[ tilesamples tilelines ], bittype)';
            end
            fseek(A, fseekeol * bytesperpix ,0); %seek to end of the samples
        end
%     else
%         fseek(A, fseekeob * bytesperpix ,0); %seek to end of band
%         V(:,:,i)=fread(A,[ ceil(samples/tilesamples)*tilesamples ceil(lines/tilelines)*tilelines ], bittype)';
%         
%     end
    fseek(A, fseekeob * bytesperpix ,0);
%     [spec2nany spec2nanx]=find(V < -2000000000000);
%     for j=1:numel(spec2nany)
%         V(spec2nany(j),spec2nanx(j))=NaN;
%     end
end

V(V <= -32768)=NaN;

V=V( (max( [1 beginline]) - (firsttile(1)-1)*tilelines+1):( min([lines endline]) -(firsttile(1)-1)*tilelines ) , (max([1 beginsample]) - (firsttile(2)-1)*tilesamples + 1):( min([samples endsample]) -(firsttile(2)-1)*tilesamples ) , : ) * multiplier + base;

image_data=V;



% 
% moonref=referenceSphere('Moon');
% moonref.Radius=1737400;
% mstruct = defaultm('pcarree');
% mstruct.geoid = moonref;

