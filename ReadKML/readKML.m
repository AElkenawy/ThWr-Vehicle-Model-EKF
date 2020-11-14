function [ cx,cy, s ] = readKML( fileName, idx )

if nargin==0, fileName='Route1.kml'; end
if nargin<2, idx=1; end

route = kml2struct(fileName);
[cx,cy] = wgs2utm(route(idx).coordMat(:,2),route(idx).coordMat(:,1));
cx=cx-cx(1);
cy=cy-cy(1);
plot(cx,cy,'.-');
grid on;
axis equal
s = [0;cumsum(sqrt(diff(cx).^2+diff(cy).^2))];

S = spline(s',[cx,cy]');
ds = 10;
sn = 1:ds:s(end);
cn = ppval(S,sn );
hold on
plot(cn(1,:),cn(2,:),'o-r');

% cx=cn(1,:);
% cy=cn(2,:);
% s = sn;
end

function kmlStruct = kml2struct(kmlFile)
% kmlStruct = kml2struct(kmlFile)
%
% Import a .kml file as a vector array of shapefile structs, with Geometry, Name,
% Description, Lon, Lat, and BoundaryBox fields.  Structs may contain a mix
% of points, lines, and polygons.
%
% .kml files with folder structure will not be presented as such, but will
% appear as a single vector array of structs.
%
% 

[FID msg] = fopen(kmlFile,'rt');

if FID<0
    error(msg)
end

txt = fread(FID,'uint8=>char')';
fclose(FID);

expr = '<Placemark.+?>.+?</Placemark>';

objectStrings = regexp(txt,expr,'match');

Nos = length(objectStrings);

for ii = 1:Nos
    % Find Object Name Field
    bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
    if isempty(bucket)
        name = 'undefined';
    else
        % Clip off flags
        name = regexprep(bucket{1},'<name.*?>\s*','');
        name = regexprep(name,'\s*</name>','');
    end
    
    % Find Object Description Field
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    if isempty(bucket)
        desc = '';
    else
        % Clip off flags
        desc = regexprep(bucket{1},'<description.*?>\s*','');
        desc = regexprep(desc,'\s*</description>','');
    end
    
    geom = 0;
    % Identify Object Type
    if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
        geom = 1;
    elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
        geom = 2;
    elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
        geom = 3;
    end
    
    switch geom
        case 1
            geometry = 'Point';
        case 2
            geometry = 'Line';
        case 3
            geometry = 'Polygon';
        otherwise
            geometry = '';
    end
    
    % Find Coordinate Field
    bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
    % Clip off flags
    coordStr = regexprep(bucket{1},'<coordinates.*?>(\s+)*','');
    coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
    % Split coordinate string by commas or white spaces, and convert string
    % to doubles
    coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
    % Rearrange coordinates to form an x-by-3 matrix
    [m,n] = size(coordMat);
    coordMat = reshape(coordMat,3,m*n/3)';
    [x,y]=ll2utm(coordMat(:,[1,2]));

    % Create structure
    kmlStruct(ii).Geometry = geometry;
    kmlStruct(ii).Name = name;
    kmlStruct(ii).Description = desc;
    kmlStruct(ii).coordMat=coordMat;
    kmlStruct(ii).c = [x,y]; 
end

end

function [x,y,f]=ll2utm(varargin)
%LL2UTM Lat/Lon to UTM coordinates precise conversion.
%	[X,Y]=LL2UTM2(LAT,LON) or LL2UTM([LAT,LON]) converts coordinates 
%	LAT,LON (in degrees) to UTM X and Y (in meters). Default datum is WGS84.
%
%	LAT and LON can be scalars, vectors or matrix. Outputs X and Y will
%	have the same size as inputs.
%
%	LL2UTM(...,DATUM) uses specific DATUM for conversion. DATUM can be one
%	of the following char strings:
%		'wgs84': World Geodetic System 1984 (default)
%		'nad27': North American Datum 1927
%		'clk66': Clarke 1866
%		'nad83': North American Datum 1983
%		'grs80': Geodetic Reference System 1980
%		'int24': International 1924 / Hayford 1909
%	or DATUM can be a 2-element vector [A,F] where A is semimajor axis (in
%	meters)	and F is flattening of the user-defined ellipsoid.
%
%	LL2UTM(...,ZONE) forces the UTM ZONE (scalar integer or same size as
%   LAT and LON) instead of automatic set.
%
%	[X,Y,ZONE]=LL2UTM(...) returns also the computed UTM ZONE (negative
%	value for southern hemisphere points).
%
%
%	XY=LL2UTM(...) or without any output argument returns a 2-column 
%	matrix [X,Y].
%
%	Note:
%		- LL2UTM does not perform cross-datum conversion.
%		- precision is near a millimeter.
%
%
%	Reference:
%		I.G.N., Projection cartographique Mercator Transverse: Algorithmes,
%		   Notes Techniques NT/G 76, janvier 1995.
%
%	Acknowledgments: Mathieu, Frederic Christen.
%
%
%	Author: Francois Beauducel, <beauducel@ipgp.fr>
%	Created: 2003-12-02
%	Updated: 2015-01-29


%	Copyright (c) 2001-2015, Franois Beauducel, covered by BSD License.
%	All rights reserved.
%
%	Redistribution and use in source and binary forms, with or without 
%	modification, are permitted provided that the following conditions are 
%	met:
%
%	   * Redistributions of source code must retain the above copyright 
%	     notice, this list of conditions and the following disclaimer.
%	   * Redistributions in binary form must reproduce the above copyright 
%	     notice, this list of conditions and the following disclaimer in 
%	     the documentation and/or other materials provided with the distribution
%	                           
%	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
%	AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
%	IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
%	ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
%	LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
%	CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
%	SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
%	INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
%	CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
%	ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
%	POSSIBILITY OF SUCH DAMAGE.

% Available datums
datums = [ ...
	{ 'wgs84', 6378137.0, 298.257223563 };
	{ 'nad83', 6378137.0, 298.257222101 };
	{ 'grs80', 6378137.0, 298.257222101 };
	{ 'nad27', 6378206.4, 294.978698214 };
	{ 'int24', 6378388.0, 297.000000000 };
	{ 'clk66', 6378206.4, 294.978698214 };
];

% constants
D0 = 180/pi;	% conversion rad to deg
K0 = 0.9996;	% UTM scale factor
X0 = 500000;	% UTM false East (m)

% defaults
datum = 'wgs84';
zone = [];

if nargin < 1
	error('Not enough input arguments.')
end

if nargin > 1 && isnumeric(varargin{1}) && isnumeric(varargin{2}) && all(size(varargin{1})==size(varargin{2}))
	lat = varargin{1};
	lon = varargin{2};
	v = 2;
elseif isnumeric(varargin{1}) && size(varargin{1},2) == 2
	lat = varargin{1}(:,1);
	lon = varargin{1}(:,2);
	v = 1;
else
	error('Single input argument must be a 2-column matrix [LAT,LON].')
end

if all([numel(lat),numel(lon)] > 1) && any(size(lat) ~= size(lon))
	error('LAT and LON must be the same size or scalars.')
end

for n = (v+1):nargin
	% LL2UTM(...,DATUM)
	if ischar(varargin{n}) || (isnumeric(varargin{n}) && numel(varargin{n})==2)
		datum = varargin{n};
	% LL2UTM(...,ZONE)
	elseif isnumeric(varargin{n}) && (isscalar(varargin{n}) || all(size(varargin{n})==size(lat)))
			zone = round(varargin{n});
	else
		error('Unknown argument #%d. See documentation.',n)
	end
end

if ischar(datum)
	% LL2UTM(...,DATUM) with DATUM as char
	if ~any(strcmpi(datum,datums(:,1)))
		error('Unkown DATUM name "%s"',datum);
	end
	k = find(strcmpi(datum,datums(:,1)));
	A1 = datums{k,2};
	F1 = datums{k,3};	
else
	% LL2UTM(...,DATUM) with DATUM as [A,F] user-defined
	A1 = datum(1);
	F1 = datum(2);
end

p1 = lat/D0;			% Phi = Latitude (rad)
l1 = lon/D0;			% Lambda = Longitude (rad)

% UTM zone automatic setting
if isempty(zone)
	F0 = round((l1*D0 + 183)/6);
else
	F0 = zone;
end

B1 = A1*(1 - 1/F1);
E1 = sqrt((A1*A1 - B1*B1)/(A1*A1));
P0 = 0/D0;
L0 = (6*F0 - 183)/D0;	% UTM origin longitude (rad)
Y0 = 1e7*(p1 < 0);		% UTM false northern (m)
N = K0*A1;

C = coef(E1,0);
B = C(1)*P0 + C(2)*sin(2*P0) + C(3)*sin(4*P0) + C(4)*sin(6*P0) + C(5)*sin(8*P0);
YS = Y0 - N*B;

C = coef(E1,2);
L = log(tan(pi/4 + p1/2).*(((1 - E1*sin(p1))./(1 + E1*sin(p1))).^(E1/2)));
z = complex(atan(sinh(L)./cos(l1 - L0)),log(tan(pi/4 + asin(sin(l1 - L0)./cosh(L))/2)));
Z = N.*C(1).*z + N.*(C(2)*sin(2*z) + C(3)*sin(4*z) + C(4)*sin(6*z) + C(5)*sin(8*z));
xs = imag(Z) + X0;
ys = real(Z) + YS;

% outputs zone if needed: scalar value if unique, or vector/matrix of the
% same size as x/y in case of crossed zones
if nargout > 2
   	f = F0.*sign(lat);
	fu = unique(f);
	if isscalar(fu)
		f = fu;
	end
end

if nargout < 2
	x = [xs(:),ys(:)];
else
	x = xs;
	y = ys;
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function c = coef(e,m)
%COEF Projection coefficients
%	COEF(E,M) returns a vector of 5 coefficients from:
%		E = first ellipsoid excentricity
%		M = 0 for transverse mercator
%		M = 1 for transverse mercator reverse coefficients
%		M = 2 for merdian arc


if nargin < 2
	m = 0;
end

switch m
	case 0
	c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
           -105/4096, 0, -45/1024, 0,  -3/32, 0, -3/8, 0, 0;
           525/16384, 0,  45/1024, 0, 15/256, 0,    0, 0, 0;
          -175/12288, 0, -35/3072, 0,      0, 0,    0, 0, 0;
          315/131072, 0,        0, 0,      0, 0,    0, 0, 0];
	  
	case 1
	c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
             1/61440, 0,   7/2048, 0,   1/48, 0,  1/8, 0, 0;
          559/368640, 0,   3/1280, 0,  1/768, 0,    0, 0, 0;
          283/430080, 0, 17/30720, 0,      0, 0,    0, 0, 0;
       4397/41287680, 0,        0, 0,      0, 0,    0, 0, 0];

	case 2
	c0 = [-175/16384, 0,   -5/256, 0,  -3/64, 0, -1/4, 0, 1;
         -901/184320, 0,  -9/1024, 0,  -1/96, 0,  1/8, 0, 0;
         -311/737280, 0,  17/5120, 0, 13/768, 0,    0, 0, 0;
          899/430080, 0, 61/15360, 0,      0, 0,    0, 0, 0;
      49561/41287680, 0,        0, 0,      0, 0,    0, 0, 0];
   
end
c = zeros(size(c0,1),1);

for i = 1:size(c0,1)
    c(i) = polyval(c0(i,:),e);
end

end


function  [x,y,utmzone,utmhemi] = wgs2utm(Lat,Lon,utmzone,utmhemi)
% -------------------------------------------------------------------------
% [x,y,utmzone] = wgs2utm(Lat,Lon,Zone)
%
% Description:
%    Convert WGS84 coordinates (Latitude, Longitude) into UTM coordinates
%    (northing, easting) according to (optional) input UTM zone and
%    hemisphere.
%
% Input:
%    Lat: WGS84 Latitude scalar, vector or array in decimal degrees.
%    Lon: WGS84 Longitude scalar, vector or array in decimal degrees.
%    utmzone (optional): UTM longitudinal zone. Scalar or same size as Lat
%       and Lon.
%    utmhemi (optional): UTM hemisphere as a single character, 'N' or 'S',
%       or array of 'N' or 'S' characters of same size as Lat and Lon.
%
% Output:
%    x: UTM easting in meters.
%    y: UTM northing in meters.
%    utmzone: UTM longitudinal zone.
%    utmhemi: UTM hemisphere as array of 'N' or 'S' characters.
%
% Author notes:
%    I downloaded and tried deg2utm.m from Rafael Palacios but found
%    differences of up to 1m with my reference converters in southern
%    hemisphere so I wrote my own code based on "Map Projections - A
%    Working Manual" by J.P. Snyder (1987). Quick quality control performed
%    only by comparing with LINZ converter
%    (www.linz.govt.nz/apps/coordinateconversions/) and Chuck Taylor's
%    (http://home.hiwaay.net/~taylorc/toolbox/geography/geoutm.html) on a
%    few test points, so use results with caution. Equations not suitable
%    for a latitude of +/- 90deg.
%
%    UPDATE: Following requests, this new version allows forcing UTM zone
%    in input.
%
% Examples:
%
%    % set random latitude and longitude arrays
%    Lat= 90.*(2.*rand(3)-1)
%    Lon= 180.*(2.*rand(3)-1)
%
%    % let the function find appropriate UTM zone and hemisphere from data
%    [x1,y1,utmzone1,utmhemi1] = wgs2utm(Lat,Lon)
%
%    % forcing unique UTM zone and hemisphere for all data entries
%    % note: resulting easting and northing are way off the usual values
%    [x2,y2,utmzone2,utmhemi2] = wgs2utm(Lat,Lon,60,'S')
%
%    % forcing different UTM zone and hemisphere for each data entry
%    % note: resulting easting and northing are way off the usual values
%    utmzone = floor(59.*rand(3))+1
%    utmhemi = char(78 + 5.*round(rand(3)))
%    [x3,y3,utmzone3,utmhemi3] = wgs2utm(Lat,Lon,utmzone,utmhemi)
%
% Author:
%   Alexandre Schimel
%   MetOcean Solutions Ltd
%   New Plymouth, New Zealand
%
% Version 2:
%   February 2011
%-------------------------------------------------------------------------

%% Argument checking
if ~sum(double(nargin==[2,4]))
    error('Wrong number of input arguments');return
end
n1=size(Lat);
n2=size(Lon);
if (n1~=n2)
    error('Lat and Lon should have same size');return
end
if exist('utmzone','var') && exist('utmhemi','var')
    n3=size(utmzone);
    n4=size(utmhemi);
    if (sort(n3)~=sort(n4))
        error('utmzone and utmhemi should have same size');return
    end
    if max(n3)~=1 && max(n3)~=max(n1)
        error('utmzone should have either same size as Lat and Long, or size=1');return
    end
end

% expand utmzone and utmhemi if needed
if exist('utmzone','var') && exist('utmhemi','var')
    n3=size(utmzone);
    n4=size(utmhemi);
    if n3==[1 1]
        utmzone = utmzone.*ones(size(Lat));
        utmhemi = char(utmhemi.*ones(size(Lat)));
    end
end

%% coordinates in radians
lat = Lat.*pi./180;
lon = Lon.*pi./180;

%% WGS84 parameters
a = 6378137;           %semi-major axis
b = 6356752.314245;    %semi-minor axis
% b = 6356752.314140;  %GRS80 value, originally used for WGS84 before refinements
e = sqrt(1-(b./a).^2); % eccentricity

%% UTM parameters
% lat0 = 0;                % reference latitude, not used here
if exist('utmzone','var')
    Lon0 = 6.*utmzone-183; % reference longitude in degrees
else
    Lon0 = floor(Lon./6).*6+3; % reference longitude in degrees
end
lon0 = Lon0.*pi./180;      % in radians
k0 = 0.9996;               % scale on central meridian

FE = 500000;              % false easting
if exist('utmhemi','var')
    FN = double(utmhemi=='S').*10000000;
else
    FN = (Lat < 0).*10000000; % false northing
end

%% Equations parameters
eps = e.^2./(1-e.^2);  % e prime square
% N: radius of curvature of the earth perpendicular to meridian plane
% Also, distance from point to polar axis
N = a./sqrt(1-e.^2.*sin(lat).^2);
T = tan(lat).^2;
C = ((e.^2)./(1-e.^2)).*(cos(lat)).^2;
A = (lon-lon0).*cos(lat);
% M: true distance along the central meridian from the equator to lat
M = a.*(  ( 1 - e.^2./4 - 3.*e.^4./64 - 5.*e.^6./256 )  .* lat         ...
    -( 3.*e.^2./8 + 3.*e.^4./32 + 45.*e.^6./1024 ) .* sin(2.*lat) ...
    +( 15.*e.^4./256 + 45.*e.^6./1024 )            .* sin(4.*lat) ...
    -(35.*e.^6./3072 )                             .* sin(6.*lat) );

%% easting
x = FE + k0.*N.*(                                  A       ...
    + (1-T+C)                      .* A.^3./6 ...
    + (5-18.*T+T.^2+72.*C-58.*eps) .* A.^5./120 );

%% northing
% M(lat0) = 0 so not used in following formula
y = FN + k0.*M + k0.*N.*tan(lat).*(                                     A.^2./2  ...
    + (5-T+9.*C+4.*C.^2)              .* A.^4./24 ...
    + (61-58.*T+T.^2+600.*C-330.*eps) .* A.^6./720 );

%% UTM zone
if exist('utmzone','var') && exist('utmhemi','var')
    utmzone = utmzone;
    utmhemi = utmhemi;
else
    utmzone = floor(Lon0./6)+31;
    utmhemi = char( 83.* (Lat < 0) + 78.* (Lat >= 0) );
end

end