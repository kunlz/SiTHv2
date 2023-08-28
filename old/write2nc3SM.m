function [] = write2nc3SM(filename, varname1, longname1,...
    varname2, longname2, varname3, longname3, VAL1, VAL2, VAL3, days, uni)
% ---------
% filename : e.g., 'xxx2012.nc'
% varname  : e.g., 'Et'
% longname : e.g., 'Actual Evapotranspiration'
% VAL      : e.g., real variable name -- X_ET
% days     : e.g., length of time series
% uni      : e.g., 'mm'
% ----------
datatype = 'int16';
spa_res  = 0.1;

time = (1:1:days)';
lat2 = -90+spa_res/2:spa_res:90-spa_res/2;
lat1 = fliplr(lat2);
lon1 = -180+spa_res/2:spa_res:180-spa_res/2;

dem_lat = 180/spa_res;
dem_lon = 360/spa_res;

% create nc file
nccreate(filename, 'time', ...
    'Dimensions', {'time', days}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);
ncwriteatt(filename, 'time', 'long_name', 'Time');
ncwriteatt(filename, 'time', 'units', 'days');

nccreate(filename, 'lat', ...
    'Dimensions', {'lat', dem_lat}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);
ncwriteatt(filename, 'lat', 'long_name', 'Latitude');
ncwriteatt(filename, 'lat', 'units', 'degrees_north');

nccreate(filename, 'lon', ...
    'Dimensions', {'lon', dem_lon}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);
ncwriteatt(filename, 'lon', 'long_name', 'Longitude');
ncwriteatt(filename, 'lon', 'units', 'degrees_east');

% SM1
nccreate(filename, varname1, ...
    'Dimensions', {'lat', dem_lat, 'lon', dem_lon, 'time', days}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);

ncwriteatt(filename, varname1, 'long_name', longname1);
ncwriteatt(filename, varname1, 'scale factor', '0.01'); 
ncwriteatt(filename, varname1, 'unit', uni); 
ncwriteatt(filename, varname1, 'resolution', '0.1 degree');

% SM2
nccreate(filename, varname2, ...
    'Dimensions', {'lat', dem_lat, 'lon', dem_lon, 'time', days}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);

ncwriteatt(filename, varname2, 'long_name', longname2);
ncwriteatt(filename, varname2, 'scale factor', '0.01'); 
ncwriteatt(filename, varname2, 'unit', uni); 
ncwriteatt(filename, varname2, 'resolution', '0.1 degree');

% SM3
nccreate(filename, varname3, ...
    'Dimensions', {'lat', dem_lat, 'lon', dem_lon, 'time', days}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);

ncwriteatt(filename, varname3, 'long_name', longname3);
ncwriteatt(filename, varname3, 'scale factor', '0.01'); 
ncwriteatt(filename, varname3, 'unit', uni); 
ncwriteatt(filename, varname3, 'resolution', '0.1 degree');

% global attribution
ncwriteatt(filename, '/', 'SiTH model', 'Simple Terrestrial Hydrosphere model');
ncwriteatt(filename, '/', 'Fill Value', '-');

% write to nc file
ncwrite(filename, varname1, VAL1);
ncwrite(filename, varname2, VAL2);
ncwrite(filename, varname3, VAL3);
ncwrite(filename, 'time', time);
ncwrite(filename, 'lat', lat1);
ncwrite(filename, 'lon', lon1);

end
