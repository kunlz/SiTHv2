function [] = write2nc(filename, varname, longname, VAL, days, uni)
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

nccreate(filename, varname, ...
    'Dimensions', {'lat', dem_lat, 'lon', dem_lon, 'time', days}, ...
    'Datatype', datatype, ...
    'DeflateLevel', 9);

ncwriteatt(filename, varname, 'long_name', longname);
ncwriteatt(filename, varname, 'scale factor', '0.01'); 
ncwriteatt(filename, varname, 'unit', uni); 
ncwriteatt(filename, varname, 'resolution', '0.1 degree');

ncwriteatt(filename, '/', 'SiTH model', 'Simple Terrestrial Hydrosphere model');
ncwriteatt(filename, '/', 'Fill Value', '-');

% write to nc file
ncwrite(filename, varname, VAL);
ncwrite(filename, 'time', time);
ncwrite(filename, 'lat', lat1);
ncwrite(filename, 'lon', lon1);

end
