import math
# DH*
import numpy
# *DH

try:
    import netCDF4
    NETCDF4_PRESENT = True
except ImportError:
    NETCDF4_PRESENT = False

from mpl_toolkits.basemap import Basemap

__all__ = ['WrfBaseMap', 'WrfBaseMapFactory', 'GeogridBaseMap',
           'get_mass_basemap']


class UnsupportedProjectionError(Exception):
    pass


class UnsupportedEarthRadiusError(Exception):
    pass


def get_mass_basemap(ncfile, coastline_resolution=None):
    '''
    A helper function to create a basemap for the mass grid from either a
    geogrid file or a wrf output file.
    '''
    if 'XLONG_M' in ncfile.variables:
        return GeogridBaseMap(ncfile, coastline_resolution)
    elif 'XLONG' in ncfile.variables:
        return WrfBaseMapFactory(ncfile,
                                 coastline_resolution).get_mass_basemap()
    else:
        raise ValueError(
            'Unable to find XLONG_M or XLONG variables in ' +
            'NetCDF file.  Suspect file is not a WRF output file or' +
            'Geogrid file.'
        )


class LambertConformalBasemapParameters(object):
    # The WRF definition of earth radius
    EARTH_RADIUS_M = 6370000
    PROJECTION_NAME = 'lcc'

    def __init__(self, lat_0=None, lon_0=None, lat_1=None, lat_2=None,
                 llcrnrlon=None, llcrnrlat=None, urcrnrlon=None,
                 urcrnrlat=None, resolution=None):
        self.lat_0 = lat_0
        self.lon_0 = lon_0
        self.lat_1 = lat_1
        self.lat_2 = lat_2
        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        self.urcrnrlon = urcrnrlon
        self.urcrnrlat = urcrnrlat
        self.resolution = resolution

    def from_basemap(self, basemap):
        """Retrieve required projection parameters from a premade Basemap
        instance."""
        if not basemap.projection == self.PROJECTION_NAME:
            raise UnsupportedProjectionError(basemap.projection)
        if not basemap.projparams['R'] == self.EARTH_RADIUS_M:
            raise UnsupportedEarthRadiusError(
                'We are expecting te same earth radius as used by WRF ' +
                '(%d).  This basemap has an earth radius of %d.' % (
                    self.EARTH_RADIUS_M, basemap.projparams['R'])
            )
        self.lat_0 = basemap.projparams['lat_0']
        self.lon_0 = basemap.projparams['lon_0']
        self.lat_1 = basemap.projparams['lat_1']
        self.lat_2 = basemap.projparams['lat_2']
        self.llcrnrlon = basemap.llcrnrlon
        self.llcrnrlat = basemap.llcrnrlat
        self.urcrnrlon = basemap.urcrnrlon
        self.urcrnrlat = basemap.urcrnrlat
        self.resolution = basemap.resolution
        return self

    def from_geogrid_file(self, ncfile):
        """gets the projection parameters from an open geogrid netcdf file"""
        self.lon_0 = ncfile.STAND_LON
        self.lat_1 = ncfile.TRUELAT1
        self.lat_2 = ncfile.TRUELAT2
        self.llcrnrlon = ncfile.corner_lons[0]
        self.llcrnrlat = ncfile.corner_lats[0]
        self.urcrnrlon = ncfile.corner_lons[2]
        self.urcrnrlat = ncfile.corner_lats[2]
        return self

    def from_wrf_file(self, ncfile, longitude_variable_name,
                      latitude_variable_name):
        """gets the projection parameters from an open WRF netcdf file"""
        for attribute_name in ['MOAD_CEN_LAT', 'STAND_LON', 'TRUELAT1',
                               'TRUELAT2', 'MAP_PROJ']:
            if not hasattr(ncfile, attribute_name):
                raise ValueError(
                    "Input file does not have global attribute %s" % (
                        attribute_name)
                )
        if ncfile.MAP_PROJ != 1:
            raise ValueError(
                "MAP_PROJ of type %d is not supported" % (ncfile.MAP_PROJ)
            )
        for variable_name in [longitude_variable_name, latitude_variable_name]:
            if not variable_name in ncfile.variables:
                raise ValueError(
                    "Input file does not have the variable %s" % (
                        variable_name)
                )
        self.lat_0 = ncfile.MOAD_CEN_LAT
        self.lon_0 = ncfile.STAND_LON
        self.lat_1 = ncfile.TRUELAT1
        self.lat_2 = ncfile.TRUELAT2
        self.llcrnrlon = ncfile.variables[longitude_variable_name][0, 0, 0]
        self.llcrnrlat = ncfile.variables[latitude_variable_name][0, 0, 0]
        self.urcrnrlon = ncfile.variables[longitude_variable_name][0, -1, -1]
        self.urcrnrlat = ncfile.variables[latitude_variable_name][0, -1, -1]
        return self

# DH*
    def from_mm5_file(self, ncfile, longitude_variable_name,
                      latitude_variable_name):
        """gets the projection parameters from an open MM5 netcdf file"""

        def xtodot(slab_crs):
            """Adaptation of the Fortran routine to put cross grid points
            on dot grid points. Fairly Fortran-y syntax here."""
            maxiy = slab_crs.shape[0] + 1
            maxjx = slab_crs.shape[1] + 1
            bot = numpy.zeros([max(maxiy, maxjx)])
            rleft = numpy.zeros([max(maxiy, maxjx)])
            slab_dot = numpy.zeros([maxiy, maxjx])
            # Extrapolate out to top and bottom edges
            j = 1
            while j < maxjx - 1:
                bot[j] = (3. * (slab_crs[0, j - 1] + slab_crs[0, j]) -
                          (slab_crs[1, j - 1] + slab_crs[1, j])) / 4.
                slab_dot[maxiy - 1, j] = \
                    (3. * (slab_crs[maxiy - 2, j - 1] + slab_crs[maxiy - 2, j]) -
                     (slab_crs[maxiy - 3, j - 1] + slab_crs[maxiy - 3, j])) / 4.
                j += 1
            # Extrapolate out to left and right edges
            i = 1
            while i < maxiy - 1:
                rleft[i] = (3. * (slab_crs[i - 1, 0] + slab_crs[i, 0]) -
                            (slab_crs[i - 1, 1] + slab_crs[i, 1])) / 4.
                slab_dot[i, maxjx - 1] = \
                    (3. * (slab_crs[i - 1, maxjx - 2] + slab_crs[i, maxjx - 2]) -
                     (slab_crs[i - 1, maxjx - 3] + slab_crs[i, maxjx - 3])) / 4.
                i += 1
            # Extrapolate out to corners
            rleft[0] = (3. * slab_crs[0, 0] - slab_crs[1, 1]) / 2.
            rleft[
                maxiy - 1] = (3. * slab_crs[maxiy - 2, 0] - slab_crs[maxiy - 3, 1]) / 2.
            bot[maxjx - 1] = (3. * slab_crs[0, maxjx - 2] -
                              slab_crs[1, maxjx - 3]) / 2.
            slab_dot[maxiy - 1, maxjx - 1] = (3. * slab_crs[maxiy - 2, maxjx - 2] -
                                              slab_crs[maxiy - 3, maxjx - 3]) / 2.
            # Interpolate in the interior
            j = maxjx - 2
            while j > 0:
                i = maxiy - 2
                while i > 0:
                    slab_dot[i, j] = .25 * (slab_crs[i - 1, j - 1] + slab_crs[i, j - 1] +
                                            slab_crs[i - 1, j] + slab_crs[i, j])
                    i -= 1
                j -= 1
            # Put "bot" and "rleft" values into slab
            j = 0
            while j <= maxjx - 1:
                slab_dot[0, j] = bot[j]
                j += 1
            i = 0
            while i <= maxiy - 1:
                slab_dot[i, 0] = rleft[i]
                i += 1
            return slab_dot

        for variable_name in ['coarse_cenlat', 'stdlon', 'stdlat_1',
                              'stdlat_2', 'map_proj_code']:
            if not variable_name in ncfile.variables:
                raise ValueError(
                    "Input file does not have the variable %s" % (
                        variable_name)
                )

        if ncfile.variables['map_proj_code'].get_value() != 1:
            raise ValueError(
                "map_proj_code of type %d is not supported" %
                (ncfile['map_proj_code'].get_value())
            )
        for variable_name in [longitude_variable_name, latitude_variable_name]:
            if not variable_name in ncfile.variables:
                raise ValueError(
                    "Input file does not have the variable %s" % (
                        variable_name)
                )
        self.lat_0 = ncfile.variables['coarse_cenlat'].get_value()
        self.lon_0 = ncfile.variables['stdlon'].get_value()
        self.lat_1 = ncfile.variables['stdlat_1'].get_value()
        self.lat_2 = ncfile.variables['stdlat_2'].get_value()

        if longitude_variable_name == "longidot":
            longicrs = ncfile.variables["longicrs"]
            longidot = xtodot(longicrs)
            latitcrs = ncfile.variables["latitcrs"]
            latitdot = xtodot(latitcrs)
            self.llcrnrlon = longidot[0, 0]
            self.llcrnrlat = latitdot[0, 0]
            self.urcrnrlon = longidot[-1, -1]
            self.urcrnrlat = latitdot[-1, -1]
        else:
            self.llcrnrlon = ncfile.variables[longitude_variable_name][0, 0]
            self.llcrnrlat = ncfile.variables[latitude_variable_name][0, 0]
            self.urcrnrlon = ncfile.variables[longitude_variable_name][-1, -1]
            self.urcrnrlat = ncfile.variables[latitude_variable_name][-1, -1]

        return self
# *DH

    def to_basemap(self):
        """returns a matplotlib Basemap object"""
        return Basemap(projection=self.PROJECTION_NAME,
                       rsphere=self.EARTH_RADIUS_M,
                       resolution=self.resolution,
                       lat_0=self.lat_0, lon_0=self.lon_0,
                       lat_1=self.lat_1, lat_2=self.lat_2,
                       llcrnrlon=self.llcrnrlon, llcrnrlat=self.llcrnrlat,
                       urcrnrlon=self.urcrnrlon, urcrnrlat=self.urcrnrlat)

    def get_geogrid_params(self):
        """returns a dictionary of geogrid namelist settings corresponding to
        these basemnap parameters"""
        return {'map_proj': 'lambert', 'ref_lat': self.lat_0,
                'ref_lon': self.lon_0, 'truelat1': self.lat_1, 'truelat2': self.lat_2,
                'stand_lon': self.lon_0}

    def output(self):
        """Return 'nice' format output suitable for printing"""
        output = 'lat_0: {}\n'.format(self.lat_0)
        output += 'lat_1: {}\n'.format(self.lat_1)
        output += 'lat_2: {}\n'.format(self.lat_2)
        output += 'llcrnrlat: {}\n'.format(self.llcrnrlat)
        output += 'llcrnrlon: {}\n'.format(self.llcrnrlon)
        output += 'lon_0: {}\n'.format(self.lon_0)
        output += 'urcrnrlat: {}\n'.format(self.urcrnrlat)
        output += 'urcrnrlon: {}\n'.format(self.urcrnrlon)
        return output

class MercatorBasemapParameters(object):
    # The WRF definition of earth radius
    EARTH_RADIUS_M = 6370000
    PROJECTION_NAME = 'merc'

    def __init__(self, lat_0, lon_0, lat_ts, llcrnrlon, llcrnrlat, urcrnrlon,
                 urcrnrlat, resolution=None):
        self.lat_0 = lat_0
        self.lon_0 = lon_0
        self.lat_ts = lat_ts
        self.llcrnrlon = llcrnrlon
        self.llcrnrlat = llcrnrlat
        self.urcrnrlon = urcrnrlon
        self.urcrnrlat = urcrnrlat
        self.resolution = resolution

    def to_basemap(self):
        """returns a matplotlib Basemap object"""
        return Basemap(
            projection=self.PROJECTION_NAME, rsphere=self.EARTH_RADIUS_M,
            resolution=self.resolution, lat_0=self.lat_0, lon_0=self.lon_0,
            lat_ts=self.lat_ts, llcrnrlon=self.llcrnrlon, llcrnrlat=self.llcrnrlat,
            urcrnrlon=self.urcrnrlon, urcrnrlat=self.urcrnrlat)

    def get_geogrid_params(self):
        """returns a dictionary of geogrid namelist settings corresponding to
        these basemnap parameters"""
        return {'map_proj': 'mercator', 'ref_lat': self.lat_0,
                'ref_lon': self.lon_0, 'truelat1': self.lat_ts}

class BaseMap(object):
    EARTH_RADIUS_M = LambertConformalBasemapParameters.EARTH_RADIUS_M
    LAMBERT_CONFORMAL_PROJECTION_NAME = \
        LambertConformalBasemapParameters.PROJECTION_NAME

    def __init__(self, basemap, sn_dimension_size, we_dimension_size):
        self.__basemap = basemap
        self.__sn_dimension_size = sn_dimension_size
        self.__we_dimension_size = we_dimension_size

    @property
    def basemap(self):
        return self.__basemap

    def interp_weights(self, lat, lon):
        '''
        Determine grid position for the given lat,lon.  Returns the i and j
        of the grid point as well as the fractional part of i and j.
        '''
        i, j = self.project(lat, lon)
        i_frac = i - math.floor(i)
        j_frac = j - math.floor(j)
        i = int(math.floor(i))
        j = int(math.floor(j))
        return i, j, i_frac, j_frac

    def __xy_to_grid(self, x, y):
        i = float(x) / self.__basemap.xmax * (self.__we_dimension_size - 1)
        j = float(y) / self.__basemap.ymax * (self.__sn_dimension_size - 1)
        return i, j

    def project(self, lat, lon):
        '''
        Given lat and lon return i and j as real numbers.
        '''
        x, y = self.__basemap(lon, lat)
        if x > self.__basemap.xmax or y > self.__basemap.ymax or \
                x < self.__basemap.xmin or y < self.__basemap.ymin:
            errstr = '(%f, %f) maps to (%f, %f) ' % (lat, lon, x, y) +\
                'which is outside of grid ' +\
                '(%f - %f , %f - %f)' % (self.__basemap.xmin,
                                         self.__basemap.xmax,
                                         self.__basemap.ymin,
                                         self.__basemap.ymax)
            raise ValueError(errstr)

        # Assumes that xmin and ymin are zero
        i, j = self.__xy_to_grid(x, y)
        return i, j

    def calculate_xy(self, igrid, jgrid):
        '''
        Given a grid index return the x,y in eastings, northings in m
        '''
        x = float(igrid) / (self.__we_dimension_size - 1) * \
            (self.__basemap.xmax)
        y = float(jgrid) / (self.__sn_dimension_size - 1) * \
            (self.__basemap.ymax)
        if x > self.__basemap.xmax or y > self.__basemap.ymax or \
                x < self.__basemap.xmin or y < self.__basemap.ymin:
            raise ValueError(
                "(%d, %d) maps to (%f, %f) " % (igrid, jgrid, x, y) +
                "which is outside of grid " +
                "(%f - %f , %f - %f)" % (self.__basemap.xmin,
                                         self.__basemap.xmax,
                                         self.__basemap.ymin,
                                         self.__basemap.ymax)
            )
        return x, y

    def unproject(self, igrid, jgrid):
        '''
        Convert a grid index into a latitude and longitude.
        '''
        x, y = self.calculate_xy(igrid, jgrid)
        lon, lat = self.__basemap(x, y, inverse=True)
        return lat, lon

    def calc_north_vector(self, igrid, jgrid):
        '''
        For a given grid point calculate a unit vector pointing towards North.
        Based on the code in RIP
        '''
        x_1, y_1 = self.calculate_xy(igrid, jgrid)
        lon1, lat1 = self.__basemap(x_1, y_1, inverse=True)
        lat2 = min(89.9999, lat1 + .1)
        x_2, y_2 = self.__basemap(lon1, lat2)
        delta_u = float(x_2 - x_1)
        delta_v = float(y_2 - y_1)
        distance_north = math.sqrt(delta_u * delta_u + delta_v * delta_v)
        u_north = delta_u / distance_north
        v_north = delta_v / distance_north
        return u_north, v_north

    def get_max_i(self):
        '''
        Get the maximum allowed i index for this basemap.
        '''
        return self.__we_dimension_size - 1

    def get_max_j(self):
        '''
        Get the maximum allowed j index for this basemap.
        '''
        return self.__sn_dimension_size - 1

    def is_in_grid(self, igrid, jgrid):
        '''
        Check if the given horizontal coordinates are within the basemap.
        '''
        return igrid >= 0 and igrid <= self.__we_dimension_size and \
            jgrid >= 0 and jgrid <= self.__sn_dimension_size

    def __str__(self):
        return "WrfBaseMap[%d, %d]" % (self.__we_dimension_size,
                                       self.__sn_dimension_size)


class GeogridBaseMap(BaseMap):

    def __init__(self, ncfile, coastline_resolution=None):
        basemap = LambertConformalBasemapParameters(
            resolution=coastline_resolution).from_geogrid_file(
                ncfile).to_basemap()
        sn_dimension_size = ncfile.dimensions['south_north']
        we_dimension_size = ncfile.dimensions['west_east']
        BaseMap.__init__(self, basemap, sn_dimension_size, we_dimension_size)


class WrfBaseMap(BaseMap):

    '''
    Builds on the matplotlib BaseMap class to create a base object for a
    WRF grid.
    Based on the bmap set of functions in pycast.
    Assumes that the grid is stationary throughout the forecast run.
    '''

    def __init__(self, nc_file, longitude_variable_name,
                 latitude_variable_name, sn_dimension_name, we_dimension_name,
                 coastline_resolution=None):
        '''
        Constructor should not be used directly, use the WrfBaseMapFactory
        class instead
        '''

        for dimension_name in [sn_dimension_name, we_dimension_name]:
            if not dimension_name in nc_file.dimensions:
                raise ValueError(
                    "Input file does not have the dimension %s" % (
                        dimension_name)
                )

        basemap = LambertConformalBasemapParameters(
            resolution=coastline_resolution).from_wrf_file(
            nc_file, longitude_variable_name,
            latitude_variable_name).to_basemap()

        if NETCDF4_PRESENT and isinstance(nc_file, netCDF4.Dataset):
            sn_dimension_size = len(nc_file.dimensions[sn_dimension_name])
            we_dimension_size = len(nc_file.dimensions[we_dimension_name])
        else:
            sn_dimension_size = nc_file.dimensions[sn_dimension_name]
            we_dimension_size = nc_file.dimensions[we_dimension_name]
        BaseMap.__init__(self, basemap, sn_dimension_size, we_dimension_size)

# DH*


class Mm5BaseMap(BaseMap):

    '''
    Builds on the matplotlib BaseMap class to create a base object for a
    MM5 grid.
    Based on the bmap set of functions in pycast.
    Assumes that the grid is stationary throughout the forecast run.
    '''

    def __init__(self, nc_file, longitude_variable_name, latitude_variable_name,
                 sn_dimension_name, we_dimension_name,
                 coastline_resolution=None):
        '''
        Constructor should not be used directly, use the WrfBaseMapFactory
        class instead
        '''

        for dimension_name in [sn_dimension_name, we_dimension_name]:
            if not dimension_name in nc_file.dimensions:
                raise ValueError(
                    "Input file does not have the dimension %s" % (
                        dimension_name)
                )

        basemap = LambertConformalBasemapParameters(
            resolution=coastline_resolution).from_mm5_file(
            nc_file, longitude_variable_name,
            latitude_variable_name).to_basemap()

        sn_dimension_size = nc_file.dimensions[sn_dimension_name]
        we_dimension_size = nc_file.dimensions[we_dimension_name]
        BaseMap.__init__(self, basemap, sn_dimension_size, we_dimension_size)
# *DH


class WrfBaseMapFactory(object):

    '''
    A helper class that lazy creates basemaps from a file.
    '''

    def __init__(self, nc_file, coastline_resolution=None):
        self.__nc_file = nc_file
        self.__resolution = coastline_resolution
        self.__mass_basemap = None
        self.__u_basemap = None
        self.__v_basemap = None

    def get_u_basemap(self):
        '''
        This method is not thread safe!
        '''
        if (self.__u_basemap is None):
            self.__u_basemap = WrfBaseMap(self.__nc_file, 'XLONG_U', 'XLAT_U',
                                          'south_north', 'west_east_stag',
                                          coastline_resolution=self.__resolution)
        return self.__u_basemap

    def get_v_basemap(self):
        '''
        This method is not thread safe!
        '''
        if (self.__v_basemap is None):
            self.__v_basemap = WrfBaseMap(self.__nc_file, 'XLONG_V', 'XLAT_V',
                                          'south_north_stag', 'west_east',
                                          coastline_resolution=self.__resolution)
        return self.__v_basemap

    def get_mass_basemap(self):
        '''
        This method is not thread safe!
        '''
        if (self.__mass_basemap is None):
            self.__mass_basemap = WrfBaseMap(self.__nc_file, 'XLONG', 'XLAT',
                                             'south_north', 'west_east',
                                             coastline_resolution=self.__resolution)
        return self.__mass_basemap

    def get_nc_file(self):
        return self.__nc_file

# DH*


class Mm5BaseMapFactory(object):

    '''
    A helper class that lazy creates basemaps from a file.
    '''

    def __init__(self, nc_file, coastline_resolution=None):
        self.__nc_file = nc_file
        self.__resolution = coastline_resolution
        self.__crs_basemap = None
        self.__dot_basemap = None

    def get_dot_basemap(self):
        '''
        This method is not thread safe!
        '''
        if (self.__dot_basemap is None):
            self.__dot_basemap = Mm5BaseMap(self.__nc_file, 'longidot',
                                            'latitdot', 'i_dot', 'j_dot',
                                            coastline_resolution=self.__resolution)
        return self.__dot_basemap

    def get_crs_basemap(self):
        '''
        This method is not thread safe!
        '''
        if (self.__crs_basemap is None):
            self.__crs_basemap = Mm5BaseMap(self.__nc_file, 'longicrs',
                                            'latitcrs', 'i_cross', 'j_cross',
                                            coastline_resolution=self.__resolution)
        return self.__crs_basemap

    def get_nc_file(self):
        return self.__nc_file
# *DH

if __name__ == "__main__":
    pass
