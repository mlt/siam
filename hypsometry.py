#!/usr/bin/env python

"""Quick & Dirty way to retrieve a hypsometry given a set of starting
points and DEM.

"""

from __future__ import print_function
import logging.config
from os.path import dirname, join, exists
import numpy as np
import argparse
import multiprocessing as mp
from copy import deepcopy
from MultiProcessingLog import QueueHandler, QueueListener
import sys
import Image
from ImageDraw import Draw
from rtree import index
from collections import defaultdict

try:
    from osgeo import gdal, osr, ogr
    from osgeo.gdalconst import GA_ReadOnly
except ImportError:
    import gdal, osr, ogr
    from gdalconst import GA_ReadOnly

class Hypsometry:
    "Hypsometry creator"

    def __init__(self, args):
        "Set up DB connection"

        self.__dict__.update(args)
        self._log = logging.getLogger(__name__).getChild(self.__class__.__name__)
        self.pts_dict = dict()  # inlet FID => feature
        self.volume = defaultdict(float)    #: inlet FID => volume
        self.area = defaultdict(float)      #: inlet FID => last slice area
        self.z = []             # inlet elevations in ascending order
        p = index.Property()
        p.dimension = 3
        self.pts_idx = index.Index(properties=p)
#        gdal.SetConfigOption('PG_USE_COPY', 'YES')

    def _read(self):
        """Read in input data into memory"""

        connstr = "PG:host={:s} port={:s} dbname={:s} user={:s}".format(self.host, self.port, self.dbname, self.user)
        connstr_dem = connstr + " table={:s} mode=2".format(self.dem_table)
        self.conn_ogr = ogr.Open(connstr)
        if getattr(self, 'find_bottom', False):
            parts = self.conn_ogr.GetLayerByName(self.dem_parts)
            self.srs = parts.GetSpatialRef()
            if self.srs is None:
                self._log.warn("Failed to get spatial reference for partitions. I hope you do use the same projection with DEM!")
            else:
                self.SRID = int(self.srs.GetAttrValue("AUTHORITY", 1))
                self._log.debug("Partitions will be read from '%s' with SRID=%d", connstr, self.SRID)
            fid = parts.GetFIDColumn()
            if fid is None or fid == '':
                self._log.warning("FID column for '%s' is unknown!! I assume it is 'gid'", self.dem_parts)
                fid = 'gid'
            connstr_dem += """ where='rid in (
select rid from {:s}, {:s}
where ST_Intersects(geom, rast::geometry)
  and {:s} = {:d})'""".format(self.dem_table, self.dem_parts, fid, self.part)
        else:
            if exists(self.layer):
                self.conn_ogr = ogr.Open(self.layer)
                self.pts = self.conn_ogr.GetLayer()
                self.layer = None
            else:
                self.pts = self.conn_ogr.GetLayerByName(self.layer)
            self.srs = self.pts.GetSpatialRef()
            if self.srs is None:
                self._log.warn("Failed to get spatial reference for inlets. I hope you do use the same projection with DEM!")
            else:
                self.SRID = int(self.srs.GetAttrValue("AUTHORITY", 1))
                self._log.info("Points will be read from '%s' with SRID=%d", connstr, self.SRID)
            if getattr(self, 'part', False):
                fid = self.pts.GetFIDColumn()
                if fid is None or fid == '':
                    self._log.critical("FID column is unknown!! I assume it is 'gid'")
                    fid = 'gid'
                where = "{:s} in (select gid from {:s} where pid={:d})".format(fid, self.parts_map, self.part)
                # it is partition table that shall be created with
                # necessary points only if those need to be restricted
                if self.pts.SetAttributeFilter(where):
                    self._log.critical("Failed to restrict features for partition #%d", self.part)
                connstr_dem += """ where='rid in (
select rid from {:s}, {:s}
where ST_Contains(geom, rast::geometry)
  and pid = {:d})'""".format(self.dem_table, self.dem_parts, self.part)

        dataset = gdal.Open(connstr_dem, GA_ReadOnly )
        if dataset is None:
            self._log.critical("Can't open %s", connstr_dem)
            raise Exception()
        self.geotransform = dataset.GetGeoTransform()
        self.invgeotransform = gdal.InvGeoTransform(self.geotransform)[1]
        self.proj = dataset.GetProjectionRef()
        srs = osr.SpatialReference(self.proj)
        SRID = int(srs.GetAttrValue("AUTHORITY", 1))
        band = dataset.GetRasterBand(1)
        self.NODATA = band.GetNoDataValue()
        self._log.debug("Raster from '%s' is used with SRID=%d and NoDATA=%f", connstr_dem, SRID, self.NODATA)
        if self.SRID != SRID:
            raise Exception("Points & DEM should use same coordinate system")
        self.raster = band.ReadAsArray()

    def _find_minmax(self):
        "An attempt to find range for partition to slice"
        masked = np.ma.masked_array(self.raster, mask=np.isclose(self.raster, self.NODATA))
        zmin = masked.min()
        self.max_height = sys.maxint
        self.z.append(zmin)

    @staticmethod
    def _transform(transform, x, y):
        xx = transform[0] + transform[1] * x + transform[2] * y
        yy = transform[3] + transform[4] * x + transform[5] * y
        return xx, yy

    def _world2pixel(self, x, y):
        return self._transform(self.invgeotransform, x, y)

    def _pixel2world(self, x, y):
        return self._transform(self.geotransform, x, y)

    def raster_value(self, geom):
        "Read raster cell value at point"
        iPixel, iLine = self._world2pixel(geom.GetX(), geom.GetY())
        return self.raster[iLine, iPixel]

    def _find_pixels(self):
        """Read raster values for POIs.

Perhaps would be better to use SQL but this way we hopefully can use a
 shapefile for non --mp version.

        """
        self.pts.ResetReading()
        for feat in self.pts:
            geom = feat.GetGeometryRef()
            x = geom.GetX()
            y = geom.GetY()
            fid = feat.GetFID()
            self._log.debug("%d => %.3f, %.3f", fid, x, y)

            try:
                val = self.raster_value(geom)
                if val != self.NODATA: # GDAL 1.9 is broken also should not(?) compare floats
                    pt = ogr.Geometry(ogr.wkbPoint25D)
                    pt.SetPoint(0, x, y, float(val))
                    self.pts_dict[fid] = pt
                    self.pts_idx.insert(fid, (x, y, val, x, y, val))
                    self.z.append(val)
            except IndexError:
                pass
        self.z.sort()

    def _mk_bottom_lyr(self):
        """Find bottom and append it to existing table created by Starter"""

        # self._log.debug('Finding lowest point in partition {:d}'.format(self.part))
        lyr = self.conn_ogr.ExecuteSQL("""
insert into {layer:s} (z, geom, pid)
select (gv).val, ST_Force_3D(ST_Centroid((gv).geom)) geom, {pid:d}
from (
  select ST_DumpAsPolygons(ST_Union(ST_Clip(rast, geom))) gv
  from {dem:s}, {poly:s}
  where ST_Intersects(rast, geom)
        and gid={pid:d}
  group by gid
) foo
where (gv).val!=(select ST_BandNoDataValue(rast,1) from {dem:s} where rid=1)
order by (gv).val
limit 1
        returning gid, z, geom;""".format(
            layer=self.layer, pid=self.part, dem=self.dem_table, poly=self.dem_parts))
        assert 1 == lyr.GetFeatureCount()
        feat = lyr.GetNextFeature()
        geom = feat.GetGeometryRef()
        pt = ogr.Geometry(ogr.wkbPoint25D)
        x = geom.GetX()
        y = geom.GetY()
        z = feat.GetField('z')
        pt.SetPoint(0, x, y, z)
        self.z.append(z)
        self.pts_dict[self.part] = pt
        self.pts_idx.insert(self.part, (x, y, z, x, y, z))
        self.conn_ogr.ReleaseResultSet(lyr)

    def _mkmem(self):
        """Set up in-memory stuff"""

        self._log.debug('Creating in-memory stuff')
        mrd = gdal.GetDriverByName('MEM')
        self.inmem = mrd.Create('', self.raster.shape[1], self.raster.shape[0], 2, gdal.GDT_Byte)
        if gdal.CE_Failure == self.inmem.SetProjection(self.proj):
            self._log.critical('Failed to set projection for in-memory raster')
        self.inmem.SetGeoTransform(self.geotransform)
        self.imb = self.inmem.GetRasterBand(1)
        self.mask = self.inmem.GetRasterBand(2)
        self.mask.WriteArray(np.logical_not(np.isclose(self.raster, self.NODATA)))
        drv = ogr.GetDriverByName('Memory')
        self.dst_ds = drv.CreateDataSource('out') # is it a must?
        self.boundary = self.dst_ds.CreateLayer('boundary', self.srs, ogr.wkbLineString)
        bndry = self.dst_ds.CreateLayer('bndry', self.srs, ogr.wkbPolygon)
        fd = ogr.FieldDefn('valid', ogr.OFTInteger)
        bndry.CreateField(fd)
        gdal.Polygonize(self.mask, None, bndry, 0)
        if bndry.SetAttributeFilter('valid > 0'):
            self._log.critical("Failed to restrict valid data boundary polygons")
        for b in bndry:
            g = b.GetGeometryRef()
            f = ogr.Feature(self.boundary.GetLayerDefn())
            f.SetGeometry(g.Boundary())
            self.boundary.CreateFeature(f)


    def _getout(self):
        """Connect to shared output layer"""
        self.out_ogr = ogr.Open(self.out, True)
        self.polys = self.out_ogr.GetLayerByName(self.table)
        if self.polys is None:
            self._log.critical('Failed to get an output layer %s', self.table)
        if getattr(self, 'find_bottom', False):
            self.pts = self.out_ogr.GetLayerByName(self.layer)
            if self.pts is None:
                self._log.critical('Failed to open bottoms layer')

    def _process(self, z):
        """Elevation steps for a 'given' inlet and possibly others"""
        accepted = 1
        while accepted>0:
            tmp = self.raster <= z
            self.imb.WriteArray(tmp.astype(np.byte))
            self._log.debug('Polygonizing...')
            lyr = self.dst_ds.CreateLayer('polys', self.srs, ogr.wkbPolygon)
            fd = ogr.FieldDefn('below', ogr.OFTInteger)
            lyr.CreateField( fd )
            gdal.Polygonize(self.imb, self.mask, lyr, 0)
            if lyr.SetAttributeFilter('below > 0'):
                self._log.critical("Failed to restrict polygons")
            accepted = sum(self._consider_polygon(z, polygon) for polygon in lyr if not self._remove_contaminated(z, polygon))
            self._log.debug('Looked through %d polygons @ %.2f, accepted %d', lyr.GetFeatureCount(), z, accepted)
#             if not accepted:
#                 self._log.debug('All polygons have been eliminated! Bailing out')
            self.dst_ds.DeleteLayer(1)
            z += self.step
        return z

    def _process_queue(self):
        """Main loop to deal with elevations queue"""
        for z_inlet in self.z:
            self._log.debug("Moving to next elevation in queue %.2f", z_inlet)
            self._process(z_inlet)


    def _remove_contaminated(self, z, feat):
        "Test whether polygon is contaminated and remove points within"
        polygon = feat.GetGeometryRef()
        self.boundary.ResetReading()
        for feature in self.boundary:
            boundary = feature.GetGeometryRef()
            if boundary.Intersects(polygon):
                self._log.debug('Reached boundary. Skipping.')
                env = polygon.GetEnvelope()
                pre = self.pts_idx.intersection((env[0], env[2], -sys.maxint, env[1], env[3], z), True)
                within = (item for item in pre if self.pts_dict[item.id].Within(polygon))
#                 idx_in = self.pts_idx.count((env[0], env[2], -sys.maxint, env[1], env[3], z))
#                 cnt = 0
                for item in within:
                    self.pts_idx.delete(item.id, item.bbox)
                    del self.pts_dict[item.id]
#                     cnt += 1
#                 idx_out = self.pts_idx.count((env[0], env[2], -sys.maxint, env[1], env[3], z))
#                 assert cnt == (idx_in - idx_out)
                return True

        return False

    def _add_minimum_bottom(self, polygon):
        """Accurately find the lowest location within a polygon.

        Quite slow. Use when having a large step to avoid centroids for oxbow shaped depressions.

        :param Geometry polygon: The polygon to use as a mask when looking for minimum
        :return: newly created OGR feature id
        :rtype: int
        :raises AssertionError: if stumbled upon NODATA. Supplied polygon is expected to be derived previously from valid data.
        """
        points = polygon.GetGeometryRef(0)
        rng = xrange(points.GetPointCount())
#         pixels = np.vstack(self._world2pixel(
#                                    np.vectorize(points.GetX, otypes=[np.float])(rng),
#                                    np.vectorize(points.GetY, otypes=[np.float])(rng))).transpose().flatten().tolist() # 67.12 sec
#         pixels = np.vstack(self._world2pixel(
#                                    np.vectorize(points.GetX)(rng),
#                                    np.vectorize(points.GetY)(rng))).transpose().flatten().tolist() # 66.8 sec
#         xx, yy = self._world2pixel(
#                                    np.vectorize(points.GetX)(rng),
#                                    np.vectorize(points.GetY)(rng))
#         pixels = zip(xx, yy) # 67.26 sec
#         pixels = [self._world2pixel(points.GetX(p), points.GetY(p)) for p in range(points.GetPointCount())] # 65.87 sec
        pixels = [self._world2pixel(points.GetX(p), points.GetY(p)) for p in rng] # 66.44 sec
        img = Image.new('L', self.raster.shape[::-1], 1)
        Draw(img).polygon(pixels)
        mask = np.fromstring(img.tostring(), 'b')
        mask.shape = self.raster.shape
        masked = np.ma.masked_array(self.raster, mask=mask)
        y, x = np.unravel_index(np.argmin(masked), mask.shape)
        z = self.raster[y, x]
        assert not np.isclose(z, self.NODATA)
        pt = ogr.Geometry(ogr.wkbPoint25D)
        easting, northing = self._pixel2world(x, y)
        pt.SetPoint(0, easting, northing, float(z))
        return self._add_bottom_point(pt)

    def _add_bottom_point(self, geom):
        """
        Add point feature (bottom) to the output layer

        :param float z: elevation of the bottom
        :param Geometry geom: bottom point :gdal:`geometry <osgeo.ogr.Geometry-class>`
        :return: feature id added
        :rtype: int
        """
        f = ogr.Feature(self.pts.GetLayerDefn())
        z = geom.GetZ()
        f.SetField('z', z)
        f.SetField('pid', self.part)
        f.SetGeometry(geom)
        self.pts.CreateFeature(f)
        k = f.GetFID()
        self._log.debug('Created %d', k)
        self.pts_dict[k] = geom
        x = geom.GetX()
        y = geom.GetY()
        self.pts_idx.insert(k, (x, y, z, x, y, z))
        return k

    def _add_centroid_bottom(self, polygon):
        """
        Add the centroid of the polygon as a bottom

        :param Geometry polygon: currently being processed polygon
        :return: feature id added
        :rtype: int
        """
        centroid = polygon.Centroid()
        z = self.raster_value(centroid)
        pt = ogr.Geometry(ogr.wkbPoint25D)
        pt.SetPoint(0, centroid.GetX(), centroid.GetY(), float(z))
        return self._add_bottom_point(pt)

    def _consider_polygon(self, z, feat):
        """Find polygons for inlets"""
        polygon = feat.GetGeometryRef()

        # Find all points within, select first lowest, link & delete others
        # TODO: there should be no point elevation check as we discover points on the go
        # This is valid however only for auto discovery but not existing starting points
        env = polygon.GetEnvelope()
        pre = self.pts_idx.intersection((env[0], env[2], -sys.maxint, env[1], env[3], z), True)
        within = [item for item in pre if self.pts_dict[item.id].Within(polygon)]

        if len(within):
            lowest = min(enumerate(within), key=lambda tup: self.pts_dict[tup[1].id].GetZ())[0]
            k = within.pop(lowest).id
            if len(within):
                gids = ','.join(str(k.id) for k in within)
                self.out_ogr.ExecuteSQL('update {bottoms:s} set merge_to={lowest:d} where gid in ({gids}) and pid={part:d}'.format(bottoms=self.layer, lowest=k, gids=gids, part=self.part))
                self._log.debug('Merging %s points into %d', gids, k)
#                 idx_in = self.pts_idx.count((env[0], env[2], -sys.maxint, env[1], env[3], z))
#                 cnt = 0
                for item in within:
                    self.pts_idx.delete(item.id, item.bbox)
                    del self.pts_dict[item.id]
#                     cnt += 1
#                 idx_out = self.pts_idx.count((env[0], env[2], -sys.maxint, env[1], env[3], z))
#                 assert cnt == (idx_in - idx_out)
        else:
            try:
                method = {'fast': self._add_centroid_bottom,
                          'rigorous': self._add_minimum_bottom}[self.find_bottom]
            except (KeyError, AttributeError):
                # No points are within and we can't add any so let's jump to next inlet above
                return False
            k = method(polygon)

        zmin = min(z, self.pts_dict[k].GetZ()) # odd small negatives

        self._log.debug('Found polygon for point %d at %.2f', k, z)
        f = ogr.Feature(self.polys.GetLayerDefn())
        f.SetField('polygon', feat.GetFID())
        f.SetField('point', k)
        # Somehow without float() it says NotImplementedError: Wrong number of arguments for overloaded function 'Feature_SetField'.
        f.SetField('z', float(z))
        stage = float(z - zmin)
        f.SetField('stage', stage)
        f.SetGeometry(polygon)
        area = polygon.GetArea()
        f.SetField('area', area)
        self.volume[k] += min(stage, self.step) * .5*(area + self.area[k])
        self.area[k] = area
        f.SetField('volume', self.volume[k])
        self.polys.CreateFeature(f)
        return True

    def run(self):
        self._read()
        if not getattr(self, 'find_bottom', False):
            self._find_pixels()
        else:
            if 'single' == self.find_bottom:
                self._mk_bottom_lyr()
            else:
                self._find_minmax()
        self._mkmem()
        self._getout()
        self._process_queue()

    # def __del__(self):
    #     "Disconnect from DB"
    #     self._log.debug("Closing ODBC/PG connection")
    #     self.conn.close()

def profile_part(args):
    import cProfile
    cProfile.runctx('run_part(args)', globals(), locals(), '%s_%d.prof' % (args['profile'], args['part'] or 0))

def run_part(args):
    """Entry point for sub-processes.

    Must be plain function on Windows. It can't be un-/bound function:("""
    _log = logging.getLogger(__name__)
    _log.setLevel(args.get('_loglevel', logging.NOTSET))
    hndlr = QueueHandler(args['_queue'])
    _log.addHandler(hndlr)
    _log.info("Working on part %d", args['part'])
    h = Hypsometry(args)
    try:
        h.run()
    except:
        # This will funnel formatted exceptions via mp.Queue
        _log.exception('Failed to run hypsometry')
    _log.info("Done working on part %d", args['part'])
    _log.removeHandler(hndlr)

class Starter:
    """Performs initial work before spawning workers.

    .. todo:: Consider OGR_TRUNCATE instead of OVERWRITE

    """

    def __init__(self, args):
        self.__dict__.update(args)
        self.args = args        # to pass around lean copy
        self._log = logging.getLogger(__name__).getChild(self.__class__.__name__)
        self._log.info("Using GDAL '%s'", gdal.VersionInfo("RELEASE_NAME"))

        manager = mp.Manager()
        self.queue = manager.Queue(-1)
        listener = QueueListener(self.queue)
        self.pool = None

    def mkout(self):
        """Set up output stuff"""

        self._log.info("Making output table '%s'", self.table)
        self.conn_ogr = ogr.Open("PG:host={:s} port={:s} dbname={:s} user={:s}".format(self.host, self.port, self.dbname, self.user))
        if exists(self.layer):
            conn = ogr.Open(self.layer)
            lyr = conn.GetLayer()
            self.layer = None
        else:
            if getattr(self, 'find_bottom', False):
                lyr = self.conn_ogr.GetLayerByName(self.dem_parts)
            else:
                lyr = self.conn_ogr.GetLayerByName(self.layer)
        self.srs = lyr.GetSpatialRef()

        self.out_ogr = ogr.Open(self.out, True)
        srid = int(self.srs.GetAttrValue("AUTHORITY", 1))
        if getattr(self, 'unlogged', False):
            sql = """drop table if exists {table};
create unlogged table {table} (
  gid serial primary key,
  geom geometry(PolygonZ, {srid}) not null,
  z real not null,
  polygon integer not null,
  point integer not null,
  stage real not null,
  area real not null,
  volume real not null);""".format(table=self.table, srid=srid)
            self.out_ogr.ExecuteSQL(sql)
        else:
            output_lyr = self.out_ogr.CreateLayer(self.table, self.srs, ogr.wkbPolygon, ['OVERWRITE=YES','GEOMETRY_NAME=geom','FID=gid','COLUMN_TYPES=z=real,stage=real,area=real,volume=real'])
            if output_lyr is None:
                self._log.critical('Failed to create an output layer %s', self.table)
            fd = ogr.FieldDefn('z', ogr.OFTReal)
            output_lyr.CreateField(fd)
            # polygon id (for a given z) + z make globally unique polygon id
            fd = ogr.FieldDefn('polygon', ogr.OFTInteger)
            output_lyr.CreateField(fd)
            fd = ogr.FieldDefn('point', ogr.OFTInteger)
            output_lyr.CreateField(fd)
            fd = ogr.FieldDefn('stage', ogr.OFTReal)
            output_lyr.CreateField(fd)
            fd = ogr.FieldDefn('area', ogr.OFTReal)
            output_lyr.CreateField(fd)
            fd = ogr.FieldDefn('volume', ogr.OFTReal)
            output_lyr.CreateField(fd)

        if getattr(self, 'find_bottom', False):
            self._log.info("Creating layer for POIs in '%s'", self.layer)
            if getattr(self, 'unlogged', False):
                sql = """drop table if exists {layer};
create unlogged table {layer} (
  gid serial primary key,
  geom geometry(PointZ, {srid}) not null,
  z real not null, {merge}
  pid integer not null);""".format(layer=self.layer, merge={False: '', True: 'merge_to integer,'}['single' != self.find_bottom], srid=srid)
                self.out_ogr.ExecuteSQL(sql)
            else:
                pts_lyr = self.conn_ogr.CreateLayer(self.layer, self.srs, ogr.wkbPoint25D, ['OVERWRITE=YES','GEOMETRY_NAME=geom','FID=gid','COLUMN_TYPES=z=real'])
                if pts_lyr is None:
                    self._log.critical('Failed to create a point layer %s', self.layer)
                fd = ogr.FieldDefn('z', ogr.OFTReal)
                pts_lyr.CreateField(fd)
                if 'single' != self.find_bottom:
                    fd = ogr.FieldDefn('merge_to', ogr.OFTInteger)
                    pts_lyr.CreateField(fd)
                # pid is somewhat useless as gid numbering is global
                fd = ogr.FieldDefn('pid', ogr.OFTInteger)
                pts_lyr.CreateField(fd)


    def _prepare(self):
        """Set up partitions for parallel processing"""

        if getattr(self, 'find_bottom', False):
            self._log.info("Counting existing partitions in '%s'", self.dem_parts)
            lyr = self.conn_ogr.GetLayerByName(self.dem_parts)
        else:
            self._log.info("Creating partitions in '%s'", self.dem_parts)
            lyr = self.conn_ogr.CreateLayer(self.dem_parts, self.srs, ogr.wkbPolygon,
                                              ['OVERWRITE=YES','GEOMETRY_NAME=geom','FID=pid','PG_USE_COPY=YES'])
            if lyr is None:
                self._log.critical("Failed to create a layer '%s' for partitions", self.dem_parts)
            # Indices will be implicitly created for this layer
            where = '' if self.where is None else 'and ' + self.where
            lyr = self.conn_ogr.ExecuteSQL("""
insert into {dem_parts:s}(pid, geom)
select row_number() over () as pid, (g).geom as geom
from (
  select ST_Dump(ST_Union(rast::geometry)) as g
  from {side_inlets:s}, {dem_table:s}
    where ST_DWithin(rast::geometry, geom, {radius:f}) {where:s}
      and not ST_BandIsNoData(rast, false)
) foo returning pid;""".format(dem_parts=self.dem_parts, side_inlets=self.layer,
                               dem_table=self.dem_table, radius=self.radius, where=where))

            self._log.info("Establishing mapping between inlets and partitions")
            where = '' if self.where is None else 'where ' + self.where
            map_lyr = self.conn_ogr.ExecuteSQL("""
drop table if exists {side_inlets_parts:s};
create table {side_inlets_parts:s} (
  pid int not null,
  gid int primary key
);
insert into {side_inlets_parts:s}(gid, pid)
select gid, pid
from (
  select gid, pid, row_number() over inlet
  from (
    select variance(ST_Value(rast, i.geom)) over part, gid, p.pid
    from {dem_parts:s} p
        join {side_inlets:s} i on ST_Dwithin(i.geom, p.geom, {radius:f})
    join {dem_table:s} on ST_Intersects(i.geom, rast) {where:s}
    window part as (partition by p.pid)
  ) foo
  window inlet as (partition by gid order by variance)
) bar
where row_number = 1;
create index on {side_inlets_parts:s}(pid);
            """.format(side_inlets_parts=self.parts_map, side_inlets=self.layer,
                       dem_parts=self.dem_parts, dem_table=self.dem_table,
                       radius=self.radius, where=where))
            # self.conn_ogr.ReleaseResultSet(map_lyr)

        lyr.ResetReading()
        feat = lyr.GetNextFeature()
        out = []
        while not feat is None:
            out.append(feat.GetFID())
            feat = lyr.GetNextFeature()

        cnt = len(out)
        self._log.debug("Having %d partitions", cnt)
        # if not self.find_bottom:
        #     self.conn_ogr.ReleaseResultSet(lyr)
        return out

    def run(self):
        if getattr(self, 'fixup', False):
            # We have to keep output table in same database to be able
            # to selectively remove records. Alternatively we can
            # fetch list of gid using WHERE and delete those in output
            # table using point field.
            pass
        else:
            self.mkout()
        if getattr(self, 'mp', False):
            parts = self._prepare()
            self.pool = mp.Pool(processes = min(len(parts), self.threads)) # is it that bad to have extra dormant workers??
            if getattr(self, 'profile', False):
                entry = profile_part
            else:
                entry = run_part
            # args = vars(args)
            def fix_part(part):
                out = deepcopy(self.args)
                out['part'] = part
                out['_queue'] = self.queue
                return out
            self.pool.map(entry, [fix_part(x) for x in parts])
        else:
            h = Hypsometry(self.args)
            h.run()
        self.add_indices()

    def add_indices(self):
        self._log.info('Building indexes')
        self.out_ogr.ExecuteSQL('create index on "%s" (point)' % self.table)
        if getattr(self, 'unlogged', False):
            self.out_ogr.ExecuteSQL('create index on "%s" using gist (geom)' % self.table)
            if getattr(self, 'find_bottom', False):
                self.out_ogr.ExecuteSQL('create index on "%s" using gist (geom)' % self.layer)

    def kill(self):
        """To be used from GUI thread to abort operations"""
        if not self.pool is None:
            self.pool.terminate()
            self.pool.join()
            # self.pool = None

if __name__ == '__main__':
    if (exists('logging.conf')):
        logging.config.fileConfig('logging.conf')
    else:
        logging.config.fileConfig(join(dirname(__file__), 'logging.conf'))
    _log = logging.getLogger(__name__)

    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--host',
                        default='localhost',
                        help='PGHOST')
    parser.add_argument('--port',
                        default='5432',
                        help='PGPORT')
    parser.add_argument('--dbname',
                        default='gis',
                        help='Database name')
    parser.add_argument('--user',
                        default='user',
                        help='User name to use while connecting to database. Use dot pgpass to supply a password.')
    parser.add_argument('--dem-table',
                        default='dem50',
                        help='Raster table name with DEM')
    parser.add_argument('--dem-parts',
                        default='dem_parts',
                        help='Table to be created with partition polygons')
    parser.add_argument('--parts-map',
                        default='side_inlets_parts',
                        help='Table with inlet-partition mapping. Used only without --find-bottom')
    parser.add_argument('--layer',
                        default='side_inlets',
                        help='Layer name with points of interest')
    parser.add_argument('--find-bottom', choices=['single', 'fast', 'rigorous'],
                        help='Create LAYER (will overwrite!) with POIs in the lowest: single) points of existing partitions, fast) depression polygon centroids, rigorous) point of depressions')
    parser.add_argument('--out',
                        default='PG:host=localhost port=5432 dbname=gis user=user',
                        help='An output recognizeable by OGR')
    parser.add_argument('--table', default='hypsometry',
                        help='Table name in DB defined by OUT. Existing table if any will be dropped!')
    parser.add_argument('--step', type=float,
                        default=.01,
                        help='Stage step')
    parser.add_argument('--radius', type=float,
                        default=500,
                        help='Clip/search radius for partitions')
    parser.add_argument('--part', type=int,
                        help='Specific partition to process. Used internally')
    parser.add_argument('--where',
                        help="SQL restriction on side inlets to consider, e.g., 'gid in (1,2)'")
    parser.add_argument('--fixup', action='store_true',
                        help="Don't overwrite output table but delete records for points matching WHERE")
    parser.add_argument('--threads', type=int,
                        default=mp.cpu_count(),
                        help='Maximum number of parallel processes')
    parser.add_argument('--mp', action='store_true',
                        help='Use multiprocess & partitioning. This is required if --find-bottom is used.')
    parser.add_argument('--unlogged', action='store_true',
                        help='Use unlogged PostgreSQL tables')
    parser.add_argument('--profile',
                        help='Prefix for cProfile output')
    # parser.add_argument('--verbose', action='store_true',
    #                     help='Drop logging level to DEBUG')
    args = parser.parse_args()
    args = vars(args)
    args['_loglevel'] = _log.getEffectiveLevel()
    s = Starter(args)
    s.run()
    _log.info("Exiting")
