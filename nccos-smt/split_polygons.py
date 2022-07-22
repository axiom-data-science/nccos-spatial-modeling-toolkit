from typing import Union, List, Generator
import json
import copy
import math
from shapely import ops, affinity
from shapely.geometry import LineString, Polygon, MultiPolygon, GeometryCollection, box, mapping

def check_crossing(lon1: float, lon2: float,
                   validate: bool = False, dlon_threshold: float = 180.0):
    """
    Assuming a minimum travel distance between two provided longitude coordinates,
    checks if the 180th meridian (antimeridian) is crossed in either direction.
    """

    if validate and any([abs(x) > 180.0 for x in [lon1, lon2]]):
        raise ValueError("longitudes must be in degrees [-180.0, 180.0]")   
    return abs(lon2 - lon1) > dlon_threshold


def translate_polygons(geometry_collection: GeometryCollection) -> Generator[
                          Union[List[dict], List[Polygon]], None, None
                       ]:
    """Shift polygons to the other side of the range if their bounds are outside -180, 180."""
    for polygon in geometry_collection.geoms:
        (minx, _, maxx, _) = polygon.bounds
        if minx < -180:
            geo_polygon = affinity.translate(polygon, xoff = 360)
        elif maxx > 180:
            geo_polygon = affinity.translate(polygon, xoff = -360)
        else: geo_polygon = polygon

        yield geo_polygon
        
            
def split_polygon(_polygon: Polygon, validate: bool = False) -> Union[
    List[dict], List[Polygon], GeometryCollection
    ]:
    """
    Given a shapely representation of a Polygon, returns a collection of
    'antimeridian-safe' constituent Polygons split at the 180th meridian, 
    ensuring compliance with GeoJSON standards (https://tools.ietf.org/html/rfc7946#section-3.1.9)
    Assumptions:
      - Any two consecutive points with over 180 degrees difference in
        longitude are assumed to cross the antimeridian
      - The polygon spans less than 360 degrees in longitude (i.e. does not wrap around the globe)
      - However, the polygon may cross the antimeridian on multiple occasions
    
    Parameters:
        geojson (dict): GeoJSON of input polygon to be split. For example:
                        {
                        "type": "Polygon",
                        "coordinates": [
                          [
                            [179.0, 0.0], [-179.0, 0.0], [-179.0, 1.0],
                            [179.0, 1.0], [179.0, 0.0]
                          ]
                        ]
                        }
        output_format (str): Available options: "geojson", "polygons", "geometrycollection"
                             If "geometrycollection" returns a Shapely GeometryCollection.
                             Otherwise, returns a list of either GeoJSONs or Shapely Polygons
        validate (bool): Checks if all longitudes are within [-180.0, 180.0]
      
    Returns:
        List[Polygon]: antimeridian-safe polygon(s)
    """
    
    geojson = json.loads(json.dumps(mapping(_polygon))) # convert to mutable geojson
    coords_shift = copy.deepcopy(geojson["coordinates"])
    shell_minx = shell_maxx = None
    split_meridians = set()
    splitter = None

    for ring_index, ring in enumerate(coords_shift):    
        if len(ring) < 1: 
            continue
        else:
            ring_minx = ring_maxx = ring[0][0]
            crossings = 0
            
        for coord_index, (lon, _) in enumerate(ring[1:], start=1):
            lon_prev = ring[coord_index - 1][0] # [0] corresponds to longitude coordinate
            if check_crossing(lon, lon_prev, validate=validate):
                direction = math.copysign(1, lon - lon_prev)
                coords_shift[ring_index][coord_index][0] = lon - (direction * 360.0)
                crossings += 1

            x_shift = coords_shift[ring_index][coord_index][0]
            if x_shift < ring_minx: ring_minx = x_shift
            if x_shift > ring_maxx: ring_maxx = x_shift

        # Ensure that any holes remain contained within the (translated) outer shell
        if (ring_index == 0): # by GeoJSON definition, first ring is the outer shell
            shell_minx, shell_maxx = (ring_minx, ring_maxx)
        elif (ring_minx < shell_minx):
            ring_shift = [[x + 360, y] for (x, y) in coords_shift[ring_index]]
            coords_shift[ring_index] = ring_shift
            ring_minx, ring_maxx = (x + 360 for x in (ring_minx, ring_maxx))
        elif (ring_maxx > shell_maxx):
            ring_shift = [[x - 360, y] for (x, y) in coords_shift[ring_index]]
            coords_shift[ring_index] = ring_shift
            ring_minx, ring_maxx = (x - 360 for x in (ring_minx, ring_maxx))

        if crossings: # keep track of meridians to split on
            if ring_minx < -180:
                split_meridians.add(-180)
            if ring_maxx > 180:
                split_meridians.add(180)

    n_splits = len(split_meridians)
    if n_splits > 1:
        raise NotImplementedError(
            """Splitting a Polygon by multiple meridians (MultiLineString) 
               not supported by Shapely"""
        )
    elif n_splits == 1:
        split_lon = next(iter(split_meridians))
        meridian = [[split_lon, -90.0], [split_lon, 90.0]]
        splitter = LineString(meridian)

    shell, *holes = coords_shift if splitter else geojson["coordinates"]
    if splitter:
        split_polygons = ops.split(Polygon(shell, holes), splitter)
    else:
        split_polygons = GeometryCollection([Polygon(shell, holes)])
        
    geo_polygons = list(translate_polygons(split_polygons))
    
    return MultiPolygon(geo_polygons)


def split_multipolygon(_multipolygon: MultiPolygon):
    "Use the split_polygon function to split up a shapely Multipolygon"
    
    split_result = []
    for _polygon in _multipolygon.geoms:
        result = split_polygon(_polygon)
        split_result.extend(result.geoms)
        
    return MultiPolygon(split_result)


def wrap_polygons(geometry_collection, verbose=False):
    """
    For a geometry collection that has mixed-sign components, try wrapping the negative
    elements to the positive side. If it makes for a smaller bounding box, keep it.
    """
    
    # Test initial bounds
    (minx, _, maxx, _) = geometry_collection.bounds
    initial_difference = maxx - minx    
    transformed = []
    
    for _polygon in geometry_collection.geoms:
        (minx, _, maxx, _) = _polygon.bounds
        if minx < 0:
            transformed.append(affinity.translate(_polygon, xoff=360))
        else:
            transformed.append(_polygon)
    
    transformed = MultiPolygon(transformed)
    (minx, _, maxx, _) = transformed.bounds
    transformed_difference = maxx - minx
    
    if transformed_difference < initial_difference:
        if verbose: print("Wrapped polygons to reduce footprint")
        return transformed
    else:
        return geometry_collection
    
def splitter(geometry_collection, verbose=False):
    
    if type(geometry_collection) == MultiPolygon:
        _split = split_multipolygon(geometry_collection)
    if type(geometry_collection) == Polygon:
        _split = split_polygon(geometry_collection)
        
    result = wrap_polygons(_split, verbose=verbose)
    
    return result