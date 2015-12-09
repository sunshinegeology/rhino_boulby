import Rhino as rh
import rhinoscriptsyntax as rs
import scriptcontext as sc
import System.Guid as sg
import copy as cp
import math as ma
import random as rd
import sys
import os

"""
TODO
- variations: skew, orientation
"""

def unif_rand(lo, hi):
    return rd.random()*(hi-lo)+lo


def rand_point(xlim, ylim):
    """Random points at nil z."""
    return rh.Geometry.Point3d(unif_rand(xlim[0], xlim[1]), unif_rand(ylim[0], ylim[1]), 0)


class intersection_data:
    """Intersection data class containing the start&endpoint, and the highest point along the line"""
    def __init__(self, spt, ept, hipt):
        self.spt = spt
        self.ept = ept
        self.hipt = hipt


def width(isdata):
    widthv = rs.VectorCreate(isdata.spt, isdata.ept)
    return rs.VectorLength(widthv)


def height(isdata):
    """Assumes end and start point have the same z-value."""
    return isdata.hipt[2]-isdata.ept[2]


def aspect_ratio(isdata):
    """Assumes end and start point have the same z-value, returns height:width."""
    return height(isdata)/width(isdata)


def output_list(fname, values):
    olines = [''.join([str(v), '\n']) for v in values]
    with open(fname, 'w') as f:
        f.writelines(olines)


def skewness(isdata):
    """Returns skewness as a measure of asymmetry as (on xy-plane): hipt-cpt/width."""
    cpt_xy = (isdata.spt+isdata.ept)/2.
    cpt_xy[2] = 0.
    hipt_xy = cp.deepcopy(isdata.hipt)
    hipt_xy[2] = 0.
    skew = rs.VectorCreate(hipt_xy, cpt_xy)
    return rs.VectorLength(skew)/width(isdata)


def analyze_intersection(cidx):
    """Returns intersection_data for curve."""
    ept = rs.CurveEndPoint(cidx)
    spt = rs.CurveStartPoint(cidx)
    hipt = rh.Geometry.Point3d(0.,0.,-sys.float_info.max)
    cpts = rs.CurvePoints(cidx)
    for cpt in cpts:
        if cpt[2] > hipt[2]:
            hipt = cpt
    return intersection_data(spt, ept, hipt)


rmin, rmax = 3.0, 5.5
mean_radius = (rmax+rmin)/2
sigma_radius = mean_radius/2.
def feature_radius(mean_radius):
    radius = 0.
    while radius < rmin or radius > rmax:
        radius = rd.gauss(mean_radius, sigma_radius)
    return radius


sigma_rotation = 1.
def feature_orientation(mean_orientation):
    return mean_orientation + rd.gauss(mean_orientation, sigma_rotation)


feature_center_min_dist = rmax*2
def non_intersecting_center_points(pt_no, xlim, ylim):
    pts = [rand_point(xlim, ylim)]
    while len(pts) < pt_no:
        new_pt = rand_point(xlim, ylim)
        cl_idx = rs.PointArrayClosestPoint(pts, new_pt)
        cl_vct = rs.VectorCreate(new_pt, pts[cl_idx])
        if rs.VectorLength(cl_vct) > feature_center_min_dist:
            pts += [new_pt]
    return pts


def delete_zero_width_elements(objects):
    for i in objects:
        if width(analyze_intersection(i)) == 0.:
            rs.DeleteObject(i)


def update_views():
    """
    views = rs.ViewNames()
    for view in views:
        rs.ViewDisplayMode(view, 'Ghosted')
    """
    sc.doc.Views.Redraw()


def main():
    # resetting
    rs.DocumentModified(False)
    rs.Command('_-New _None')
    update_views()
    
    # settings    
    half_length = 100.
    xlim, ylim = (-half_length,half_length), (-half_length,half_length)
    feat_no = 50
    shear_angle = 30.
    
    # adding layers
    feat_lname = 'domes'
    l = rs.AddLayer(feat_lname)
    rs.CurrentLayer(feat_lname)
    
    # adding features
    feat_ids = []
    feat_radii = [feature_radius(mean_radius) for i in range(feat_no)]
    rs.CurrentLayer(feat_lname)
    feat_centers = non_intersecting_center_points(feat_no, xlim, ylim)
    for i, c in enumerate(feat_centers):
        f = rh.Geometry.Sphere(c, feat_radii[i])
        i = sc.doc.Objects.AddSphere(f)
        feat_ids += [i]
    
    # shearing features
    rs.CurrentView('Front')
    for i,(c, idx) in enumerate(zip(feat_centers, feat_ids)):
        ref_pt = cp.deepcopy(c)
        ref_pt[2] = feat_radii[i]
        rs.ShearObject(idx, c, ref_pt, shear_angle)
        
    # rotating features
    mean_orientation = 0.
    axis = [0,0,1]
    for c, idx in zip(feat_centers, feat_ids):
        angle = feature_orientation(mean_orientation)
        rs.RotateObject(idx, c, angle, axis)
    
    # creating a  z = 0 cutplane
    cutplane_pts = []
    cutplane_pts += [rh.Geometry.Point3d(xlim[0]-2*mean_radius, ylim[0]-2*mean_radius, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[1]+2*mean_radius, ylim[0]-2*mean_radius, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[1]+2*mean_radius, ylim[1]+2*mean_radius, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[0]-2*mean_radius, ylim[1]+2*mean_radius, 0)]
    aux_lname = 'auxilliary'
    rs.AddLayer(aux_lname) 
    rs.CurrentLayer(aux_lname)
    cplane_id = rs.AddSrfPt(cutplane_pts)
    
    # cutting features, delete cutplane
    rs.CurrentLayer(feat_lname)
    for i in feat_ids:
        split_ids = rs.SplitBrep(i, cplane_id, True)
        rs.DeleteObject(split_ids[1])
    rs.DeleteObject(cplane_id)
    
    # placing observation walls
    double_wall_dist = 7. # distance between two opposing walls
    wall_height = 2*mean_radius
    walls_east_no = 6
    dl_walls = (xlim[1]-xlim[0])/(walls_east_no-1)
    walls_east_lname = 'walls_east'
    rs.AddLayer(walls_east_lname) 
    rs.CurrentLayer(walls_east_lname)
    walls_east_ids = []
    for i in range(walls_east_no):
        dwds = [double_wall_dist/2., -double_wall_dist/2.]
        for dwd in dwds:
            pts = []
            pts += [rh.Geometry.Point3d(xlim[0]+i*dl_walls+dwd, ylim[0], -wall_height)]
            pts += [rh.Geometry.Point3d(xlim[0]+i*dl_walls+dwd, ylim[0], wall_height)]
            pts += [rh.Geometry.Point3d(xlim[0]+i*dl_walls+dwd, ylim[1], wall_height)]
            pts += [rh.Geometry.Point3d(xlim[0]+i*dl_walls+dwd, ylim[1], -wall_height)]
            walls_east_ids += [rs.AddSrfPt(pts)]
    walls_north_no = walls_east_no
    dl_walls = (ylim[1]-ylim[0])/(walls_north_no-1)
    walls_north_lname = 'walls_north'
    rs.AddLayer(walls_north_lname) 
    rs.CurrentLayer(walls_north_lname)
    walls_north_ids = []
    for i in range(walls_north_no):
        dwds = [double_wall_dist/2., -double_wall_dist/2.]
        for dwd in dwds:
            pts = []
            pts += [rh.Geometry.Point3d(xlim[0], ylim[0]+i*dl_walls+dwd, -wall_height)]
            pts += [rh.Geometry.Point3d(xlim[0], ylim[0]+i*dl_walls+dwd, wall_height)]
            pts += [rh.Geometry.Point3d(xlim[1], ylim[0]+i*dl_walls+dwd, wall_height)]
            pts += [rh.Geometry.Point3d(xlim[1], ylim[0]+i*dl_walls+dwd, -wall_height)]
            walls_north_ids += [rs.AddSrfPt(pts)]
    
    # intersections
    walls_north_int_lname = 'walls_north_int'
    rs.AddLayer(walls_north_int_lname) 
    rs.CurrentLayer(walls_north_int_lname)
    objects = rs.ObjectsByLayer(walls_north_lname)
    objects += rs.ObjectsByLayer(feat_lname)
    rs.SelectObjects(objects)
    rs.Command('_Intersect')
    rs.UnselectAllObjects()
    walls_east_int_lname = 'walls_east_int'
    rs.AddLayer(walls_east_int_lname) 
    rs.CurrentLayer(walls_east_int_lname)
    objects = rs.ObjectsByLayer(walls_east_lname)
    objects += rs.ObjectsByLayer(feat_lname)
    rs.SelectObjects(objects)
    rs.Command('_Intersect')
    rs.UnselectAllObjects()
    update_views()
    
    # deleting circular intersections (should be very few)
    delete_zero_width_elements(rs.ObjectsByLayer(walls_east_int_lname))
    delete_zero_width_elements(rs.ObjectsByLayer(walls_north_int_lname))
    
            
    # analyze
    objects = rs.ObjectsByLayer(walls_east_int_lname)
    walls_east_int_no = len(objects)
    walls_east_int_data = [analyze_intersection(iobj) for iobj in objects]    
    objects = rs.ObjectsByLayer(walls_north_int_lname)
    walls_north_int_no = len(objects)
    walls_north_int_data = [analyze_intersection(iobj) for iobj in objects]
    
    walls_east_int_skewness = [skewness(idata) for idata in walls_east_int_data]
    walls_north_int_skewness = [skewness(idata) for idata in walls_north_int_data]
    
    walls_east_int_aratio = [aspect_ratio(idata) for idata in walls_east_int_data]
    walls_north_int_aratio = [aspect_ratio(idata) for idata in walls_north_int_data]
    
    walls_east_int_width = [width(idata) for idata in walls_east_int_data]
    walls_north_int_width = [width(idata) for idata in walls_north_int_data]
    
    output_list('skewness_east.txt', walls_east_int_skewness)
    output_list('skewness_north.txt', walls_north_int_skewness)  
    output_list('aratio_east.txt', walls_east_int_aratio)
    output_list('aratio_north.txt', walls_north_int_aratio)
    output_list('width_east.txt', walls_east_int_width)
    output_list('width_north.txt', walls_north_int_width)
    
    print 'North wall intersections:', walls_north_int_no
    print 'East wall intersections:', walls_north_int_no


if __name__== '__main__':
    rd.seed(0)
    for m in range(1):
        main()

