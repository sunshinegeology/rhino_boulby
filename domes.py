
import Rhino as rh
import rhinoscriptsyntax as rs
import scriptcontext as sc
import System.Guid as sg
import copy as cp
import math as ma
import random as rd
import sys
import json
import os


def Ellipsoid(c, radius, aratio):
    """c is centerpoint, radius is short side along in x-dir, aratio is long side to short side or radius in y-dir."""
    pt0, pt1, pt2 = cp.deepcopy(c), cp.deepcopy(c), cp.deepcopy(c)
    pt0[1] += radius*aratio
    pt1[0] += radius
    pt2[2] += radius

    cmd = "_Ellipsoid " + str(c) + " " + str(pt0) + " " + str(pt1) + " " + str(pt2)    
    blnResult = rs.Command(cmd, echo=False)

    if blnResult == True:
        return rs.LastCreatedObjects()[0]
    raise IOError('Cmdline construction of ellipsoid didnt work')


def bend_selected(beg_spine, end_spine, bend_pt):
    cmd = "_Bend " + str(beg_spine) + " " + str(end_spine) + " Symmetric=Yes " + str(bend_pt)
    rs.Command(cmd, echo=False)

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


def feature_radius(rmin, rmax, rsigma):
    mean_radius = (rmax+rmin)/2
    radius = 0.
    while radius < rmin or radius > rmax:
        radius = rd.gauss(mean_radius, rsigma)
    return radius


def feature_orientation(mean_orientation, sigma_rotation):
    return mean_orientation + rd.gauss(mean_orientation, sigma_rotation)


def non_intersecting_center_points(pt_no, xlim, ylim, feature_center_min_dist):
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
    

def finalize_views():
    views = rs.ViewNames()
    for view in views:
        rs.ViewDisplayMode(view, 'Ghosted')
    sc.doc.Views.Redraw()


def main(settings):
    # resetting
    rs.DocumentModified(False)
    rs.Command('_-New _None')
    update_views()
    
    # settings    
    half_length = settings['half length']
    xlim, ylim = (-half_length,half_length), (-half_length,half_length)
    feat_no = int(settings['density']*(half_length*2.)**2)
    shear_angle = settings['shear angle']
    rmin, rmax, rsigma = settings['rmin'], settings['rmax'], settings['rsigma']
    sigma_rotation = settings['orientation sigma']
    mean_orientation = settings['orientation mean']
    walls_no = settings['walls']
    double_wall_dist = settings['opposing walls distance']
    height_scaling_factor = settings['height scaling factor']
    is_domes = settings['dome']
    feature_center_min_dist = rmax*2
    if not is_domes:
        ridges_aratio = settings['ridges aspect ratio']
        feature_center_min_dist = rmax*ridges_aratio*2
        ridges_bangle = settings['ridges bending angle']
    # adding layers
    feat_lname = 'domes'
    l = rs.AddLayer(feat_lname)
    rs.CurrentLayer(feat_lname)
    
    # adding features
    feat_ids = []
    feat_radii = [feature_radius(rmin, rmax, rsigma) for i in range(feat_no)]
    rs.CurrentLayer(feat_lname)
    feat_centers = non_intersecting_center_points(feat_no, xlim, ylim, feature_center_min_dist)
    for i, c in enumerate(feat_centers):
        if is_domes:
            f = rh.Geometry.Sphere(c, feat_radii[i])
            i = sc.doc.Objects.AddSphere(f)  
        else:
            i = Ellipsoid(c, feat_radii[i], ridges_aratio)
        feat_ids += [i]
    
    # shearing features
    rs.CurrentView('Front')
    for i,(c, idx) in enumerate(zip(feat_centers, feat_ids)):
        ref_pt = cp.deepcopy(c)
        ref_pt[2] = feat_radii[i]
        rs.ShearObject(idx, c, ref_pt, shear_angle)
    
    # bending ellipses
    pi = 3.14159
    if not is_domes:
        rs.CurrentView('Top')
        rs.UnselectAllObjects()
        for i in range(len(feat_centers)):
            # orientation of ellipse hardcoded here
            cpoint = feat_centers[i]
            beg_spine, end_spine = cp.deepcopy(cpoint), cp.deepcopy(cpoint)
            beg_spine[1] += feat_radii[i]*ridges_aratio
            end_spine[1] -= feat_radii[i]*ridges_aratio
            bend_pt = cp.deepcopy(beg_spine)
            r = feat_radii[i]*ridges_aratio*2
            t = 3.*pi/2.+ridges_bangle*pi/180.
            bend_pt[0] += r * ma.cos(t)
            bend_pt[1] += r * ma.sin(t)
            rs.SelectObjects(feat_ids[i])
            bend_selected(beg_spine, end_spine, bend_pt)
            rs.UnselectAllObjects()
            
        
    # rotating features
    axis = [0,0,1]
    for c, idx in zip(feat_centers, feat_ids):
        angle = feature_orientation(mean_orientation, sigma_rotation)
        rs.RotateObject(idx, c, angle, axis)
    
    # creating a  z = 0 cutplane
    cutplane_pts = []
    cutplane_pts += [rh.Geometry.Point3d(xlim[0]-2*rmax, ylim[0]-2*rmax, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[1]+2*rmax, ylim[0]-2*rmax, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[1]+2*rmax, ylim[1]+2*rmax, 0)]
    cutplane_pts += [rh.Geometry.Point3d(xlim[0]-2*rmax, ylim[1]+2*rmax, 0)]
    aux_lname = 'auxilliary'
    rs.AddLayer(aux_lname) 
    rs.CurrentLayer(aux_lname)
    cplane_id = rs.AddSrfPt(cutplane_pts)
    
    # scaling features in z direction only
    rs.CurrentView('Front')
    for c, idx in zip(feat_centers, feat_ids):
        rs.ScaleObject(idx, c, [1., height_scaling_factor, 1.])
        
    # cutting features, delete cutplane
    rs.CurrentLayer(feat_lname)
    for i in feat_ids:
        split_ids = rs.SplitBrep(i, cplane_id, True)
        if split_ids:
            rs.DeleteObject(split_ids[1])
    rs.DeleteObject(cplane_id)
    
    # placing observation walls
    double_wall_dist = 7. # distance between two opposing walls
    wall_height = 2*rmax
    walls_east_no = walls_no
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
    
    walls_east_int_height = [height(idata) for idata in walls_east_int_data]
    walls_north_int_height = [height(idata) for idata in walls_north_int_data]
    
    # north alog x-axis, east along y-axis
    output_list('skewness_east.txt', walls_east_int_skewness)
    output_list('skewness_north.txt', walls_north_int_skewness)  
    output_list('aratio_east.txt', walls_east_int_aratio)
    output_list('aratio_north.txt', walls_north_int_aratio)
    output_list('width_east.txt', walls_east_int_width)
    output_list('width_north.txt', walls_north_int_width)
    output_list('height_east.txt', walls_east_int_height)
    output_list('height_north.txt', walls_north_int_height)
    
    finalize_views()
    
    print 'Feature number', feat_no
    print 'North wall intersections:', walls_north_int_no
    print 'East wall intersections:', walls_north_int_no


if __name__== '__main__': 
    with open('rhino_settings.json', 'r') as f:
        settings = json.load(f)
    """
    settings = dict()
    settings['half length'] = 150. # [m] half length of model so that xmin/ymin = -half length and xmax/ymax = half length
    settings['density'] = 0.0001#0.0015 # features per square meter
    settings['shear angle'] = 10. # [deg]
    settings['rmin'] = 1. # [m] minimum dome radius
    settings['rmax'] = 2. # [m] maximum dome radius
    settings['rsigma'] = 1.5 # [m] standard deviation of radii distribution
    settings['orientation sigma'] = 10. # [deg] standard deviation of orientation distribuition
    settings['orientation mean'] = 0. # [deg] mean orientation of features w.r.t. x-axis
    settings['walls'] = 10 # number of double walls east and north
    settings['opposing walls distance'] = 7. # [m] distance between two opposing walls
    settings['height scaling factor'] = 0.29 # 0-1, where 1 means no scaling
    settings['dome'] = False #True # True for domes, False for ridges
    settings['ridges aspect ratio'] = 20.
    settings['ridges bending angle'] = 30.
    """
    
    rd.seed(0)
    main(settings)

