vert_ids = []

import bpy
import bmesh

def find_vert(vert_id):
    file_variable = open('C:/Users/Clara/Documents/ssq1/models/high_sphere.obj')
    all_lines_variable = file_variable.readlines()
    line = all_lines_variable[vert_id + 3]
    segs = line.split()
    loc = (float(segs[1]), float(segs[2]), float(segs[3]))
    print("loc: ", loc)
    return loc

mat = bpy.data.materials.get("mat_single_pt")
if not mat:
    mat = bpy.data.materials.new("mat_single_pt")
    mat.diffuse_color = (1.0, 0.647, 0)

bpyscene = bpy.context.scene

for v in vert_ids:
    mesh = bpy.data.meshes.new('Basic_Sphere')
    basic_sphere = bpy.data.objects.new("v_" + str(v), mesh)
    bpyscene.objects.link(basic_sphere)
    bpyscene.objects.active = basic_sphere
    basic_sphere.select = True
    basic_sphere.location = find_vert(v)
    basic_sphere.active_material = mat
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=1, diameter=0.02)
    bm.to_mesh(mesh)
    bm.free()
