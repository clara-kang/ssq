nodes_locs = []

import bpy
import bmesh

mat = bpy.data.materials.get("node_mat")
if not mat:
    mat = bpy.data.materials.new("node_mat")
    mat.diffuse_color = (0.65882353, 0.19607843, 0.46666667)

bpyscene = bpy.context.scene

# create empty object as parent
points_holder = bpy.data.objects.new('points_holder', None)
bpyscene.objects.link(points_holder)

for loc in nodes_locs:
    mesh = bpy.data.meshes.new('Basic_Sphere')
    basic_sphere = bpy.data.objects.new("Basic_Sphere", mesh)
    bpyscene.objects.link(basic_sphere)
    bpyscene.objects.active = basic_sphere
    basic_sphere.select = True
    basic_sphere.parent = points_holder
    basic_sphere.location = loc
    basic_sphere.active_material = mat
    bm = bmesh.new()
    bmesh.ops.create_icosphere(bm, subdivisions=1, diameter=0.02)
    bm.to_mesh(mesh)
    bm.free()

# for fldr in range (0, 3):
#     for i in range (0, file_nums[fldr]):
#         file_loc = 'C:\\Users\\Clara\\Documents\\ssq1\\models\\' + dir_names[fldr] + '\\'
#         imported_object = bpy.ops.import_scene.obj(filepath=file_loc+"icosphere" + str(i) + ".obj")
#         obj_object = bpy.context.selected_objects[0] ####<--Fix
#         obj_object.parent = points_holder
#         obj_object.active_material = mats[fldr]
