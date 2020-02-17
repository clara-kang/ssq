max_locs = [None] * 6
sdl_locs = [None] * 10
min_locs = [None] * 6

import bpy
import bmesh

mats = [None, None, None]

cols = [(1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, 0.0, 1.0)]
for i in range (0, 3) :
    mats[i] = bpy.data.materials.get("mat_" + str(i))
    if not mats[i]:
        mats[i] = bpy.data.materials.new("mat_" + str(i))
        mats[i].diffuse_color = cols[i]

locs = [max_locs, sdl_locs, min_locs]

bpyscene = bpy.context.scene

# create empty object as parent
points_holder = bpy.data.objects.new('points_holder', None)
bpyscene.objects.link(points_holder)

for i in range (0, 3) :
    for loc in locs[i]:
        mesh = bpy.data.meshes.new('Basic_Sphere')
        basic_sphere = bpy.data.objects.new("Basic_Sphere", mesh)
        bpyscene.objects.link(basic_sphere)
        bpyscene.objects.active = basic_sphere
        basic_sphere.select = True
        basic_sphere.parent = points_holder
        basic_sphere.location = loc
        basic_sphere.active_material = mats[i]
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
