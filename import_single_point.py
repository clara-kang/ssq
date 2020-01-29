loc = (1.538041, 0.000000, -0.382683)

import bpy
import bmesh

mat = bpy.data.materials.get("mat_single_pt")
if not mat:
    mat = bpy.data.materials.new("mat_single_pt")
    mat.diffuse_color = (1.0, 0.647, 0)

bpyscene = bpy.context.scene

mesh = bpy.data.meshes.new('Basic_Sphere')
basic_sphere = bpy.data.objects.new("sp", mesh)
bpyscene.objects.link(basic_sphere)
bpyscene.objects.active = basic_sphere
basic_sphere.select = True
basic_sphere.location = loc
basic_sphere.active_material = mat
bm = bmesh.new()
bmesh.ops.create_icosphere(bm, subdivisions=1, diameter=0.02)
bm.to_mesh(mesh)
bm.free()
