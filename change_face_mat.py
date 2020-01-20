import bpy

line_mat = bpy.data.materials.get("line_mat")
if not line_mat:
    line_mat = bpy.data.materials.new("line_mat")
    line_mat.diffuse_color = (0.0, 1.0, 0.0)
    line_mat.type = 'WIRE'

cube_obj = bpy.data.objects['Cube']

if line_mat.name not in cube_obj.data.materials:
    cube_obj.data.materials.append( line_mat )

for i in range (0, len(cube_obj.data.polygons), 2):
    cube_obj.data.polygons[i].material_index = 1
