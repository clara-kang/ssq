import bpy
import random

def getRandomColor():
    return (random.random(), random.random(), random.random())

mats = [None] * len(patch_ids)

patch_num = max(patch_ids)+1

# init mats
for mat_id in range (patch_num+1) : # create one extra, since mat0 is default
    mats[mat_id] = bpy.data.materials.get("patch_mat_" + str(mat_id))
    if not mats[mat_id]:
        mats[mat_id] = bpy.data.materials.new("patch_mat_" + str(mat_id))
        mats[mat_id].diffuse_color = getRandomColor()

obj = bpy.data.objects['Sphere3']

for mat_id in range (patch_num+1) :
    if mats[mat_id].name not in obj.data.materials:
        obj.data.materials.append( mats[mat_id] )

# clear
for i in range (0, len(obj.data.polygons)):
    poly = obj.data.polygons[i]
    poly.material_index = 0

for i in range (0, len(obj.data.polygons)):
    poly = obj.data.polygons[i]
    for v_indx in poly.vertices:
        if v_indx in sl_verts:
            continue
        else:
            patch_id = patch_ids[v_indx] + 1
            poly.material_index = patch_id
            break
