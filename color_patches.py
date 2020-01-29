patch_ids = [17, 17, 17, 19, 19, 19, 17, 17, 5, 5, 17, 17, 5, 5, 17, 17, 5, 5, 17, 17, 5, 5, 19, 19, 19, 17, 17, 17, 19, 19, 17, 19, 17, 19, 17, 19, 17, 5, 17, 5, 17, 5, 17, 5, 19, 17, 19, 17, 17, 5, 5, 17, 10, 19, 17, 19, 17, 19, 17, 19, 17, 5, 19, 17, 19, 17, 5, 17, 10, 17, 10, 17, 10, 19, 10, 19, 17, 10, 10, 17, 10, 19, 10, 19, 8, 19, 8, 19, 17, 10, 17, 17, 8, 19, 19, 8, 10, 17, 10, 17, 10, 17, 10, 19, 10, 19, 8, 17, 10, 10, 19, 10, 19, 8, 19, 8, 19, 8, 10, 17, 10, 17, 8, 19, 10, 17, 10, 17, 10, 17, 10, 19, 10, 19, 8, 19, 8, 19, 9, 10, 10, 19, 8, 19, 8, 19, 8, 10, 9, 10, 9, 10, 11, 9, 9, 10, 9, 10, 9, 10, 19, 10, 19, 8, 19, 8, 19, 10, 10, 19, 19, 8, 19, 9, 19, 9, 11, 9, 11, 9, 10, 9, 9, 9, 11, 9, 11, 9, 11, 19, 11, 19, 9, 19, 9, 19, 11, 19, 11, 9, 19, 9, 19, 9, 11, 9, 11, 9, 11, 9, 11, 19, 11, 11, 9, 11, 9, 11, 19, 11, 19, 9, 19, 9, 19, 9, 9, 11, 9, 19, 19, 9, 19, 9, 11, 9, 11, 9, 11, 9, 11, 19, 9, 9, 11, 9, 11, 19, 11, 19, 9, 19, 9, 19, 9, 11, 11, 9, 9, 19, 19, 9, 14, 9, 14, 9, 14, 9, 14, 19, 14, 19, 9, 9, 14, 19, 14, 12, 14, 12, 14, 14, 9, 14, 9, 14, 14, 14, 9, 14, 9, 14, 9, 14, 9, 14, 12, 14, 12, 14, 12, 14, 14, 12, 9, 14, 12, 14, 12, 14, 14, 9, 14, 9, 14, 9, 14, 14, 14, 14, 9, 14, 9, 14, 9, 14, 12, 14, 12, 14, 12, 14, 12, 9, 14, 14, 12, 14, 12, 14, 12, 14, 14, 9, 14, 9, 14, 12, 14, 14, 14, 17, 14, 17, 14, 12, 14, 12, 14, 12, 14, 14, 18, 12, 18, 12, 18, 12, 18, 12, 18, 18, 17, 18, 17, 18, 12, 18, 17, 12, 18, 17, 18, 17, 18, 12, 18, 12, 18, 12, 18, 12, 18, 12, 18, 12, 18, 12, 17, 18, 17, 18, 17, 18, 12, 18, 12, 17, 17, 18, 18, 17, 18, 12, 18, 12, 18, 12, 18, 12, 18, 18, 12, 18, 18, 18, 18, 18, 17, 18, 17, 18, 17, 18, 12, 18, 12, 18, 18, 17, 17, 18, 12, 18, 12, 18, 12, 18, 18, 17, 18, 17, 18, 17, 19, 18, 17, 18, 17, 18, 17, 18, 17, 18, 19, 18, 19, 17, 17, 18, 19, 18, 19, 17, 19, 17, 19, 17, 18, 17, 18, 18, ]


sl_verts = [54, 52, 70, 68, 89, 94, 116, 123, 54, 76, 63, 88, 97, 119, 125, 54, 74, 82, 104, 111, 132, 54, 45, 27, 26, 471, 465, 452, 190, 188, 186, 178, 152, 145, 123, 190, 192, 173, 153, 148, 125, 190, 168, 161, 140, 132, 190, 198, 219, 227, 248, 270, 260, 259, 271, 275, 23, 260, 253, 232, 224, 202, 195, 174, 166, 145, 123, 260, 262, 264, 266, 268, 270, 260, 281, 287, 310, 318, 342, 348, 371, 370, 351, 353, 25, 129, 127, 125, 370, 349, 341, 317, 309, 301, 280, 278, 270, 370, 368, 366, 364, 362, 375, 373, 371, 370, 377, 396, 420, 426, 447, 452, 476, 454, 466, 464, 23, 476, 3, 5, 35, 59, 64, 87, 94, 116, 123, 476, 455, 448, 425, 421, 399, 380, 371, 476, 478, 480, 481, 469, 471, 465, 452, ]

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
            patch_id = patch_ids[v_indx]
            poly.material_index = patch_id
            break
