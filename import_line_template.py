import bpy
import bmesh

scene = bpy.context.scene

line_holder = bpy.data.objects.new('lines_holder', None)
scene.objects.link(line_holder)

line_mat = bpy.data.materials.new("mat_line")
line_mat.diffuse_color = (0.0, 1.0, 0.0)
line_mat.type = 'WIRE'
line_mat.use_shadeless

for line in range(0, len(verts)):
    mesh = bpy.data.meshes.new("mesh" + str(line))  # add a new mesh
    obj = bpy.data.objects.new("line" + str(line), mesh)  # add a new object using the mesh
    scene.objects.link(obj)  # put the object into the scene (link)
    scene.objects.active = obj  # set as the active object in the scene
    obj.select = True  # select object
    obj.parent = line_holder
    obj.active_material = line_mat
    mesh = bpy.context.object.data
    bm = bmesh.new()

    bm_verts = [None] * len(verts[line])
    cnt = 0
    for v in verts[line]:
        bm_verts[cnt] = bm.verts.new(v)  # add a new vert
        cnt += 1

    for v_indx in range(0, len(bm.verts) - 1):
        bm.edges.new((bm_verts[v_indx], bm_verts[v_indx+1]))

    # make the bmesh the object's mesh
    bm.to_mesh(mesh)
    bm.free()  # always do this when finished
