import bpy

def color_vertices(obj):
    """Paints a single vertex where vert is the index of the vertex
    and color is a tuple with the RGB values."""

    mesh = obj.data
    scn = bpy.context.scene

    #check if our mesh already has Vertex Colors, and if not add some... (first we need to make sure it's the active object)
    scn.objects.active = obj
    obj.select = True
    if mesh.vertex_colors:
        vcol_layer = mesh.vertex_colors.active
    else:
        vcol_layer = mesh.vertex_colors.new()

    for poly in mesh.polygons:
        for loop_index in poly.loop_indices:
            loop_vert_index = mesh.loops[loop_index].vertex_index
            vcol_layer.data[loop_index].color = cols[loop_vert_index]

#example usage
# for i in range (len(cols)):
color_vertices(bpy.context.scene.objects["sphere1"])
