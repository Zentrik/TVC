import numpy as np
import bpy
#from os import listdir

scene = bpy.context.scene

path = bpy.path.abspath("//src/")
#path = PATH + "/" + sorted(listdir(PATH))[-1]

#iteration_folders = filter(lambda s: s.isdigit(),listdir(path))
#path += "/" + sorted(iteration_folders)[-1]
print(path)

FPS = scene.render.fps
scene.frame_current = 0

try:
    X = np.load(path + "x.npy") 
    U = np.load(path + "u.npy") 
    t = np.load(path + "t.npy") 
    t = t[0]
except OSError:
    print("Data not found.")
    #continue

K = X.shape[0]
scene.frame_end = int(t[-1] * FPS) #we will have time * fps frames

body_ob = bpy.data.objects.get("Starship")
eng1_ob = bpy.data.objects.get("Engine1")
eng2_ob = bpy.data.objects.get("Engine2")
eng3_ob = bpy.data.objects.get("Engine3")
fir1_ob = bpy.data.objects.get("Fire1")
fir2_ob = bpy.data.objects.get("Fire2")
fir3_ob = bpy.data.objects.get("Fire3")
light_ob = bpy.data.lights.get("Point")
pad_ob = bpy.data.objects.get("Ship")

body_ob.animation_data_clear()
eng1_ob.animation_data_clear()
eng2_ob.animation_data_clear()
eng3_ob.animation_data_clear()
fir1_ob.animation_data_clear()
fir2_ob.animation_data_clear()
fir3_ob.animation_data_clear()
light_ob.animation_data_clear()
pad_ob.animation_data_clear()

#scene.frame_current = 0
#pad_ob.location = np.append(X[-1, 0:2] / 10, 1.327) # moves landing pad to where the rocket landed as we are not targeting a location

T_max = np.max(np.linalg.norm(U[:, 1:3], axis=1))

for k in range(K):
    scene.frame_current = int(t[k] * FPS)
    x = X[k]
    print(k)
    u = U[k]
    
    body_ob.location = x[0:3]
    
    body_ob.rotation_quaternion = x[6:10]
    
    body_ob.keyframe_insert(data_path='location')
#    for fcurve in body_ob.animation_data.action.fcurves:
#            kf = fcurve.keyframe_points[-1]
#            kf.interpolation = 'QUAD'
            
    body_ob.keyframe_insert(data_path='rotation_quaternion')

    rx = np.arctan(-u[1] / u[2])
    if np.isnan(rx):
        rx = 0.0
        
    ry = np.arctan(u[0] / u[2])
    if np.isnan(ry):
        ry = 0.0
        
    min_l = 0.0
    l = min_l + (1.-min_l) * (np.linalg.norm(u[0:3]) / T_max)
    l *= 1
    
    for eng in [eng1_ob, eng2_ob, eng3_ob]:
        eng.rotation_euler = (rx, ry, 0)
        eng.keyframe_insert(data_path='rotation_euler')
        
    for fir in [fir1_ob, fir2_ob, fir3_ob]:
        fir.scale[2] = l
        fir.keyframe_insert(data_path='scale')
        for fcurve in fir.animation_data.action.fcurves:
            kf = fcurve.keyframe_points[-1]
            kf.interpolation = 'LINEAR'
        
    light_ob.energy = fir1_ob.scale[2] * 50
    light_ob.keyframe_insert(data_path='energy')
        
scene.frame_current += FPS // 4

for fir in [fir1_ob, fir2_ob, fir3_ob]:
    fir.scale[2] = 0
    fir.keyframe_insert(data_path='scale')

light_ob.energy = 0
light_ob.keyframe_insert(data_path='energy')

scene.frame_current = 0
