from matplotlib import pyplot as plt
from matplotlib import animation
from os import listdir
from os.path import isfile, join
import sys
import numpy as np

if len(sys.argv) < 3:
    print('usage: visualize.py <N> <M>')
    exit()
M, N = int(sys.argv[1]), int(sys.argv[2])

rounds_dir = 'rounds'
round_files = [f for f in listdir(rounds_dir) if isfile(join(rounds_dir, f))]
scenes = []
for f in sorted(round_files):
    Z = []
    with open(join(rounds_dir, f)) as fin:
        row = []
        for line in fin:
            ll = line.split()
            row = [int(e) for e in ll]
            Z += row
    scenes.append(Z)

print(len(scenes))

def scene_generator(scenes):
    for scene in scenes:
        yield scene

sg = scene_generator(scenes)
scene = next(sg)
Z = np.array(scene)
Z = Z.reshape((N, M))
print(Z.shape)

plt.style.use('ggplot')

fig = plt.figure()
im = plt.imshow(Z, interpolation='nearest', animated=True, cmap=plt.get_cmap('terrain'))
plt.axis('off')

def updatefig(*args):
    global Z
    scene = next(sg)
    Z = np.array(scene)
    Z = Z.reshape((N, M))
    im.set_array(Z)
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=100, blit=True, frames=len(scenes), repeat=False)
plt.show()
