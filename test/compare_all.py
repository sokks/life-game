import sys
from os import listdir
from os.path import isfile, join


if len(sys.argv) < 5:
    print('usage: visualize.py <N> <M> <folder1> <folder2>')
    exit()
M, N = int(sys.argv[1]), int(sys.argv[2])
folder1, folder2 = sys.argv[3], sys.argv[4]


def read_scenes_in_dir(fld):
    round_files = [f for f in listdir(fld) if isfile(join(fld, f))]
    scenes = []
    for f in sorted(round_files):
        Z = []
        with open(join(fld, f)) as fin:
            row = []
            for line in fin:
                ll = line.split()
                row = [int(e) for e in ll]
                Z += row
        scenes.append(Z)
    return scenes


scenes1 = read_scenes_in_dir(folder1)
scenes2 = read_scenes_in_dir(folder2)
assert len(scenes1) == len(scenes2), 'different number of files'
for i in range(len(scenes1)):
    assert scenes1[i] == scenes2[i], 'different scenes: {}'.format(i)

print('TESTS:     OK')
