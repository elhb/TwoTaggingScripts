#! /bin/env python

small = False

import os
import sys
#get files

files = []
for filename in os.listdir(sys.argv[1]+'/fwd'):
    if filename.endswith(".fa"): files.append(sys.argv[1]+'/fwd/'+filename)
for filename in os.listdir(sys.argv[1]+'/rev'):
    if filename.endswith(".fa"): files.append(sys.argv[1]+'/rev/'+filename)

cores = 16

os.mkdir(sys.argv[1]+'/sbatch')
tmpCounter = 0
breakFlag=False
while True:
    print 'creating sbatch for '+str(tmpCounter)+' out of '+str(len(files)/cores)
    try:
        current_block = files[:cores]
        files=files[cores:]
    except IndexError:
        current_block = files
        break
    f = open(sys.argv[1]+'/sbatch/sbatch.'+str(tmpCounter)+'.sh','w')
    f.write('#! /bin/bash -l'+'\n')
    f.write('#SBATCH -A b2011168\n')
    f.write('#SBATCH -n 8 -p node'+'\n')
    if not small:
            f.write('#SBATCH -t 24:00:00'+'\n')
    else:
            f.write('#SBATCH -t 1:00:00'+'\n')
    f.write('#SBATCH -J muscle_'+str(tmpCounter)+'\n')
    f.write('#SBATCH -e '+os.path.abspath(sys.argv[1])+'/sbatch/sbatch.'+str(tmpCounter)+'.stderr.txt'+'\n')
    f.write('#SBATCH -o '+os.path.abspath(sys.argv[1])+'/sbatch/sbatch.'+str(tmpCounter)+'.stdout.txt'+'\n')
    f.write('#SBATCH --mail-type=All'+'\n')
    f.write('#SBATCH --mail-user=erik.borgstrom@scilifelab.se'+'\n')
    f.write('echo "$(date) Running on: $(hostname)"'+'\n')
    f.write('cd '+os.getcwd()+'\n')
    f.write('module load bioinfo-tools muscle/3.8.31'+'\n')
    for filename in current_block:f.write('muscle -in '+filename+' -out '+filename+'.aligned -quiet &\n')
    f.write('wait\n')
    f.close()
    if len(files) == 0: break
    tmpCounter+=1