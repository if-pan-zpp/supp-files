#!/usr/bin/python3

REPS = 5

import os, sys
import subprocess
import time
from shutil import copy

test = sys.argv[1]
binName = sys.argv[2]
args = sys.argv[3:]

assert os.path.exists(test)
assert os.path.isdir(test)
assert os.path.exists(binName)
root = os.getcwd()

label = input('Enter program label: ')

with open('results.txt', 'a') as out:
    copy(binName, os.path.join(root, test, os.path.basename(binName)))

    os.chdir(os.path.join(root, test))

    for num_threads in [1, 2, 4, 8]:
        timings = []
        for rep in range(1, REPS + 1):
            print('Running test', test + '...',
                    str(rep) + '/' + str(REPS),
                    str(num_threads), 'thread(s)')

            start_time = time.time()
            new_env = os.environ.copy()
            new_env['OMP_NUM_THREADS'] = str(num_threads)
            subprocess.run(['./' + binName] +
                            args,
                            capture_output=True,
                            env = new_env)
            wall_time = time.time() - start_time
            timings.append(wall_time)

        out.write('{}, {}, {}, {:.3f}\n'.format(label, test, num_threads, min(timings)))
        out.write('# ' + str(['{:.3f}'.format(time) for time in timings]) + '\n')
    out.write("-----\n")


