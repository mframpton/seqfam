import os
from seqfam.sge import SGE
import sys
#import pathlib


def get_map_task_l(chr_l):

    map_task_l, map_task_exec_l = [],[]
    for i in range(len(chr_l)):
        map_task = " ".join(["python","1_map_task.py","--chr",chr_l[i]])
        map_task_l.append(map_task)
        #The user can replace the if condition with one confirming output for the map task is absent, and hence it should be in map_task_exec_l. 
        if i % 2 == 0:
            map_task_exec_l.append(map_task)

    return [map_task_l,map_task_exec_l]


'''Set the script and data directory.'''
script_dir = os.path.abspath(os.path.join("..","..","data","sge"))
#pathlib.Path(script_dir).mkdir(parents=True, exist_ok=True) 
if not os.path.exists(script_dir):
    os.makedirs(script_dir)
data_dir = os.path.abspath(os.path.join("..","..","data","sge"))

chr_l = [str(chrom) for chrom in range(1,23)] + ["X","Y"]
print("Making map tasks...")
[map_task_l, map_task_exec_l] = get_map_task_l(chr_l)
print("Making reduce tasks...")
reduce_task_l = [" ".join(["python","2_merge_results.py"])]
reduce_task_l.append(" ".join(["python","3_summarise_results.py"]))
reduce_task = "\n".join(reduce_task_l)
#print(map_task_l)
#print(map_task_exec_l)
#print(reduce_task)

print("Writing job scripts...")
sge = SGE(script_dir)
sge.make_map_reduce_jobs("test", map_task_l, reduce_task, map_task_exec_l)
