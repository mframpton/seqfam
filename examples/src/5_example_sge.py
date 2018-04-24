import os
from seqfam.sge import SGE


def get_map_task_l(chr_l):

    map_task_l, map_task_exec_l = [],[]
    for i in range(len(chr_l)):
        map_task = " ".join(["python","1_map_task.py","--chr",chr_l[i]])
        map_task_l.append(map_task)
        #The user can replace the if condition with one confirming output for the map task is absent, and hence it should be in map_task_exec_l. 
        if i % 2 == 0:
            map_task_exec_l.append(map_task)

    return [map_task_l,map_task_exec_l]


print("Making map and reduce tasks...")
chr_l = [str(chrom) for chrom in range(1,23)] + ["X","Y"]
[map_task_l, map_task_exec_l] = get_map_task_l(chr_l)
reduce_tasks = "\n".join(["python 2_merge_results.py","python 3_summarise_results.py"])

print("Writing job scripts...")
script_dir = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data","sge")) 
if not os.path.exists(script_dir):
    os.makedirs(script_dir)
sge = SGE(script_dir)
sge.make_map_reduce_jobs("test", map_task_l, reduce_tasks, map_task_exec_l)
