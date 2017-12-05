import os
from sge.sge import SGE
import sys


def get_map_task_l(chr_l, data_dir):

    map_task_l, map_task_exec_l = [],[]
    for i in xrange(len(chr_l)):
        map_task = " ".join(["python","1_map_task.py","--chr",chr_l[i]])
        map_task_l.append(map_task)
        if os.path.isfile(os.path.join(data_dir,".".join(["map_task",chr_l[i],"out"]))) == False:
            map_task_exec_l.append(map_task)

    return [map_task_l,map_task_exec_l]


'''Set the script and data directory.'''
script_dir = os.path.abspath(os.path.join("..","..","data","sge"))
data_dir = os.path.abspath(os.path.join("..","..","data","sge"))

chr_l = [str(chr) for chr in xrange(1,23)] + ["X","Y"]
print "Making map tasks..."
[map_task_l, map_task_exec_l] = get_map_task_l(chr_l,data_dir)
print "Making reduce tasks..."
reduce_task_l = [" ".join(["python","2_summarise_results.py"])]
reduce_task_l.append(" ".join(["python","3_summarise_results.py"]))
reduce_task = "\n".join(reduce_task_l)
#print map_task_l
#print map_task_exec_l
#print reduce_task

print "Writing job scripts..."
sge = SGE(script_dir)
sge.make_map_reduce_jobs("test", map_task_l, reduce_task, map_task_exec_l)
