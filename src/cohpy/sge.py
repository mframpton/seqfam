import pandas as pd
import numpy as np
from datetime import datetime
import os
import re
import glob
import sys
from natsort import natsorted
import subprocess
from misc import Logger


class SGE(object):
	
	def __init__(self, scripts_dir):
		self.logger = Logger()
		self.scripts_dir = scripts_dir
		
		
	def make_map_reduce_jobs(self, prep, map_task_l, reduce_task, map_task_exec_l=[], mem="14G"):

		'''Delete existing .sh files in the scripts dir.'''
		for f in glob.glob(os.path.join(self.scripts_dir,"*.sh")):
			os.remove(f)
	
		'''Make map.sh'''
		map_file = os.path.join(self.scripts_dir,".".join([prep,"map","sh"]))
		map_stream = open(map_file,'w')
		map_stream.write(self.get_sge_job_txt(prep, map_cmd_n=len(map_task_exec_l), mem=mem))
		map_stream.close()
	
		'''Make reduce.sh'''
		reduce_file = os.path.join(self.scripts_dir,".".join([prep,"reduce","sh"]))
		reduce_stream = open(reduce_file,'w')
		reduce_stream.write(self.get_sge_job_txt(prep, reduce_task=reduce_task, mem=mem))
		reduce_stream.close()
		
		'''Make bash task scripts for each map task and a text file containing references to them.'''
		self.logger.log("# of map tasks to execute: {0}".format(len(map_task_exec_l)))
		all_map_task_stream = open(os.path.join(self.scripts_dir,".".join([prep,"all_map_task","txt"])), 'w') 	
		for i in xrange(len(map_task_l)):
			task_script = os.path.join(self.scripts_dir,".".join([prep,str(i+1),"sh"]))
			task_stream = open(task_script, 'w')
			task_stream.write("#!/bin/bash\n\n")
			task_stream.write(map_task_l[i])
			task_stream.close()
			if not map_task_exec_l == False:
				if map_task_l[i] in map_task_exec_l:
					all_map_task_stream.write("bash " + task_script + "\n")
			else:
				all_map_task_stream.write("bash " + task_script + "\n") 
		all_map_task_stream.close()	
	
		'''Make the submit script.'''
		submit_stream = open(os.path.join(self.scripts_dir,"submit_map_reduce.sh"),'w')
		submit_stream_l = ["#!/bin/bash"]
		submit_stream_l.append("\n")
		if len(map_task_exec_l) == 0:
			submit_stream_l.append(" ".join(["qsub",reduce_file]))
		else:
			submit_stream_l.append(" ".join(["qsub","-N",prep+".map","-cwd",map_file]))
			submit_stream_l.append(" ".join(["qsub","-hold_jid",prep+".map","-cwd",reduce_file]))
		submit_stream.write("\n".join(submit_stream_l))
		submit_stream.close()
		if os.name != "nt":
			subprocess.call(['chmod', 'u+x', os.path.join(self.scripts_dir,"submit_map_reduce.sh")])
		
	
	def get_sge_job_txt(self, prep, map_cmd_n=None, reduce_task=None, mem="14G"):

		oe_txt_l = ["mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err","exec >${scriptname}.qsub.out/${scriptname}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${JOB_ID}.err"]
		sge_job_txt_l = ["#!/bin/bash","#$ -S /bin/bash","#$ -o /dev/null","#$ -e /dev/null","#$ -cwd","#$ -V","#$ -l tmem="+mem+",h_vmem="+mem,"#$ -l h_rt=24:0:0","set -u","set -x","\n"]
		if map_cmd_n != None:
			sge_job_txt_l.insert(-3,"#$ -t 1:"+str(map_cmd_n))
			sge_job_txt_l.extend(["scriptname=" + prep + ".$SGE_TASK_ID"] + oe_txt_l + ['CMDS=`sed -n -e "$SGE_TASK_ID p" ' + prep + '.all_map_task.txt`',"$CMDS"])
		elif reduce_task != None:
			sge_job_txt_l.extend(["scriptname=" + prep + ".reduce"] + oe_txt_l + [reduce_task])
	
		sge_job_txt_str = "\n".join(sge_job_txt_l)
	
		return sge_job_txt_str	
