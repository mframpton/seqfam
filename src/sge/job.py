import pandas as pd
import numpy as np
from datetime import datetime
import os
import re
import glob
import sys
from natsort import natsorted
import subprocess


def make_all_job_scripts(scripts_dir, job_name_l, job_script_l, cmd_l, merge_cmd):

	#print scripts_dir
	#print job_name_l
	#print job_script_l
	#print cmd_l
	#print make_merge_script
	for i in xrange(len(job_name_l)):
		make_job_script(job_script_l[i], job_name_l[i], cmd_l[i])
	master_script_stream = open(os.path.join(scripts_dir,"run_all_jobs.sh"), 'w')
   	master_script_stream.write("#!/bin/bash\n\n")
   	for i in xrange(len(job_name_l)):
		master_script_stream.write("qsub -N " + job_name_l[i] + " -cwd " + job_script_l[i] + "\n")
	if merge_cmd != None:
		make_job_script(os.path.join(scripts_dir,"merge_results.sh"), "merge", merge_cmd)
		master_script_stream.write("qsub -hold_jid " + ",".join(job_name_l) + " -cwd merge_results.sh")
	master_script_stream.close()


def make_job_script(job_script, job_name, command, mem="14G"):

    job_stream = open(job_script, 'w')
    job_stream.write("#!/bin/bash\n")
    job_stream.write("#$ -S /bin/bash\n")
    job_stream.write("#$ -o /dev/null\n")
    job_stream.write("#$ -e /dev/null\n")
    job_stream.write("#$ -cwd\n")
    job_stream.write("#$ -V\n")
    job_stream.write("#$ -l tmem="+mem+",h_vmem="+mem+"\n")
#job_stream.write("#$ -l tmem=14G,h_vmem=14G\n")
    job_stream.write("#$ -l h_rt=24:0:0\n")
    job_stream.write("set -u\n")
    job_stream.write("set -x\n\n")
    job_stream.write("scriptname=" + job_name + "\n")
    job_stream.write("mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err\n")
    job_stream.write("exec >${scriptname}.qsub.out/${scriptname}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${JOB_ID}.err\n\n")
    job_stream.write(command + "\n")
    job_stream.close()


def merge_results(data_dir, in_file_pat, out_file_prefix, index_col_l, skiprows=1, sep=",", concat_axis=0, sort_col_l=None, sort_col_numeric_l=None, descending=False, out_compression=False, write=True):

	df_l = []
	csv_l = glob.glob(os.path.join(data_dir, in_file_pat))
	csv_l = filter(lambda x: "all" not in os.path.basename(x), csv_l)
	#print csv_l
	csv_l = natsorted(csv_l)
	#print csv_l
	for csv in csv_l:
		print csv
		this_df = pd.read_csv(csv, skiprows=skiprows, index_col=index_col_l, dtype="str", na_filter=False, sep=sep)
        #print this_df.head()
		if this_df.empty:
			continue
		df_l += [this_df]
	df = pd.concat(df_l, axis=concat_axis)
	df.reset_index(inplace=True)

	if sort_col_l != None:
		for i in xrange(len(sort_col_l)):
			if sort_col_numeric_l[i] == True:
				df[sort_col_l[i]] = df[sort_col_l[i]].apply(pd.to_numeric, errors="coerce")
		ascending = not descending
		df.sort_values(by=sort_col_l, inplace=True, ascending=ascending)
	
	df.set_index(index_col_l, inplace=True)

	if write == True:
		out_path = os.path.join(data_dir, ".".join(["all",out_file_prefix]))
		out_path = out_path + ".csv" if sep == "," else out_path + ".tsv"
		print "Writing to " + out_path + "..."
		if out_compression == True:
			df.to_csv(out_path + ".gz", compression="gzip", sep=sep)
		else:
			df.to_csv(out_path, sep=sep)

	return df


def make_map_reduce_jobs(scripts_dir, prep, map_task_l, reduce_task, map_task_exec_l=[], mem="14G"):

	'''Delete existing .sh files in the scripts dir.'''
	for f in glob.glob(os.path.join(scripts_dir,"*.sh")):
    	os.remove(f)

	'''Make map.sh'''
	map_file = os.path.join(scripts_dir,".".join([prep,"map","sh"]))
	map_stream = open(map_file,'w')
	map_stream.write(get_sge_job_txt(prep, map_cmd_n=len(map_task_exec_l), mem=mem))
	map_stream.close()

	'''Make reduce.sh'''
	reduce_file = os.path.join(scripts_dir,".".join([prep,"reduce","sh"]))
	reduce_stream = open(reduce_file,'w')
	reduce_stream.write(get_sge_job_txt(prep, reduce_task=reduce_task, mem=mem))
	reduce_stream.close()
	
	'''Make bash task scripts for each map task and a text file containing references to them.'''
	print "# of map tasks to execute: " + str(len(map_task_exec_l))
	all_map_task_stream = open(os.path.join(scripts_dir,".".join([prep,"all_map_task","txt"])), 'w') 	
	for i in xrange(len(map_task_l)):
		task_script = os.path.join(scripts_dir,".".join([prep,str(i+1),"sh"]))
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
	submit_stream = open(os.path.join(scripts_dir,"submit_map_reduce.sh"),'w')
	submit_stream_l = ["#!/bin/bash"]
	submit_stream_l.append("\n")
	if len(map_task_exec_l) == 0:
		submit_stream_l.append(" ".join(["qsub",reduce_file]))
	else:
		submit_stream_l.append(" ".join(["qsub","-N",prep+".map","-cwd",map_file]))
		submit_stream_l.append(" ".join(["qsub","-hold_jid",prep+".map","-cwd",reduce_file]))
	submit_stream.write("\n".join(submit_stream_l))
	submit_stream.close()
	subprocess.call(['chmod', 'u+x', os.path.join(scripts_dir,"submit_map_reduce.sh")])


def get_sge_job_txt(prep, map_cmd_n=None, reduce_task=None, mem="14G"):

	oe_txt_l = ["mkdir -p ${scriptname}.qsub.out ${scriptname}.qsub.err","exec >${scriptname}.qsub.out/${scriptname}_${JOB_ID}.out 2>${scriptname}.qsub.err/${scriptname}_${JOB_ID}.err"]
	sge_job_txt_l = ["#!/bin/bash","#$ -S /bin/bash","#$ -o /dev/null","#$ -e /dev/null","#$ -cwd","#$ -V","#$ -l tmem="+mem+",h_vmem="+mem,"#$ -l h_rt=24:0:0","set -u","set -x","\n"]
	if map_cmd_n != None:
		sge_job_txt_l.insert(-3,"#$ -t 1:"+str(map_cmd_n))
		sge_job_txt_l.extend(["scriptname=" + prep + ".$SGE_TASK_ID"] + oe_txt_l + ['CMDS=`sed -n -e "$SGE_TASK_ID p" ' + prep + '.all_map_task.txt`',"$CMDS"])
	elif reduce_task != None:
		sge_job_txt_l.extend(["scriptname=" + prep + ".reduce"] + oe_txt_l + [reduce_task])

	sge_job_txt_str = "\n".join(sge_job_txt_l)

	return sge_job_txt_str

