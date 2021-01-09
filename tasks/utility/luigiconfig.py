#! /usr/bin/env python3
import os
import re
import shutil
from collections import OrderedDict
import optparse
from sys import exit
import psutil
import subprocess
import pandas as pd
import luigi

def run_cmd(cmd):
	p = subprocess.Popen(cmd, bufsize=-1,
				 shell=True,
				 universal_newlines=True,
				 stdout=subprocess.PIPE,
				 executable='/bin/bash')
	output = p.communicate()[0]
	return output

def createFolder(directory):
	try:
		if not os.path.exists(directory):
			os.makedirs(directory)
	except OSError:
		print ('Error: Creating directory. ' + directory)


totalcpus = psutil.cpu_count()
threads = str(totalcpus -1)

mem = psutil.virtual_memory()
maxMemory= int((mem.available/1073741824) -1)

class configureProject(luigi.Task):
	
	cpus = luigi.Parameter(default=f'{threads}')
	maxMemory=luigi.Parameter(default=f'{maxMemory}')
	domain=luigi.Parameter()
	rawDataDir=luigi.Parameter()
	genomeName=luigi.Parameter()
	genomeSize=luigi.Parameter()
	projectName=luigi.Parameter()
	schedulerPort=luigi.Parameter(default="8082")
	userEmail=luigi.Parameter()
	symLinkDir=luigi.Parameter(default="raw_data_symlink")
	
	def require(self):
		return[]

	def output(self):
		current_folder=os.path.join(os.getcwd())
		return {'out': luigi.LocalTarget(current_folder + ".luigi.cfg.tmp")}

	def run (self):
		current_folder=os.path.join(os.getcwd())
		excludeDir=os.path.abspath(os.path.join(symLinkDir, "ex" ))
		genomeDir=os.path.abspath(os.path.join(symLinkDir, "geno" ))
		adapter=os.path.abspath(os.path.join(os.getcwd(),"tasks","utility",'adapters.fasta.gz'))
		paired_end_read_dir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.symLinkDir, "pe" ))
		mate_pair_read_dir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.symLinkDir,"mp"))
		pac_read_dir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.symLinkDir,"pac"))
		ont_read_dir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.symLinkDir,"ont"))
		projectDir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.projectName))
		rawDataDir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.rawDataDir))
		symLinkDir=os.path.abspath(os.path.join(os.path.join(os.getcwd()),self.symLinkDir))

		createFolder(symLinkDir)
		createFolder(excludeDir)
		createFolder(genomeDir)
		createFolder(paired_end_read_dir)
		createFolder(mate_pair_read_dir)
		createFolder(pac_read_dir)
		createFolder(ont_read_dir)
		createFolder(projectDir)
		createFolder(projectDir)
		createFolder("config")

		if os.path.isdir(rawDataDir):
			files = [f for f in os.listdir(rawDataDir) if os.path.isfile(os.path.join(rawDataDir, f))]
			keys = []
			fileList = re.compile(r'^(.+?).(gff|gtf|fasta|fna|ffn|fa|fastq|fq|fastq\.gz|fq\.gz)?$')
			for file in files:
				if fileList.search(file):
					keys.append(file)

		dicts = OrderedDict ()
		#keys = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

		for i in keys:
			
			accepted_values="pe mp ont pac ex geno".split()
			
			val = str(input("Enter Data Type of {data}: \tchoose from [pe:paired-end, mp: mate-pair, ont:nanopore, pac:pacbio geno:reference assembly, ex: exclude file]: ".format(data=i)))
			if val in accepted_values:
				dicts[i] = val
			else:
				print(f'{val} is not a valid option. \tchoose from [pe mp ont pac ex geno]: ')
				val = str(input("Enter Data Type of {data}: \tchoose from [pe mp ont pac ex geno]: ".format(data=i)))

		#print(dicts)

		for key, val in dicts.items():
			if not os.path.exists(os.path.join(symLinkDir, val)):
				os.mkdir(os.path.join(symLinkDir, val))

		##ln -nsf method
		for key, val in dicts.items():
			dest = (os.path.join(symLinkDir,val,key))
			src = (os.path.join(inputDataDir,key))
			source = os.path.abspath(src)
			destination = os.path.abspath(dest)
			escape="\'"
			print("Source:\t {skip}{src}{skip}".format(skip=escape,src=source))
			print("Destination:\t {skip}{dest}{skip}".format(skip=escape,dest=destination))
			#print("Desitnation:", '\'destination\'')

			link_cmd = "ln -nsf "
			create_symlink = "{link_cmd} {source} {destination}".format(link_cmd=link_cmd,source=source,destination=destination)
			print("****** NOW RUNNING COMMAND ******: " + create_symlink)
			print (run_cmd(create_symlink))

		###########################################
		def paired_end_samples(pe_dir):
			pe_read_list=os.listdir(paired_end_read_dir)

			sample_list=[]
			for read in pe_read_list:
				pe_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
				if any (read.endswith(ext) for ext in pe_allowed_extn):
					
					sample_name=read.split("_R1.f",1)[0]
					sample_list.append(sample_name)

					file_extn=read.split('.',1)[1]

			with open ((os.path.join(os.getcwd(),"config",'pe_samples.lst')),'w') as file:
				for sample in sample_list:
					file.write("%s\n" % sample)
			file.close()

			return file_extn

		#############################################
		def mate_pair_samples(mp_dir):
			mp_read_list=os.listdir(mate_pair_read_dir)

			sample_list=[]
			for read in mp_read_list:
				mp_allowed_extn=["_R1.fq","_R1.fastq","_R1.fq.gz","_R1.fastq.gz"]
				if any (read.endswith(ext) for ext in mp_allowed_extn):
					
					sample_name=read.split("_R1.f",1)[0]
					sample_list.append(sample_name)

					file_extn=read.split('.',1)[1]

			with open ((os.path.join(os.getcwd(),"config",'mp_samples.lst')),'w') as file:
				for sample in sample_list:
					file.write("%s\n" % sample)
			file.close()

			return file_extn

		#################################################################################

		def pacbio_samples(pb_dir):
			raw_pb_read_list=os.listdir(pac_read_dir)
			sample_list=[]
			for read in raw_pb_read_list:
				raw_pb_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
				
				if any (read.endswith(ext) for ext in raw_pb_allowed_extn):
					
					sample_name=read.split(".",1)[0]
					sample_list.append(sample_name)

					file_extn=read.split('.',1)[1]


			with open ((os.path.join(os.getcwd(),"config",'pac_samples.lst')),'w') as file:
				for sample in sample_list:
					file.write("%s\n" % sample)
			file.close()

			return file_extn

		################################################################################
		def ont_samples(ont_raw_dir):
			corr_ont_read_list=os.listdir(ont_read_dir)
			sample_list=[]
			for read in corr_ont_read_list:
				corr_ont_allowed_extn=[".fq",".fastq",".fq.gz",".fastq.gz"]
				
				if any (read.endswith(ext) for ext in corr_ont_allowed_extn):
					
					sample_name=read.split(".",1)[0]
					sample_list.append(sample_name)
					file_extn=read.split('.',1)[1]


			with open ((os.path.join(os.getcwd(),"config",'ont_samples.lst')),'w') as file:
				for sample in sample_list:
					file.write("%s\n" % sample)
			file.close()

			return file_extn

		#################################################################################
		if not (len(os.listdir(genomeDir))==0):
			genome=genome(genomeDir)
			genomeName=genome[0]
			genomeExtn=genome[1]
		else:
			genomeName="NA"
			genomeExtn="NA"
		
		
		if (len(os.listdir(paired_end_read_dir))!=0):
			paired_end_read_suffix=paired_end_samples(paired_end_read_dir)
		else:
			paired_end_read_suffix="NA"

		if ((len(os.listdir(mate_pair_read_dir))!=0)):
			matepair_read_suffix=mate_pair_samples(mate_pair_read_dir)
		else:
			matepair_read_suffix="NA"

		if ((len(os.listdir(ont_read_dir))!=0)):
			ont_read_suffix=ont_samples(ont_read_dir)
		else:
			ont_read_suffix="NA"

		if ((len(os.listdir(pac_read_dir))!=0)):
			pac_read_suffix=pacbio_samples(pac_read_dir)
		else:
			pac_read_suffix="NA"

			

		with open('.luigi.cfg.tmp', 'w') as config:
			config.write('[core]\n')
			config.write('default-scheduler-port:{port}\n'.format(port=self.schedulerPort))
			config.write('error-email={email}\n\n'.format(email=self.userEmail))
			config.write('[GlobalParameter]\n')
			config.write('projectName={projectName}\n'.format(projectName=self.projectName))
			config.write('assembly_name={assembly_name}\n'.format(assembly_name=self.genomeName))
			config.write('projectDir={projectDir}/\n'.format(projectDir=projectDir))
			config.write('domain={domain}\n'.format(domain=self.domain))
			config.write('adapter={adapter}\n'.format(adapter=adapter))
			config.write('genome_size={genome_size}\n'.format(genome_size=self.genomeSize))

			config.write('genome_dir={genomeDir}/\n'.format(genomeDir=genomeDir))
			config.write('genome_name={genomeName}\n'.format(genomeName=genomeName))
			config.write('genome_suffix={genomeExtn}\n'.format(genomeExtn=genomeExtn))

			
			config.write('pe_read_dir={paired_end_read_dir}/\n'.format(paired_end_read_dir=paired_end_read_dir))
			config.write('pe_read_suffix={paired_end_read_suffix}\n'.format(paired_end_read_suffix=paired_end_read_suffix))

			config.write('mp_read_dir={mate_pair_read_dir}/\n'.format(mate_pair_read_dir=mate_pair_read_dir))
			config.write('mp_read_suffix={matepair_read_suffix}\n'.format(matepair_read_suffix=matepair_read_suffix))

			config.write('ont_read_dir={ont_read_dir}/\n'.format(ont_read_dir=ont_read_dir))
			config.write('ont_read_suffix={ont_read_suffix}\n'.format(ont_read_suffix=ont_read_suffix))

			config.write('pac_read_dir={pac_read_dir}/\n'.format(pac_read_dir=pac_read_dir))
			config.write('pac_read_suffix={pac_read_suffix}\n'.format(pac_read_suffix=pac_read_suffix))


			#if all( [cond1 == 'val1', cond2 == 'val2', cond3 == 'val3', cond4 == 'val4'] ):

			#PE READ
			if all ( [(len(os.listdir(paired_end_read_dir))!=0),
					  (len(os.listdir(ont_read_dir))==0),
					  (len(os.listdir(pac_read_dir))==0),
					  (len(os.listdir(mate_pair_read_dir))==0)]):
				config.write('seq_platforms=pe\n')
				
				
			if all ([(len(os.listdir(ont_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))==0),
				(len(os.listdir(paired_end_read_dir))==0),
				(len(os.listdir(pac_read_dir))==0)]):

				config.write('seq_platforms=ont\n')

				
			if all ([(len(os.listdir(pac_read_dir))!=0),
					(len(os.listdir(ont_read_dir))==0),
					(len(os.listdir(mate_pair_read_dir))==0),
					(len(os.listdir(paired_end_read_dir))==0)]):
				config.write('seq_platforms=pac\n')

				
			#PE and MP RAED
			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))!=0),
				(len(os.listdir(ont_read_dir))==0),
				(len(os.listdir(pac_read_dir))==0)]):
				config.write('seq_platforms=pe-mp\n')

			#PE MP ONT
			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))!=0),
				(len(os.listdir(ont_read_dir))!=0),
				(len(os.listdir(pac_read_dir))==0)]):
				
				config.write('seq_platforms=pe-mp-ont\n')

		
			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))!=0),
				(len(os.listdir(pac_read_dir))!=0),
				(len(os.listdir(ont_read_dir))==0)]):
				
				config.write('seq_platforms=pe-mp-pac\n')


			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(ont_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))==0),
				(len(os.listdir(pac_read_dir))==0)]):
				
				config.write('seq_platforms=pe-ont\n')		

			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(pac_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))==0),
				(len(os.listdir(ont_read_dir))==0)]):
				
				config.write('seq_platforms=pe-pac\n')


			if ([(len(os.listdir(paired_end_read_dir))!=0),
				(len(os.listdir(mate_pair_read_dir))!=0),
				(len(os.listdir(pac_read_dir))!=0),
				(len(os.listdir(ont_read_dir))!=0)]):
				
				config.write('seq_platforms=pe-mp-pac-ont\n')

			
			config.write('threads={cpus}\n'.format(cpus=self.cpus)) 
			config.write('maxMemory={memory}\n'.format(memory=self.maxMemory))
			config.close()

			print("the luigi config file generated")

			rename_config_cmd="mv .luigi.cfg.tmp luigi.cfg"
			print(run_cmd(rename_config_cmd))

if __name__ == "__main__":
	luigi.run()

