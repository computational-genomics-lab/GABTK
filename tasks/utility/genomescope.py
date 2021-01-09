import os
import subprocess
import re
import luigi

from tasks.readCleaning.preProcessReads import cleanFastq
from tasks.readCleaning.reFormatReads import reformat

class GlobalParameter(luigi.Config):
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	projectName = luigi.Parameter()
	domain=luigi.Parameter()
	assembly_name=luigi.Parameter()
	pe_read_dir=luigi.Parameter()

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

genomescope_folder=os.path.join(os.getcwd(),GlobalParameter().projectName,"GenomeAssembly", "GenomeScope")
createFolder(genomescope_folder)

class GenomeScope(luigi.Task):
	projectName = GlobalParameter().projectName
	domain=GlobalParameter().domain
	assembly_name = GlobalParameter().assembly_name
	kmer=luigi.IntParameter(default=21,description="Minimal Kmer Length for assembly. [--kmer 21]")
	ploidy=luigi.IntParameter(default=2,description="Genome Ploidy. [--ploidy 2]")
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	threads=GlobalParameter().threads  
	seq_platforms = luigi.Parameter(default="pe")
	threads=GlobalParameter().threads

	def requires(self):

		if self.seq_platforms == "pe" and self.pre_process_reads=="yes":
			return [cleanFastq(seq_platforms="pe",sampleName=i) 
						for i in [line.strip() for line in 
							open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]


		if self.seq_platforms == "pe" and self.pre_process_reads=="no":
			return [reformat(seq_platforms="pe",sampleName=i)
						for i in [line.strip() for line in 
							open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
					
	def output(self):
		GenomeScope_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"GenomeAssembly", "GenomeScope", GlobalParameter().assembly_name+ "/")
		return {'out': luigi.LocalTarget(GenomeScope_folder, "summary.txt")}

	def run(self):

		GenomeScope_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"GenomeAssembly", "GenomeScope", GlobalParameter().assembly_name+ "/")
		GenomeScope = os.path.join(os.getcwd(), GlobalParameter().projectName,"GenomeAssembly", "GenomeScope"+ "/")
		GenomeScope_log_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"log", "GenomeAssembly", "GenomeScope", GlobalParameter().assembly_name + "/")
		temp=os.path.join(os.getcwd(), GlobalParameter().projectName,"GenomeAssembly", "GenomeScope", "temp")


		if self.pre_process_reads=="yes":

			def genomescope_cleanFastq(samplefile):

				with open(samplefile) as fh:
					sample_name_list = fh.read().splitlines()
					left_read_name_suffix = '_R1.fastq'
					right_read_name_suffix = '_R2.fastq'
					pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","CleanedReads","PE-Reads" + "/")
					left_read_name_list = [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
					right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
					#left_reads = ' '.join(left_read_name_list)
					#right_reads = ' '.join(right_read_name_list)
					pe_reads = left_read_name_list + right_read_name_list

					with open((os.path.join(os.getcwd(),GlobalParameter().projectName,"GenomeAssembly", "GenomeScope",GlobalParameter().assembly_name+".lst")), 'w') as fh:
						fh.writelines("%s\n"  % read for read in pe_reads)

		if self.pre_process_reads=="no":

			def genomescope_reformat(samplefile):

				with open(samplefile) as fh:
					sample_name_list = fh.read().splitlines()
					left_read_name_suffix = '_R1.fastq'
					right_read_name_suffix = '_R2.fastq'
					pe_cleaned_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","VerifiedReads","PE-Reads" + "/")
					left_read_name_list = [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
					right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
					#left_reads = ' '.join(left_read_name_list)
					#right_reads = ' '.join(right_read_name_list)
					pe_reads = left_read_name_list + right_read_name_list

					with open((os.path.join(os.getcwd(),GlobalParameter().projectName,"GenomeAssembly", "GenomeScope",GlobalParameter().assembly_name+".lst")), 'w') as fh:
						fh.writelines("%s\n"  % read for read in pe_reads)


		if self.pre_process_reads=="yes":
			genomescope_cleanFastq((os.path.join(os.getcwd(),"sample_list",'pe_samples.lst')))
		if self.pre_process_reads=="no":
			genomescope_reformat((os.path.join(os.getcwd(),"sample_list",'pe_samples.lst')))



		read_list=os.path.join(os.getcwd(),GlobalParameter().projectName,"GenomeAssembly", "GenomeScope",GlobalParameter().assembly_name+".lst")

		run_kmc = "[ -d  {temp} ] || mkdir -p {temp}; " \
				  "[ -d  {GenomeScope} ] || mkdir -p {GenomeScope}; cd {GenomeScope};" \
				  "kmc -k{kmer} -t{threads}  -m6 -ci1 -cs10000 @{read_list} reads {temp}".format(temp=temp,
				   kmer=self.kmer, threads=GlobalParameter().threads, read_list=read_list, GenomeScope=GenomeScope)

		run_hist= "[ -d  {GenomeScope} ] || mkdir -p {GenomeScope}; cd {GenomeScope};" \
				  "kmc_tools transform reads histogram reads.histo -cx10000".format(GenomeScope=GenomeScope)

		run_estimate="[ -d  {GenomeScope} ] || mkdir -p {GenomeScope}; cd {GenomeScope};" \
					"genomescope.R -i reads.histo -o {GenomeScope_folder} -k{kmer} -p{ploidy} ".format(GenomeScope=GenomeScope,
						GenomeScope_folder=GenomeScope_folder,kmer=self.kmer,ploidy=self.ploidy)



		print("****** NOW RUNNING COMMAND ******: " + run_kmc)
		run_cmd(run_kmc)

		print("****** NOW RUNNING COMMAND ******: " + run_hist)
		run_cmd(run_hist)

		print("****** NOW RUNNING COMMAND ******: " + run_estimate)
		run_cmd(run_estimate)

		'''
		command = ("kmc -k21 -t4 -m6 -ci1 -cs10000 @read_listxx.txt reads tmp".format(tmp=tempPath))
		hist_command = ("kmc_tools transform reads histogram reads.histo -cx10000")
		estimate_command = ("genomescope.R -i reads.histo -o genomescope -k 21")


		command = "[ -d  {kmergenie_folder} ] || mkdir -p {kmergenie_folder}; " \
						   "mkdir -p {kmergenie_log_folder}; cd {kmergenie_folder}; " \
						   "/usr/bin/time -v kmergenie {kmergenie_folder}{assembly_name}.lst " \
						   "2>&1 | tee {kmergenie_log_folder}kmergenni.log " \
			.format(kmergenie_folder=kmergenie_folder,assembly_name=GlobalParameter().assembly_name,
					kmergenie_log_folder=kmergenie_log_folder)

		print("Estimating Optimal Kmers using kmergenie")
		print("-----------------------------------------")
		print("Command: ", command)
		print("\n")

		proc = subprocess.Popen(command,
								shell=True,
								stdin=subprocess.PIPE,
								stdout=subprocess.PIPE,
								stderr=subprocess.PIPE)

		stdout_value, stderr_value = proc.communicate()

		parse_output = stdout_value.strip().decode("utf-8")

		p = re.compile(r'^best k:\s+(\d{2})', re.M)
		optimal_k = ' '.join(re.findall(p, parse_output))

		return optimal_k



class GlobalParameter(luigi.Config):
	pairedendDir = luigi.Parameter()
	matepairDir = luigi.Parameter()
	singleendDir = luigi.Parameter()
	longreadDir=luigi.Parameter()
	suffix = luigi.Parameter()
	lrsuffix=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()

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


kmergenie = os.path.join(os.getcwd(), "GenomeAssembly", "KmerGenie" + "/")
createFolder(kmergenie)



def kmergenie_formater(samplefile):
	kmergenie_folder = os.path.join(os.getcwd(), "GenomeAssembly", "KmerGenie" + "/")

	with open(samplefile) as fh:
		sample_name_list = fh.read().splitlines()
		left_read_name_suffix = '_R1.fastq'
		right_read_name_suffix = '_R2.fastq'
		pe_cleaned_read_folder = os.path.join(os.getcwd(), "ReadCleaning", "Cleaned_PE_Reads" + "/")
		left_read_name_list = [ pe_cleaned_read_folder + x + left_read_name_suffix for x in sample_name_list]
		right_read_name_list =[pe_cleaned_read_folder + x + right_read_name_suffix for x in sample_name_list]
		#left_reads = ' '.join(left_read_name_list)
		#right_reads = ' '.join(right_read_name_list)
		pe_reads = left_read_name_list + right_read_name_list

		with open((os.path.join(os.getcwd(),"GenomeAssembly", "KmerGenie","kmergenni_pe.lst")), 'w') as fh:
			fh.writelines("%s\n"  % read for read in pe_reads)


class kmergenie(luigi.Task):
	projectName = luigi.Parameter(default="GenomeAssembly")


	def requires(self):
		return []

	def output(self):
		kmergenie_folder = os.path.join(os.getcwd(), "GenomeAssembly", "KmerGenie" + "/")
		return {'out': luigi.LocalTarget(kmergenie_folder + "histograms.dat")}

	def run(self):
		kmergenie_folder = os.path.join(os.getcwd(), "GenomeAssembly", "KmerGenie" + "/")

		kmergenie_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "KmerGenie" +  "/")


		pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")

		kmergenie_formater(pe_sample_list)


		kmergenie_pe_cmd = "[ -d  {kmergenie_folder} ] || mkdir -p {kmergenie_folder}; " \
						"mkdir -p {kmergenie_log_folder}; cd {kmergenie_folder}; " \
						"kmergenie {kmergenie_folder}kmergenni_pe.lst " \
						"2>&1 | tee {kmergenie_log_folder}kmergenni.log " \
			.format(kmergenie_folder=kmergenie_folder,
					kmergenie_log_folder=kmergenie_log_folder)

		print("****** NOW RUNNING COMMAND ******: " + kmergenie_pe_cmd)
		run_cmd(kmergenie_pe_cmd)
'''



