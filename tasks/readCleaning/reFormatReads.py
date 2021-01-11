import luigi
import time
import os
import subprocess
from tasks.readCleaning.cleanedReadQC import *

class GlobalParameter(luigi.Config):
	pe_read_dir=luigi.Parameter()	
	mp_read_dir=luigi.Parameter()
	pac_read_dir=luigi.Parameter()
	ont_read_dir=luigi.Parameter()
	pe_read_suffix=luigi.Parameter()		
	mp_read_suffix=luigi.Parameter()
	pac_read_suffix=luigi.Parameter()
	ont_read_suffix=luigi.Parameter()
	projectName=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	adapter = luigi.Parameter()
	

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

createFolder("task_logs")

class reformat(luigi.Task):
	paired_end_read_dir = GlobalParameter().pe_read_dir
	mate_pair_read_dir = GlobalParameter().mp_read_dir
	nanopore_read_dir = GlobalParameter().ont_read_dir
	pacbio_read_dir = GlobalParameter().pac_read_dir

	paired_end_read_suffix = GlobalParameter().pe_read_suffix
	mate_pair_read_suffix = GlobalParameter().mp_read_suffix
	nanopore_read_suffix = GlobalParameter().ont_read_suffix
	pacbio_read_suffix = GlobalParameter().pac_read_suffix

	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	projectName = GlobalParameter().projectName

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")

	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)

	

	def output(self):
		pe_verified_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","VerifiedReads","PE-Reads" + "/")
		mp_verified_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","VerifiedReads","MP-Reads" + "/")
		ont_verified_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","VerifiedReads","ONT-Reads" + "/")
		pac_verified_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","VerifiedReads","PAC-Reads" + "/")


		###############################################################################################
		if self.seq_platforms == "pe":
			return {'out1': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "ont":
			return {'out1': luigi.LocalTarget(ont_verified_read_folder + self.sampleName + ".fastq")}

		if self.seq_platforms == "pac":
			return {'out1': luigi.LocalTarget(pac_verified_read_folder + self.sampleName + ".fastq")}


		if self.seq_platforms == "mp":
			return {'out1': luigi.LocalTarget(mp_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_verified_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "pe-mp":
			return {'out1': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_verified_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "pe-ont":
			return {'out1': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(ont_verified_read_folder + self.sampleName + ".fastq")
					}

		if self.seq_platforms == "pe-pac":
			return {'out1': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_verified_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(pac_verified_read_folder + self.sampleName + ".fastq")
					}

		
	def run(self):
		pe_verified_read_folder = os.path.join(os.getcwd(),self.projectName,"ReadQC","VerifiedReads","PE-Reads" + "/")
		mp_verified_read_folder = os.path.join(os.getcwd(),self.projectName,"ReadQC","VerifiedReads","MP-Reads" + "/")
		ont_verified_read_folder = os.path.join(os.getcwd(),self.projectName,"ReadQC","VerifiedReads","ONT-Reads" + "/")
		pac_verified_read_folder = os.path.join(os.getcwd(),self.projectName,"ReadQC","VerifiedReads","PAC-Reads" + "/")


		pe_verification_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "VerifiedReads" ,"PE-Reads" + "/")
		mp_verification_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "VerifiedReads","MP-Reads" + "/")
		ont_verification_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "VerifiedReads" ,"ONT-Reads"+ "/")
		pac_verification_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "VerifiedReads" ,"PAC-Reads" + "/")

		
		cmd_verify_pe ="[ -d  {pe_verified_read_folder} ] || mkdir -p {pe_verified_read_folder}; mkdir -p {pe_verification_log_folder}; " \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "tossbrokenreads=t " \
					   "verifypaired=t " \
					   "in1={pe_read_dir}{sampleName}_R1.{pe_read_suffix} " \
					   "in2={pe_read_dir}{sampleName}_R2.{pe_read_suffix} " \
					   "out={pe_verified_read_folder}{sampleName}_R1.fastq " \
					   "out2={pe_verified_read_folder}{sampleName}_R2.fastq " \
					   " 2>&1 | tee {pe_verification_log_folder}{sampleName}_pe_reformat_run.log "\
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					pe_read_dir=GlobalParameter().pe_read_dir,
					pe_read_suffix=GlobalParameter().pe_read_suffix,
					sampleName=self.sampleName,
					pe_verified_read_folder=pe_verified_read_folder,
					pe_verification_log_folder=pe_verification_log_folder)

		##################
		cmd_verify_mp = "[ -d  {mp_verified_read_folder} ] || mkdir -p {mp_verified_read_folder}; mkdir -p {mp_verification_log_folder}; " \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "verifypaired=t " \
					   "tossbrokenreads=t " \
					   "in1={mp_read_dir}{sampleName}_R1.{mp_read_suffix} " \
					   "in2={mp_read_dir}{sampleName}_R2.{mp_read_suffix} " \
					   "out={mp_verified_read_folder}{sampleName}_R1.fastq " \
					   "out2={mp_verified_read_folder}{sampleName}_R2.fastq " \
					   " 2>&1 | tee {mp_verification_log_folder}{sampleName}_pe_reformat_run.log " \
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					mp_read_dir=self.mate_pair_read_dir,
					mp_read_suffix=self.mate_pair_read_suffix,
					sampleName=self.sampleName,
					mp_verified_read_folder=mp_verified_read_folder,
					mp_verification_log_folder=mp_verification_log_folder)
		##################


		cmd_verify_ont = "[ -d  {ont_verified_read_folder} ] || mkdir -p {ont_verified_read_folder}; mkdir -p {ont_verification_log_folder};" \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "tossbrokenreads=t " \
					   "in1={ont_read_dir}{sampleName}.{ont_read_suffix} " \
					   "out={ont_verified_read_folder}{sampleName}.fastq " \
					   " 2>&1 | tee {ont_verification_log_folder}{sampleName}_reformat_run.log " \
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					ont_read_dir=GlobalParameter().ont_read_dir,
					ont_read_suffix=GlobalParameter().ont_read_suffix,
					sampleName=self.sampleName,
					ont_verified_read_folder=ont_verified_read_folder,
					ont_verification_log_folder=ont_verification_log_folder)

		cmd_verify_pac = "[ -d  {pac_verified_read_folder} ] || mkdir -p {pac_verified_read_folder}; mkdir -p {pac_verification_log_folder};" \
					   "reformat.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "tossbrokenreads=t " \
					   "in1={pac_read_dir}{sampleName}.{pac_read_suffix} " \
					   "out={pac_verified_read_folder}{sampleName}.fastq " \
					   " 2>&1 | tee {pac_verification_log_folder}{sampleName}_reformat_run.log " \
			.format(Xmx=GlobalParameter().maxMemory,
					cpu=GlobalParameter().threads,
					pac_read_dir=GlobalParameter().pac_read_dir,
					pac_read_suffix=GlobalParameter().pac_read_suffix,
					sampleName=self.sampleName,
					pac_verified_read_folder=pac_verified_read_folder,
					pac_verification_log_folder=pac_verification_log_folder)

		if self.seq_platforms == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))

		if self.seq_platforms == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp)
			print(run_cmd(cmd_verify_mp))

		if self.seq_platforms == "ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_ont)
			print(run_cmd(cmd_verify_ont))

		if self.seq_platforms == "pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pac)
			print(run_cmd(cmd_verify_pac))


		if self.seq_platforms == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_mp)
			print(run_cmd(cmd_verify_mp))

		if self.seq_platforms == "pe-ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_ont)
			print(run_cmd(cmd_verify_ont))

		if self.seq_platforms == "pe-pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pe)
			print(run_cmd(cmd_verify_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_verify_pac)
			print(run_cmd(cmd_verify_pac))





class reformatReads(luigi.Task):

	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio",
                                             choices=["pe", "mp", "pe-mp", "pe-ont", "pe-pac"], var_type=str)

	#seq_platforms=GlobalParameter().seq_platforms
	def requires(self):

		if self.seq_platforms == "pe":
			return [

					[reformat(seq_platforms=self.seq_platforms,sampleName=i)
										for i in [line.strip() for line in
							  					open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
			        ]

		if self.seq_platforms == "mp":
			return [

					[reformat(seq_platforms=self.seq_platforms,sampleName=i)
										for i in [line.strip() for line in
							  					open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
			        ]

		if self.seq_platforms == "pe-mp":

			return [
						[reformat(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[reformat(seq_platforms="mp",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","mp_samples.lst")))]]
				   ]

		if self.seq_platforms == "pe-ont":

			return [
						[reformat(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[reformat(seq_platforms="ont",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","ont_samples.lst")))]]
				   ]

		if self.seq_platforms == "pe-pac":

			return [
						[reformat(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[reformat(seq_platforms="pac",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pac_samples.lst")))]]
				   ]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.validate.reads.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('read validation finished at {t}'.format(t=timestamp))