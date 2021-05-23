import luigi
import os
import subprocess
from tasks.readCleaning.preProcessReads import cleanFastq
from tasks.readCleaning.preProcessReads import filtlong

class GlobalParameter(luigi.Config):
	genome_size=luigi.Parameter()
	threads = luigi.Parameter()
	maxMemory = luigi.Parameter()
	assembly_name = luigi.Parameter()
	projectName = luigi.Parameter()

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


class canu(luigi.Task):
	projectName = GlobalParameter().projectName
	pre_process_reads = luigi.ChoiceParameter(choices=["yes", "no"], var_type=str)
	assembly_name = GlobalParameter().assembly_name
	threads=GlobalParameter().threads  

	seq_platform = luigi.ChoiceParameter(description="Choose From['pacbio raw, pacbio corrected, nanopore raw, "
												"nanopore corrected']",
									choices=["pacbio-raw", "pacbio-corr", "nano-raw", "nano-corr"], var_type=str)

	genome_size = GlobalParameter().genome_size

	threads=GlobalParameter().threads



	def requires(self):
		if self.seq_platform=="pacbio-raw" and self.pre_process_reads=="yes":
			return [filtlong(seq_platforms="pac",sampleName=i)
							for i in [line.strip() for line in
								open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]

		if self.seq_platform=="nano-raw" and self.pre_process_reads=="yes":
			return [filtlong(seq_platforms="ont",sampleName=i)
							for i in [line.strip() for line in
								open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]

		if self.seq_platform=="pacbio-corr" and self.pre_process_reads=="no":
			return [reformat(seq_platforms="pac",sampleName=i)
							for i in [line.strip()for line in
								open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]

		if self.seq_platform=="nano-corr" and self.pre_process_reads=="no":
			return [reformat(seq_platforms="ont",sampleName=i)
							for i in [line.strip()for line in
								open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]


	

	def output(self):
		canu_assembly_folder = os.path.join(os.getcwd(), self.projectName,"GenomeAssembly", "CANU_" + self.seq_platform, self.assembly_name + "/")
		return {'out': luigi.LocalTarget(canu_assembly_folder,"contigs.fasta")}

	def run(self):
		canu_assembly_folder = os.path.join(os.getcwd(), self.projectName,"GenomeAssembly", "CANU_" + self.seq_platform, self.assembly_name + "/")
		
		canu_assembly_log_folder = os.path.join(os.getcwd(), "log", "GenomeAssembly", "canu_" + self.seq_platform +  "/")		
		

		ont_sample_list = os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")
		pac_sample_list = os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")

		verified_ont_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "ONT-Reads" + "/")
		verified_pac_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "PAC-Reads" + "/")
		cleaned_ont_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "ONT-Reads" + "/")
		cleaned_pac_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "PAC-Reads" + "/")


		def canu_formater(lrfile,inputDir):
			with open(lrfile) as fh:
				sample_name_list = fh.read().splitlines()
				read_name_suffix = '.fastq'
				read_name_list = [inputDir + x + read_name_suffix for x in sample_name_list]
				lr_parse_string = ' '.join(read_name_list)
				return lr_parse_string

		if self.pre_process_reads=="yes" and self.seq_platform=="pacbio-raw":
			input_long_reads = canu_formater(pac_sample_list, cleaned_pac_read_folder)
			platform="-pacbio-raw"
			#processing='-corrected'

		if self.pre_process_reads=="yes" and self.seq_platform=="nano-raw":
			input_long_reads = canu_formater(pac_sample_list, cleaned_ont_read_folder)
			platform="-nanopore-raw"
			#processing='-corrected'


		if self.pre_process_reads=="no" and self.seq_platform=="nano-corr":
			input_long_reads = canu_formater(pac_sample_list, verified_ont_read_folder)
			platform="-nanopore"
			#processing='-corrected'

		if self.pre_process_reads=="no" and self.seq_platform=="pacbio-corr":
			input_long_reads = canu_formater(pac_sample_list, verified_pac_read_folder)
			platform="-pacbio"
			#processing='-corrected'



		canu_cmd = "[ -d  {canu_assembly_folder} ] || mkdir -p {canu_assembly_folder}; " \
						"mkdir -p {canu_assembly_log_folder}; cd {canu_assembly_folder}; " \
						"/usr/bin/time -v canu -correct " \
						"-p {assembly_name} " \
						"-d {canu_assembly_folder} " \
						"genomeSize={genome_size} " \
						"maxThreads={threads} " \
						"{platform} " \
						"{input_long_reads} " \
						"2>&1 | tee {canu_assembly_folder}canu_assembly.log " \
			.format(canu_assembly_folder=canu_assembly_folder,
					assembly_name=self.assembly_name,
					genome_size=self.genome_size,
					platform=platform,
					input_long_reads=input_long_reads,
					canu_assembly_log_folder=canu_assembly_log_folder,
					threads=GlobalParameter().threads)

		print("****** NOW RUNNING COMMAND ******: " + canu_cmd)
		run_cmd(canu_cmd)


 
