import luigi
import os
import subprocess
from tasks.readCleaning.preProcessReads import cleanFastq
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat
from tasks.readCleaning.reFormatReads import reformatReads


class GlobalParameter(luigi.Config):
	assembly_name=luigi.Parameter()
	projectName=luigi.Parameter()
	pac_read_dir=luigi.Parameter()
	pac_read_suffix=luigi.Parameter()
	genome_size=luigi.Parameter()
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


class correctONT(luigi.Task):
	projectName = GlobalParameter().projectName
	seq_platforms = luigi.Parameter(default="ont")
	genome_size = GlobalParameter().genome_size
	assembly_name = GlobalParameter().assembly_name
	threads=GlobalParameter().threads
	pre_process_reads = luigi.ChoiceParameter(choices=["yes","no"],var_type=str)

	
	def requires(self):
		if self.pre_process_reads=="yes":
			return [filtlong(seq_platforms="ont",sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]


		if self.pre_process_reads=="no":
			return [reformat(seq_platforms="ont",sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]


	def output(self):
		necat_correct_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","CorrectReads","ONT-Reads" + "/")
		return {'out': luigi.LocalTarget(necat_correct_folder,self.assembly_name+"_corrected.fasta")}

	def run(self):
		necat_correct_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","CorrectReads","ONT-Reads" + "/")
		read_correction_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "CorrectReads", "ONT-Reads"+ "/")

		necat_assembly_folder=os.path.join(os.getcwd(), self.projectName,"GenomeAssembly", "NECAT" + "/")
		necat_corrected_fasta = os.path.join(os.getcwd(), self.projectName,"GenomeAssembly","NECAT",self.assembly_name,"1-consensus","cns_final.fasta"+"/")

		if self.pre_process_reads=="no":
			long_clean_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","VerifiedReads","ONT-Reads" + "/")
			
		if self.pre_process_reads=="yes":
			long_clean_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","CleanedReads","ONT-Reads" + "/")


				
		def necat_config_generator(lrfile):
			if self.pre_process_reads=="no":
				long_clean_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","VerifiedReads","ONT-Reads" + "/")
			
			if self.pre_process_reads=="yes":
				long_clean_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName,"ReadQC","CleanedReads","ONT-Reads" + "/")


			createFolder(long_clean_read_folder)						

			long_read_list= os.path.join(long_clean_read_folder, "long_read.lst")

			with open(lrfile) as fh:
				
				with open(long_read_list, "w") as fo:
					sample_name_list = fh.read().splitlines()
					read_name_suffix = '.fastq'
					read_name_list = [long_clean_read_folder + x + read_name_suffix for x in sample_name_list]
					lr_parse_string = ' '.join(read_name_list)
					fo.write(lr_parse_string)
		
		assembly_name=self.assembly_name
		lr_sample_list = os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")
		necat_config_generator(lr_sample_list)

		genome_size=self.genome_size
		threads=GlobalParameter().threads
		config_folder_path=os.path.join(os.getcwd(), "configuration")
		necat_correct_config_path=os.path.join(os.getcwd(), "configuration", "necat_correct.config")
		createFolder(config_folder_path)
		
		with open (necat_correct_config_path,'w') as configfile:
			configfile.write('PROJECT={assembly_name}\n'.format(assembly_name=self.assembly_name))
			configfile.write('ONT_READ_LIST={long_clean_read_folder}long_read.lst\n'.format(long_clean_read_folder=long_clean_read_folder))
			configfile.write('GENOME_SIZE={genome_size}\n'.format(genome_size=genome_size))
			configfile.write('THREADS={threads}\n'.format(threads=GlobalParameter().threads))
			configfile.write('MIN_READ_LENGTH=3000/\n')
			configfile.write('OVLP_FAST_OPTIONS="-n 500 -z 20 -b 2000 -e 0.5 -j 0 -u 1 -a 1000"\n')
			configfile.write('OVLP_SENSITIVE_OPTIONS="-n 500 -z 10 -e 0.5 -j 0 -u 1 -a 1000"\n')
			configfile.write('CNS_FAST_OPTIONSC"-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0"\n')
			configfile.write('CNS_SENSITIVE_OPTIONS="-a 2000 -x 4 -y 12 -l 1000 -e 0.5 -p 0.8 -u 0"\n')
			configfile.write('TRIM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 1 -a 400"\n')
			configfile.write('ASM_OVLP_OPTIONS="-n 100 -z 10 -b 2000 -e 0.5 -j 1 -u 0 -a 400"\n')
			configfile.write('NUM_ITER=2\n')
			configfile.write('CNS_OUTPUT_COVERAGE=45\n')
			configfile.write('CLEANUP=0\n')
			configfile.write('USE_GRID=false\n')
			configfile.write('GRID_NODE=0\n')
			configfile.write('FSA_OL_FILTER_OPTIONS="--max_overhang=-1 --min_identity=-1 --coverage=40"\n')
			configfile.write('FSA_ASSEMBLE_OPTIONS=""\n')
			configfile.write('FSA_CTG_BRIDGE_OPTIONS="--dump --read2ctg_min_identity=80 --read2ctg_min_coverage=4 --read2ctg_max_overhang=500 --read_min_length=5000 --ctg_min_length=1000 --read2ctg_min_aligned_length=5000 --select_branch=best"\n')


		necat_correct_cmd = "[ -d  {necat_assembly_folder} ] || mkdir -p {necat_assembly_folder}; " \
						"mkdir -p {read_correction_log_folder}; cd {necat_assembly_folder}; " \
						"/usr/bin/time -v necat.pl correct {necat_correct_config_path} " \
						"2>&1 | tee {read_correction_log_folder}necat_assembly.log " \
			.format(necat_assembly_folder=necat_assembly_folder,necat_correct_config_path=necat_correct_config_path,
					read_correction_log_folder=read_correction_log_folder)

		print("****** NOW RUNNING COMMAND ******: " + necat_correct_cmd)
		run_cmd(necat_correct_cmd)


		cp_corrected_fasta_cmd="[ -d  {necat_correct_folder} ] || mkdir -p {necat_correct_folder}; " \
							   "cp {necat_corrected_fasta} {necat_correct_folder}".format(necat_correct_folder=necat_correct_folder,necat_corrected_fasta=necat_corrected_fasta)

		print("****** NOW RUNNING COMMAND ******: " + cp_corrected_fasta_cmd)
		run_cmd(cp_corrected_fasta_cmd)