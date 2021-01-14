import luigi
import time
import os
import subprocess
from tasks.readCleaning.cleanedReadQC import *
from tasks.readCleaning.rawReadQC import readqc

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
	seq_platforms=luigi.Parameter()
	

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

#createFolder("task_logs")

class cleanFastq(luigi.Task):
	pe_read_dir = GlobalParameter().pe_read_dir
	mp_read_dir = GlobalParameter().mp_read_dir
	ont_read_dir = GlobalParameter().ont_read_dir
	pac_read_dir = GlobalParameter().pac_read_dir

	pe_read_suffix = GlobalParameter().pe_read_suffix
	mp_read_suffix = GlobalParameter().mp_read_suffix
	ont_read_suffix = GlobalParameter().ont_read_suffix
	pac_read_suffix = GlobalParameter().pac_read_suffix

	
	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	projectName = GlobalParameter().projectName

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")

	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)


	adapter = GlobalParameter().adapter

	kmer_length=luigi.OptionalParameter(default="21",description="Kmer length used for finding contaminants. Contaminants "
												 "shorter than kmer length will not be found.Default: 21")

	corret_error=luigi.BoolParameter(default=False,description="Perform Error Correction Or Not")

	k_trim = luigi.ChoiceParameter(default="r",description="Trimming protocol to remove bases matching reference kmers from reads. "
											   "Choose From['f: dont trim','r: trim to right','l: trim to left]",
								   choices=["f", "r", "l"], var_type=str)

	quality_trim=luigi.ChoiceParameter(default="lr",description="Trimming protocol to remove bases with quality below the minimum "
												   "average region quality from read ends. Performed after looking for kmers."
												   " If enabled, set also 'Average quality below which to trim region'. "
												   "Choose From['f: trim neither end', 'rl: trim both end','r: trim only right end','l: trim only left end]",
								   choices=["f", "lr", "r","l"], var_type=str)

	trim_quality = luigi.IntParameter(description="Average quality below which to trim region ",default=6)

	

	min_length = luigi.OptionalParameter(default="40",description="reads shorter than min_length will be discarded. Default: "
													"min_length=40")

	min_average_quality = luigi.OptionalParameter(default="10", description="Reads with average quality (after trimming) below "
																"this will be discarded. Default: min_average_quality=10")
	min_base_quality = luigi.OptionalParameter(default="0",description="Reads with any base below this quality (after trimming) will be discarded. "
																"Default: min_base_quality=0")
	min_GC = luigi.OptionalParameter(default="0.0", description="Discard reads with GC content below this. Default: min_gc=0.0")
	max_GC = luigi.OptionalParameter(default="1.0", description="Discard reads with GC content below this. Default: max_gc=1")
	kmer = luigi.OptionalParameter(default="13", description="Kmer length used for finding contaminants. Default: kmer=13")

	trim_front = luigi.Parameter(default="0", description="trimming how many bases in front for read. Default: "
														  "trim_front=0")
	trim_tail = luigi.Parameter(default="0", description="trimming how many bases in tail for read. Default: "
														  "trim_tail=0")
	max_n=luigi.IntParameter(default=-1,description="Maximum number of Ns after trimming [maxns=-1]. "
												   "If non-negative, reads with more Ns than this (after trimming) will be discarded.")

	trim_by_overlap=luigi.Parameter(default="f",description="Trim adapters based on where paired-end reads overlap [tbo]")
	trim_pairs_evenly=luigi.Parameter(default="f", description="Trim both sequences of paired-end reads to the minimum length of either sequence ["
															   "tpe]")
	
	long_read_min_length = luigi.OptionalParameter(default="1000", description="This parameter is specific for long read only ["
																			   "seq_platforms='pac or ont']. Reads shorter than "
																			   "min_length will be discarded. Default: "
																			   "long_read_min_length=1000")

	long_read_mean_quality = luigi.IntParameter(default="80", description="This parameter is specific for long read only. [seq_platforms='pac or ont']"
																		  "The mean quality is the mean read identity as indicated by the Phread Quality Score. Default long_read_mean_quality=80")

	long_read_keep_percent = luigi.IntParameter(default="90", description="This parameter is specific for long read only. [seq_platforms='pac or ont'] "
																		  "The percentage of the best reads to be retained. Default: keep_percent=90")
	
	def output(self):


		pe_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "PE-Reads" + "/")
		mp_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "MP-Reads" + "/")
		ont_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "ONT-Reads" + "/")
		pac_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "PAC-Reads" + "/")
		
		if self.seq_platforms == "pe":
			return {'out1': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "ont":
			return {'out1': luigi.LocalTarget(ont_clean_read_folder + self.sampleName + ".fastq")}

		if self.seq_platforms == "pac":
			return {'out1': luigi.LocalTarget(pac_clean_read_folder + self.sampleName + ".fastq")}


		if self.seq_platforms == "mp":
			return {'out1': luigi.LocalTarget(mp_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(mp_clean_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "pe-mp":
			return {'out1': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(mp_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out4': luigi.LocalTarget(mp_clean_read_folder + self.sampleName + "_R2.fastq")}

		if self.seq_platforms == "pe-ont":
			return {'out1': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(ont_clean_read_folder + self.sampleName + ".fastq")
					}

		if self.seq_platforms == "pe-pac":
			return {'out1': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R1.fastq"),
					'out2': luigi.LocalTarget(pe_clean_read_folder + self.sampleName + "_R2.fastq"),
					'out3': luigi.LocalTarget(pac_clean_read_folder + self.sampleName + ".fastq")
					}

		
	def run(self):
		pe_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","PE-Reads" + "/")
		mp_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","MP-Reads" + "/")

		ont_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "ONT-Reads" + "/")
		pac_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads", "PAC-Reads" + "/")
		
		read_clean_log_folder = os.path.join(os.getcwd(), "log","ReadCleaning" + "/")

		cleanFastq_pe_clean_stat_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "CleanedReads","Cleaned_PE_Reads_STAT" + "/")
		cleanFastq_mp_clean_stat_folder = os.path.join(os.getcwd(), self.projectName, "ReadQC","CleanedReads","Cleaned_MP_Reads_STAT" + "/")

		


		
		cmd_clean_pe = "[ -d  {pe_clean_read_folder} ] || mkdir -p {pe_clean_read_folder}; " \
					   "mkdir -p {cleanFastq_pe_clean_stat_folder}; mkdir -p {read_clean_log_folder}; bbduk.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "ecco={corret_error} " \
					   "minlength={min_length} " \
					   "minavgquality={min_average_quality} " \
					   "minbasequality={min_base_quality} " \
					   "trimq={trim_quality} " \
					   "qtrim={quality_trim} " \
					   "ftl={trim_front} " \
					   "ftr2={trim_tail} " \
					   "mingc={min_GC} " \
					   "maxgc={max_GC} " \
					   "maxns={max_n} " \
					   "tbo={trim_by_overlap} " \
					   "tpe={trim_pairs_evenly} " \
					   "in1={paired_end_read_dir}{sampleName}_R1.{paired_end_read_suffix} " \
					   "in2={paired_end_read_dir}{sampleName}_R2.{paired_end_read_suffix} " \
					   "out={pe_clean_read_folder}{sampleName}_R1.fastq " \
					   "out2={pe_clean_read_folder}{sampleName}_R2.fastq " \
					   "outs={pe_clean_read_folder}{sampleName}.fastq " \
					   "ziplevel=9 " \
					   "ref={adapter} " \
					   "stats={cleanFastq_pe_clean_stat_folder}{sampleName}.stat " \
					   "bqhist={cleanFastq_pe_clean_stat_folder}{sampleName}.qual.hist " \
					   "gchist={cleanFastq_pe_clean_stat_folder}{sampleName}.gc.hist " \
					   " 2>&1 | tee {read_clean_log_folder}{sampleName}_pe_cleanFastq_run.log "\
			.format(Xmx=self.maxMemory,
					cpu=self.threads,
					paired_end_read_dir=self.pe_read_dir,
					paired_end_read_suffix=self.pe_read_suffix,
					sampleName=self.sampleName,
					corret_error=self.corret_error,
					adapter=self.adapter,
					pe_clean_read_folder=pe_clean_read_folder,
					min_length=self.min_length,
					min_average_quality=self.min_average_quality,
					min_base_quality=self.min_base_quality,
					trim_quality=self.trim_quality,
					quality_trim=self.quality_trim,
					trim_front=self.trim_front,
					trim_tail=self.trim_tail,
					min_GC=self.min_GC, max_n=self.max_n,
					max_GC=self.max_GC,
					kmer=self.kmer,
					trim_by_overlap=self.trim_by_overlap,
					trim_pairs_evenly=self.trim_pairs_evenly,
					cleanFastq_pe_clean_stat_folder=cleanFastq_pe_clean_stat_folder,
					read_clean_log_folder=read_clean_log_folder)

		##################
		cmd_clean_mp = "[ -d  {mp_clean_read_folder} ] || mkdir -p {mp_clean_read_folder}; " \
					   "mkdir -p {cleanFastq_mp_clean_stat_folder}; mkdir -p {read_clean_log_folder}; bbduk.sh " \
					   "-Xmx{Xmx}g " \
					   "threads={cpu} " \
					   "ecco={corret_error} " \
					   "minlength={min_length} " \
					   "minavgquality={min_average_quality} " \
					   "minbasequality={min_base_quality} " \
					   "trimq={trim_quality} " \
					   "qtrim={quality_trim} " \
					   "ftl={trim_front} " \
					   "ftr2={trim_tail} " \
					   "mingc={min_GC} " \
					   "maxgc={max_GC} " \
					   "maxns={max_n} " \
					   "tbo={trim_by_overlap} " \
					   "tpe={trim_pairs_evenly} " \
					   "in1={mate_pair_read_dir}{sampleName}_R1.{mate_pair_read_suffix} " \
					   "in2={mate_pair_read_dir}{sampleName}_R2.{mate_pair_read_suffix} " \
					   "out={mp_clean_read_folder}{sampleName}_R1.fastq " \
					   "out2={mp_clean_read_folder}{sampleName}_R2.fastq " \
					   "ziplevel=9 " \
					   "ref={adapter} " \
					   "stats={cleanFastq_mp_clean_stat_folder}{sampleName}.stat " \
					   "bqhist={cleanFastq_mp_clean_stat_folder}{sampleName}.qual.hist " \
					   "gchist={cleanFastq_mp_clean_stat_folder}{sampleName}.gc.hist " \
					   " 2>&1 | tee {read_clean_log_folder}{sampleName}_mp_cleanFastq_run.log " \
			.format(Xmx=self.maxMemory,
					cpu=self.threads,
					mate_pair_read_dir=self.mp_read_dir,
					mate_pair_read_suffix=self.mp_read_suffix,
					sampleName=self.sampleName,
					corret_error=self.corret_error,
					adapter=self.adapter,
					mp_clean_read_folder=mp_clean_read_folder,
					min_length=self.min_length,
					min_average_quality=self.min_average_quality,
					min_base_quality=self.min_base_quality,
					trim_quality=self.trim_quality,quality_trim=self.quality_trim,
					trim_front=self.trim_front,
					trim_tail=self.trim_tail,
					min_GC=self.min_GC, max_n=self.max_n,
					max_GC=self.max_GC,
					kmer=self.kmer,
					trim_by_overlap=self.trim_by_overlap,
					trim_pairs_evenly=self.trim_pairs_evenly,
					cleanFastq_mp_clean_stat_folder=cleanFastq_mp_clean_stat_folder,
					read_clean_log_folder=read_clean_log_folder)
		##################

		
		cmd_clean_ont = "[ -d  {ont_clean_read_folder} ] || mkdir -p {ont_clean_read_folder}; " \
					   "cd {ont_clean_read_folder}; " \
					   "filtlong --min_length {long_read_min_length} " \
					   "--keep_percent {long_read_keep_percent} " \
					   "--min_mean_q {long_read_mean_quality} " \
					   "{ont_read_dir}{sampleName}.{ont_read_suffix} > {ont_clean_read_folder}{sampleName}.fastq " \
			.format(ont_clean_read_folder=ont_clean_read_folder,
					ont_read_suffix=self.ont_read_suffix,
					sampleName=self.sampleName,
					ont_read_dir=GlobalParameter().ont_read_dir,
					long_read_keep_percent=self.long_read_keep_percent,
					long_read_mean_quality=self.long_read_mean_quality,
					long_read_min_length=self.long_read_min_length)



		cmd_clean_pac = "[ -d  {pac_clean_read_folder} ] || mkdir -p {pac_clean_read_folder}; " \
					   "cd {pac_clean_read_folder}; " \
					   "filtlong --min_length {long_read_min_length} " \
					   "--keep_percent {long_read_keep_percent} " \
					   "--min_mean_q {long_read_mean_quality} " \
					   "{pac_read_dir}{sampleName}.{pac_read_suffix} > {pac_clean_read_folder}{sampleName}.fastq " \
			.format(pac_clean_read_folder=pac_clean_read_folder,
					pac_read_suffix=self.pac_read_suffix,
					sampleName=self.sampleName,
					pac_read_dir=GlobalParameter().pac_read_dir,
					long_read_keep_percent=self.long_read_keep_percent,
					long_read_mean_quality=self.long_read_mean_quality,
					long_read_min_length=self.long_read_min_length)


		if self.seq_platforms == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print(run_cmd(cmd_clean_pe))

		if self.seq_platforms == "pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pac)
			print(run_cmd(cmd_clean_pac))
			
		if self.seq_platforms == "ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_ont)
			print(run_cmd(cmd_clean_ont))

		if self.seq_platforms == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_mp)
			print(run_cmd(cmd_clean_mp))


		if self.seq_platforms == "pe-mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print(run_cmd(cmd_clean_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_mp)
			print(run_cmd(cmd_clean_mp))

		if self.seq_platforms == "pe-ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print(run_cmd(cmd_clean_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_ont)
			print(run_cmd(cmd_clean_ont))

		if self.seq_platforms == "pe-pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pe)
			print(run_cmd(cmd_clean_pe))
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pac)
			print(run_cmd(cmd_clean_pac))



class cleanReads(luigi.Task):
	
	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)
	def requires(self):

		if self.seq_platforms == "pe":
			return [

					[cleanFastq(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

					[readqc(seq_platforms=self.seq_platforms,
									sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]
			       ]


		if self.seq_platforms == "mp":
			return [[cleanFastq(seq_platforms=self.seq_platforms,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(),"sample_list", "mp_samples.lst")))]],

					[readqc(seq_platforms=self.seq_platforms,
							sampleName=i)
					 for i in [line.strip()
							   for line in
							   open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
					]


		if self.seq_platforms == "pe-mp":

			return [
						[cleanFastq(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[cleanFastq(seq_platforms="mp", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","mp_samples.lst")))]],

						[readqc(seq_platforms="pe",
									sampleName=i)
								for i in [line.strip()
						   				 for line in
						  						open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

						[readqc(seq_platforms="mp",
									sampleName=i)
				 				for i in [line.strip()
						   				 for line in
						   						open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
				  ]

		if self.seq_platforms == "ont":
			return [
				[cleanFastq(seq_platforms="ont", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]],
				[readqc(seq_platforms="ont",
						sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]
			    ]


		if self.seq_platforms == "pac":
			return [
				[cleanFastq(seq_platforms="pac", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]],
				[readqc(seq_platforms="pac",
						sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]
			    ]

		

		if self.seq_platforms == "pe-ont":

			return [
						[cleanFastq(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],
						[cleanFastq(seq_platforms="ont", sampleName=i)
				 				for i in [line.strip()
						   					for line in
						   						open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]],
						[readqc(seq_platforms="pe",
									sampleName=i)
								for i in [line.strip()
						   				 for line in
						  						open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],
						[readqc(seq_platforms="ont",
									sampleName=i)
				 				for i in [line.strip()
						   				for line in
						  						 open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]
				  ]


		if self.seq_platforms == "pe-pac":

			return [
						[cleanFastq(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],
						[cleanFastq(seq_platforms="pac", sampleName=i)
				 				for i in [line.strip()
						   					for line in
						   						open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]],
						[readqc(seq_platforms="pe",
									sampleName=i)
								for i in [line.strip()
						   				 for line in
						  						open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],
						[readqc(seq_platforms="pac",
									sampleName=i)
				 				for i in [line.strip()
						   				for line in
						  						 open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]
				  ]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.clean.shortread.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('short read processing finished at {t}'.format(t=timestamp))


class filtlong(luigi.Task):

	threads=GlobalParameter().threads
	ont_read_dir = GlobalParameter().ont_read_dir
	ont_read_suffix = GlobalParameter().ont_read_suffix
	pac_read_dir = GlobalParameter().pac_read_dir
	pac_read_suffix = GlobalParameter().pac_read_suffix

	platform=luigi.ChoiceParameter(description="Choose From['ont': nanopore','pac':pacbio]",
                                             choices=["ont","pac"], var_type=str)
	
	projectName = GlobalParameter().projectName

	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	long_read_min_length = luigi.IntParameter(default="1000")
	long_read_mean_quality = luigi.Parameter(default="80", description="The mean quality is the mean read identity as indicated by the "
													 "Phred quality scores. Example: example, consider a read where all the fastq quality characters are +. "
													 "The qscores for each base are 10 which equates to a 90 percentage hance of being correct. "
													 "This read would then have a mean quality score of 90. "
													 "Read mean qualities are converted to a z-score and scaled to the range 0-100 to make the mean quality score. "
													 "This means that the read with the worst mean quality in the input set will get a mean quality score of 0 "
													 "and the read with the best mean quality will get a mean quality score of 100."
													 " Default: meanQ=80")

	long_read_keep_percent = luigi.Parameter(default="90", description="The percentage of the best reads to be retained."
													  " Default: keep_percent=90")

	def output(self):
		ont_clean_read_folder = os.path.join(os.getcwd(), self.projectName, "ReadQC", "CleanedReads", "ONT-Reads" + "/")
		pac_clean_read_folder = os.path.join(os.getcwd(), self.projectName, "ReadQC", "CleanedReads", "PAC-Reads" + "/")
 
		if self.platform == "ont":
			return {'out': luigi.LocalTarget(ont_clean_read_folder + self.sampleName + ".fastq")}
		if self.platform == "pac":
			return {'out': luigi.LocalTarget(pac_clean_read_folder + self.sampleName + ".fastq")}


	def run(self):
		ont_clean_read_folder = os.path.join(os.getcwd(), self.projectName, "ReadQC", "CleanedReads", "ONT-Reads" + "/")
		pac_clean_read_folder = os.path.join(os.getcwd(), self.projectName, "ReadQC", "CleanedReads", "PAC-Reads" + "/")

		cmd_clean_ont = "[ -d  {ont_clean_read_folder} ] || mkdir -p {ont_clean_read_folder}; " \
					   "cd {ont_clean_read_folder}; " \
					   "filtlong --min_length {long_read_min_length} " \
					   "--keep_percent {long_read_keep_percent} " \
					   "--min_mean_q {long_read_mean_quality} " \
					   "{ont_read_dir}{sampleName}.{ont_read_suffix} > {ont_clean_read_folder}{sampleName}.fastq "\
		.format(ont_read_suffix=GlobalParameter().ont_read_suffix,
				sampleName=self.sampleName,
				ont_read_dir=GlobalParameter().ont_read_dir,
				ont_clean_read_folder=ont_clean_read_folder,
				long_read_keep_percent=self.long_read_keep_percent,
				long_read_mean_quality=self.long_read_mean_quality,
				long_read_min_length=self.long_read_min_length)




		cmd_clean_pac = "[ -d  {pac_clean_read_folder} ] || mkdir -p {pac_clean_read_folder}; " \
					   "cd {pac_clean_read_folder}; " \
					   "filtlong --min_length {long_read_min_length} " \
					   "--keep_percent {long_read_keep_percent} " \
					   "--min_mean_q {long_read_mean_quality} " \
					   "{pac_read_dir}{sampleName}.{pac_read_suffix} > {pac_clean_read_folder}{sampleName}.fastq "\
		.format(pac_read_suffix=GlobalParameter().pac_read_suffix,
				sampleName=self.sampleName,
				pac_read_dir=GlobalParameter().pac_read_dir,
				pac_clean_read_folder=pac_clean_read_folder,
				long_read_keep_percent=self.long_read_keep_percent,
				long_read_mean_quality=self.long_read_mean_quality,
				long_read_min_length=self.long_read_min_length)

		if self.seq_platforms=="ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_ont)
			print(run_cmd(cmd_clean_ont))

		if self.seq_platforms=="pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_clean_pac)
			print(run_cmd(cmd_clean_pac))

'''
class cleanLongReads(luigi.Task):
	platform = luigi.ChoiceParameter(description="Choose From[ont: nanopore, pac: pacbio]",
                                     choices=["ont","pac"], var_type=str)

	def requires(self):
		if self.platform=="ont":
			return [
				[filtlong(sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]],

				[readqc(seq_platforms="ont", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]
				]

		if self.platform=="pac":
			return [
				[filtlong(sampleName=i)
				for i in [line.strip()
						  for line in
						  open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]],

				[readqc(seq_platforms="pac", sampleName=i)
				 for i in [line.strip()
						   for line in
						   open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]
				]


	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.clean.longread.complete.{t}'.format(t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('long read processing finished at {t}'.format(t=timestamp))
'''