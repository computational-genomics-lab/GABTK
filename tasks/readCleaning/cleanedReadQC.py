import luigi
import time
import os
import subprocess

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

	#seq_platforms=luigi.Parameter()

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

class postreadqc(luigi.Task):
	#projectName = luigi.Parameter(default="ReadQC")
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)


	def output(self):
		pe_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PostQC", "PE-Reads" + "/")
		mp_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","MP-Reads" + "/")
		ont_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","ONT-Reads" + "/")
		pac_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","PACBIO-Reads" + "/")

		
		if self.seq_platforms == "pe":
			return {'out1': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R2_fastqc.html")}

		if self.seq_platforms == "ont":
			return {'out1': luigi.LocalTarget(ont_readQC_folder + self.sampleName + "_nanoQC.html")}


		if self.seq_platforms == "mp":
			return {'out1': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R1_fastqc.html")}


		if self.seq_platforms == "pe-mp":
			return {'out1': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R2_fastqc.html"),
					'out3': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out4': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R2_fastqc.html")}

		if self.seq_platforms == "pe-ont":
			return {'out1': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R2_fastqc.html"),
					'out3': luigi.LocalTarget(ont_readQC_folder + self.sampleName + "_nanoQC.html")}

		if self.seq_platforms == "pe-pac":
			return {'out1': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R2_fastqc.html"),
					'out3': luigi.LocalTarget(pac_readQC_folder + self.sampleName + "_nanoQC.html")}

		

	def run(self):
		pe_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PostQC", "PE-Reads" + "/")
		mp_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","MP-Reads" + "/")
		ont_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","ONT-Reads" + "/")
		pac_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC", "PostQC","PACBIO-Reads" + "/")

				
		pe_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","PE-Reads" + "/")
		mp_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","MP-Reads" + "/")
		ont_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","ONT-Reads" + "/")
		pac_clean_read_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","CleanedReads","PACBIO-Reads" + "/")


		read_QC_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "CleanedReads" + "/")


		cmd_cleaned_pe_qc = "[ -d  {pe_readQC_folder} ] || mkdir -p {pe_readQC_folder}; mkdir -p {read_QC_log_folder}; " \
					   "/usr/bin/time -v fastqc " \
						"-t {cpu} " \
						"{pe_clean_read_folder}{sampleName}_R1.fastq " \
						"{pe_clean_read_folder}{sampleName}_R2.fastq " \
						"-o {pe_readQC_folder} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_cleaned_pe_fastqc.log".format(

													   sampleName=self.sampleName,
													   pe_readQC_folder=pe_readQC_folder,
													   cpu=GlobalParameter().threads,
													   pe_clean_read_folder=pe_clean_read_folder,
													   read_QC_log_folder=read_QC_log_folder)

		cmd_cleaned_mp_qc = "[ -d  {mp_readQC_folder} ] || mkdir -p {mp_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"fastqc " \
						"-t {cpu} " \
						"{mp_clean_read_folder}{sampleName}_R1.fastq " \
						"{mp_clean_read_folder}{sampleName}_R2.fastq " \
						"-o {mp_readQC_folder} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_cleaned_mp_fastqc.log".format(
													   sampleName=self.sampleName,
													   mp_readQC_folder=mp_readQC_folder,
													   cpu=GlobalParameter().threads,
													   read_QC_log_folder=read_QC_log_folder,
													   mp_clean_read_folder=mp_clean_read_folder)

		cmd_cleaned_ont_qc = "[ -d  {ont_readQC_folder} ] || mkdir -p {ont_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"nanoQC -o {ont_readQC_folder} " \
						"{ont_clean_read_folder}{sampleName}.fastq " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_cleaned_lr_nanoqc.log".format(sampleName=self.sampleName,
													   ont_readQC_folder=ont_readQC_folder,
													   read_QC_log_folder=read_QC_log_folder,
													   ont_clean_read_folder=ont_clean_read_folder)

		cmd_cleaned_pac_qc = "[ -d  {pac_readQC_folder} ] || mkdir -p {pac_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"nanoQC -o {pac_readQC_folder} " \
						"{pac_clean_read_folder}{sampleName}.fastq " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_cleaned_lr_nanoqc.log".format(sampleName=self.sampleName,
													   pac_readQC_folder=pac_readQC_folder,
													   read_QC_log_folder=read_QC_log_folder,
													   pac_clean_read_folder=pac_clean_read_folder)




		cmd_mv_ont_qc = "cd {ont_readQC_folder};  " \
						"mv nanoQC.html {sampleName}_nanoQC.html ".format(sampleName=self.sampleName,
													   ont_readQC_folder=ont_readQC_folder)

		cmd_mv_pac_qc = "cd {pac_readQC_folder};  " \
						"mv nanoQC.html {sampleName}_nanoQC.html ".format(sampleName=self.sampleName,
													   pac_readQC_folder=pac_readQC_folder)

		if self.seq_platforms == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pe_qc)
			print (run_cmd(cmd_cleaned_pe_qc))

		if self.seq_platforms == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_mp_qc)
			print (run_cmd(cmd_cleaned_mp_qc))

		if self.seq_platforms == "ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_ont_qc)
			print (run_cmd(cmd_cleaned_ont_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_ont_qc)
			print(run_cmd(cmd_mv_ont_qc))

		if self.seq_platforms == "pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pac_qc)
			print (run_cmd(cmd_cleaned_pac_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_pac_qc)
			print(run_cmd(cmd_mv_pac_qc))

		if self.seq_platforms == "pe-mp" :
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pe_qc)
			print(run_cmd(cmd_cleaned_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_mp_qc)
			print(run_cmd(cmd_cleaned_mp_qc))
		

		if self.seq_platforms == "pe-ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pe_qc)
			print(run_cmd(cmd_cleaned_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_ont_qc)
			print(run_cmd(cmd_cleaned_ont_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_ont_qc)
			print(run_cmd(cmd_mv_ont_qc))


		if self.seq_platforms == "pe-pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pe_qc)
			print(run_cmd(cmd_cleaned_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_cleaned_pac_qc)
			print(run_cmd(cmd_cleaned_pac_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_pac_qc)
			print(run_cmd(cmd_mv_pac_qc))


class cleanReadsQC(luigi.Task):
	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)

	def requires(self):

		if self.seq_platforms == "pe":
			return [postreadqc(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

		if self.seq_platforms == "ont":
			return [postreadqc(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]

		if self.seq_platforms == "pac":
			return [postreadqc(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]


		if self.seq_platforms == "mp":
			return [postreadqc(seq_platforms=self.seq_platforms,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(),"sample_list", "mp_samples.lst")))]]


		if self.seq_platforms == "pe-mp":

			return [
						[postreadqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[postreadqc(seq_platforms="mp", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","mp_samples.lst")))]]
				  ]


		if self.seq_platforms == "pe-ont":

			return [
						[postreadqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[postreadqc(seq_platforms="lr", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","lr_samples.lst")))]]
				  ]


		if self.seq_platforms == "pe-pac":

			return [
						[postreadqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[postreadqc(seq_platforms="lr", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pac_samples.lst")))]]
				  ]

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.clean.read.qc.analysis.complete.{t}'.format(
			t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Cleaned Read QC Assessment finished at {t}'.format(t=timestamp))










		





