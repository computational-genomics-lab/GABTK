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

class readqc(luigi.Task):
	paired_end_read_dir = GlobalParameter().pe_read_dir
	mate_pair_read_dir = GlobalParameter().mp_read_dir
	nanopore_read_dir = GlobalParameter().ont_read_dir
	pacbio_read_dir = GlobalParameter().pac_read_dir

	paired_end_read_suffix = GlobalParameter().pe_read_suffix
	mate_pair_read_suffix = GlobalParameter().mp_read_suffix
	nanopore_read_suffix = GlobalParameter().ont_read_suffix
	pacbio_read_suffix = GlobalParameter().pac_read_suffix

	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)

	threads = GlobalParameter().threads
	maxMemory = GlobalParameter().maxMemory
	projectName = GlobalParameter().projectName
	sampleName = luigi.Parameter(description="name of the sample to be analyzed. (string)")
	

	def output(self):
		pe_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "PE-Reads" + "/")
		mp_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "MP-Reads" + "/")
		ont_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "ONT-Reads" + "/")
		pac_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "PACBIO-Reads" + "/")
		
		if self.seq_platforms == "pe":
			return {'out1': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(pe_readQC_folder + self.sampleName + "_R2_fastqc.html")}

		if self.seq_platforms == "mp":
			return {'out1': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R1_fastqc.html"),
					'out2': luigi.LocalTarget(mp_readQC_folder + self.sampleName + "_R1_fastqc.html")}
							

		if self.seq_platforms == "ont":
			return {'out1': luigi.LocalTarget(ont_readQC_folder + self.sampleName + "_nanoQC.html")}


		if self.seq_platforms == "pac":
			return {'out1': luigi.LocalTarget(pacbio_readQC_folder + self.sampleName + "_nanoQC.html")}


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
		pe_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "PE-Reads" + "/")
		mp_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "MP-Reads" + "/")
		ont_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "ONT-Reads" + "/")
		pac_readQC_folder = os.path.join(os.getcwd(), self.projectName,"ReadQC","PreQC", "PACBIO-Reads" + "/")

		read_QC_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "ReadQC", "PreQC" + "/")



		cmd_raw_pe_qc = "[ -d  {pe_readQC_folder} ] || mkdir -p {pe_readQC_folder}; mkdir -p {read_QC_log_folder}; " \
					   "/usr/bin/time -v fastqc " \
						"-t {cpu} " \
						"{paired_end_read_dir}{sampleName}_R1.{paired_end_read_suffix} " \
						"{paired_end_read_dir}{sampleName}_R2.{paired_end_read_suffix} " \
						"-o {pe_readQC_folder} 2>&1 | tee  {read_QC_log_folder}{sampleName}_pe_fastqc.log".format(
													   sampleName=self.sampleName,
													   paired_end_read_suffix=self.paired_end_read_suffix,
													   pe_readQC_folder=pe_readQC_folder,
													   cpu=GlobalParameter().threads,
													   paired_end_read_dir=self.paired_end_read_dir,
													   read_QC_log_folder=read_QC_log_folder)

		cmd_raw_mp_qc = "[ -d  {mp_readQC_folder} ] || mkdir -p {mp_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"fastqc " \
						"-t {cpu} " \
						"{mate_pair_read_dir}{sampleName}_R1.{mate_pair_read_suffix} " \
						"{mate_pair_read_dir}{sampleName}_R2.{mate_pair_read_suffix} " \
						"-o {mp_readQC_folder} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_mp_fastqc.log".format(
													   sampleName=self.sampleName,
													   mate_pair_read_suffix=self.mate_pair_read_suffix,
													   mp_readQC_folder=mp_readQC_folder,
													   cpu=GlobalParameter().threads,
													   read_QC_log_folder=read_QC_log_folder,
													   mate_pair_read_dir=self.mate_pair_read_dir)

		
		cmd_raw_ont_qc = "[ -d  {ont_readQC_folder} ] || mkdir -p {ont_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"nanoQC -o {ont_readQC_folder} " \
						"{ont_read_dir}{sampleName}.{ont_read_suffix} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_ont_nanoqc.log".format(sampleName=self.sampleName,
													   ont_readQC_folder=ont_readQC_folder,
													   ont_read_suffix=self.ont_read_suffix,
													   ont_read_dir=GlobalParameter().ont_read_dir)		

		cmd_raw_pac_qc = "[ -d  {pac_readQC_folder} ] || mkdir -p {pac_readQC_folder};  mkdir -p {read_QC_log_folder}; " \
						"nanoQC -o {pac_readQC_folder} " \
						"{pac_read_dir}{sampleName}.{pac_read_suffix} " \
						"2>&1 | tee  {read_QC_log_folder}{sampleName}_pac_nanoqc.log".format(sampleName=self.sampleName,
													   pac_readQC_folder=pac_readQC_folder,
													   pac_read_suffix=self.pac_read_suffix,
													   pac_read_dir=GlobalParameter().pac_read_dir)					

		cmd_mv_ont_qc = "cd {ont_readQC_folder};  " \
						"mv nanoQC.html {sampleName}_nanoQC.html ".format(sampleName=self.sampleName,
													   ont_readQC_folder=ont_readQC_folder)

		cmd_mv_pac_qc = "cd {pac_readQC_folder};  " \
						"mv nanoQC.html {sampleName}_nanoQC.html ".format(sampleName=self.sampleName,
													   pac_readQC_folder=pac_readQC_folder)

		

		if self.seq_platforms == "pe":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print (run_cmd(cmd_raw_pe_qc))

		if self.seq_platforms == "mp":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_mp_qc)
			print (run_cmd(cmd_raw_mp_qc))

		if self.seq_platforms == "ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_ont_qc)
			print (run_cmd(cmd_raw_ont_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_ont_qc)
			print(run_cmd(cmd_mv_ont_qc))

		if self.seq_platforms == "pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pac_qc)
			print (run_cmd(cmd_raw_pac_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_pac_qc)
			print(run_cmd(cmd_mv_pac_qc))


		if self.seq_platforms == "pe-mp" :
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print(run_cmd(cmd_raw_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_mp_qc)
			print(run_cmd(cmd_raw_mp_qc))
		

		if self.seq_platforms == "pe-ont":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print(run_cmd(cmd_raw_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_ont_qc)
			print(run_cmd(cmd_raw_ont_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_ont_qc)
			print(run_cmd(cmd_mv_ont_qc))


		if self.seq_platforms == "pe-pac":
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pe_qc)
			print(run_cmd(cmd_raw_pe_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_raw_pac_qc)
			print(run_cmd(cmd_raw_pac_qc))
			print("****** NOW RUNNING COMMAND ******: " + cmd_mv_pac_qc)
			print(run_cmd(cmd_mv_pac_qc))



class rawReadsQC(luigi.Task):
	seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                             choices=["pe", "mp","pe-mp", "pe-ont", "pe-pac","ont","pac"], var_type=str)

	def requires(self):

		if self.seq_platforms == "pe":
			return [readqc(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

		if self.seq_platforms == "ont":
			return [readqc(seq_platforms=self.seq_platforms,
						sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]


		if self.seq_platforms == "mp":
			return [readqc(seq_platforms=self.seq_platforms,
						   sampleName=i)
					for i in [line.strip()
							  for line in
							  open((os.path.join(os.getcwd(),"sample_list", "mp_samples.lst")))]]


		if self.seq_platforms == "pe-mp":

			return [
						[readqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[readqc(seq_platforms="mp", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","mp_samples.lst")))]]
				  ]


		if self.seq_platforms == "pe-ont":

			return [
						[readqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[readqc(seq_platforms="lr", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","ont_samples.lst")))]]
				  ]


		if self.seq_platforms == "pe-pac":

			return [
						[readqc(seq_platforms="pe",sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pe_samples.lst")))]],

						[readqc(seq_platforms="lr", sampleName=i)
								for i in [line.strip()
										  for line in
												open((os.path.join(os.getcwd(), "sample_list","pac_samples.lst")))]]
					]
		

	def output(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		return luigi.LocalTarget(os.path.join(os.getcwd(),"task_logs",'task.read.qc.analysis.complete.{t}'.format(
			t=timestamp)))

	def run(self):
		timestamp = time.strftime('%Y%m%d.%H%M%S', time.localtime())
		with self.output().open('w') as outfile:
			outfile.write('Read QC Assessment finished at {t}'.format(t=timestamp))










		





