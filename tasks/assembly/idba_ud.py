#!/usr/bin/env python3
import os
import luigi
import subprocess
from tasks.readCleaning.preProcessReads import cleanFastq
from tasks.readCleaning.preProcessReads import filtlong
from tasks.readCleaning.reFormatReads import reformat
from tasks.readCleaning.reFormatReads import reformatReads
from tasks.readCleaning.necat_correct import correctONT
from tasks.readCleaning.mecat2_correct import correctPAC

class GlobalParameter(luigi.Config):
    threads = luigi.Parameter()
    maxMemory = luigi.Parameter()
    projectName = luigi.Parameter()
    domain=luigi.Parameter()
    assembly_name=luigi.Parameter()
    #seq_platforms=luigi.Parameter()
    pe_read_dir=luigi.Parameter()


def run_cmd(cmd):
    p = subprocess.Popen(cmd, bufsize=-1,
                         shell=True,
                         universal_newlines=True,
                         stdout=subprocess.PIPE,
                         executable='/bin/bash')
    output = p.communicate()[0]
    return output


class idba(luigi.Task):
    projectName = GlobalParameter().projectName
    domain=GlobalParameter().domain
    assembly_name = GlobalParameter().assembly_name
    pre_process_reads = luigi.ChoiceParameter(choices=["yes","no"],var_type=str)
    seq_platforms = luigi.ChoiceParameter(description="Choose From['pe: paired-end','pe-mp: paired-end and mate-pair',pe-ont: paired-end and nanopore, pe-pac: paired-end and pacbio, ont: nanopore, pac: pacbio]",
                                          choices=["pe", "pe-mp", "pe-ont", "pe-pac"], var_type=str)

    def requires(self):
        #Paired-end
        if self.seq_platforms == "pe" and self.pre_process_reads=="yes":
            return [cleanFastq(seq_platforms="pe",sampleName=i)
                for i in [line.strip()
                          for line in
                          open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

        if self.seq_platforms == "pe" and self.pre_process_reads=="no":
            return [reformat(seq_platforms="pe",sampleName=i)
                for i in [line.strip()
                          for line in
                          open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]]

        # Paired-end with Mate-pair
        if self.seq_platforms == "pe-mp" and self.pre_process_reads =="yes":
            return [
                        [cleanFastq(seq_platforms="pe",sampleName=i)
                            for i in [line.strip()  for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [cleanFastq(seq_platforms="mp",sampleName=i)
                            for i in [line.strip()  for line in
                                open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
                    ]

        if self.seq_platforms == "pe-mp" and self.pre_process_reads =="no":
            return [
                        [reformat(seq_platforms="pe",sampleName=i)
                            for i in [line.strip()  for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [reformat(seq_platforms="mp",sampleName=i)
                            for i in [line.strip()  for line in
                                open((os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")))]]
                    ]
        #Paired-end with Pacbio
        if self.seq_platforms == "pe-pac" and self.pre_process_reads == "yes":
            return [
                        [cleanFastq(seq_platforms="pe", sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [filtlong(seq_platforms="pac",sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]
                    ]


        if self.seq_platforms == "pe-pac" and self.pre_process_reads == "no":
            return [
                        [reformat(seq_platforms="pe", sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [reformat(seq_platforms="pac",sampleName=i)
                            for i in [line.strip()for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")))]]
                    ]


        #######
        #Paired-end with Nanopore
        if self.seq_platforms == "pe-ont" and self.pre_process_reads == "yes":
            return [
                        [cleanFastq(seq_platforms="pe", sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [filtlong(seq_platforms="pac",sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]
                    ]


        if self.seq_platforms == "pe-ont" and self.pre_process_reads == "no":
            return [
                        [reformat(seq_platforms="pe", sampleName=i)
                            for i in [line.strip() for line in
                                open((os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")))]],

                        [reformat(seq_platforms="ont",sampleName=i)
                            for i in [line.strip()for line in
                                open((os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")))]]
                    ]

    def output(self):
        idba_assembly_folder = os.path.join(os.getcwd(), self.projectName,"GenomeAssembly", "IDBA_" + self.seq_platforms, self.assembly_name + "/")
        return {'out': luigi.LocalTarget(idba_assembly_folder + "scaffold.fa")}

    def run(self):

        idba_assembly_folder = os.path.join(os.getcwd(), self.projectName,"GenomeAssembly", "IDBA_" + self.seq_platforms, self.assembly_name + "/")
        idba_assembly_log_folder = os.path.join(os.getcwd(), self.projectName,"log", "GenomeAssembly", "IDBA",self.assembly_name + "/")

        if self.seq_platforms=="pe":
            pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")

        if self.seq_platforms=="pe-mp":
            pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")
            mp_sample_list = os.path.join(os.getcwd(), "sample_list", "mp_samples.lst")

        if self.seq_platforms=="pe-ont":
            pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")
            ont_sample_list = os.path.join(os.getcwd(), "sample_list", "ont_samples.lst")
            
        if self.seq_platforms=="pe-pac":
            pe_sample_list = os.path.join(os.getcwd(), "sample_list", "pe_samples.lst")
            pac_sample_list = os.path.join(os.getcwd(), "sample_list", "pac_samples.lst")



        verified_pe_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "PE-Reads" + "/")
        verified_mp_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "MP-Reads" + "/")
        verified_ont_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "ONT-Reads" + "/")
        verified_pac_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "VerifiedReads", "PAC-Reads" + "/")

        cleaned_pe_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "PE-Reads" + "/")
        cleaned_mp_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "MP-Reads" + "/")
        cleaned_ont_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "ONT-Reads" + "/")
        cleaned_pac_read_folder = os.path.join(os.getcwd(), GlobalParameter().projectName, "ReadQC", "CleanedReads", "PAC-Reads" + "/")


        def idba_illumina(samplefile,inputDir):
            with open(samplefile) as fh:
                sample_name_list = fh.read().splitlines()
                left_read_name_suffix = '_R1.fastq'
                right_read_name_suffix = '_R2.fastq'
                left_read_name_list = [x + left_read_name_suffix for x in sample_name_list]
                right_read_name_list = [x + right_read_name_suffix for x in sample_name_list]
                pe_cleaned_read_folder = inputDir

                left_reads_list=[pe_cleaned_read_folder + x for x in left_read_name_list]
                right_reads_list=[pe_cleaned_read_folder + x for x in right_read_name_list]
                left_reads=' '.join(left_reads_list)
                right_reads=' '.join(right_reads_list)

                return left_reads,right_reads

        def idba_longread(samplefile,inputDir):
            with open(samplefile) as fh:
                sample_name_list = fh.read().splitlines()
                read_name_suffix = '.fastq'
                read_name_list = [x + read_name_suffix for x in sample_name_list]
                input_read_folder = inputDir
                reads_list = [input_read_folder + x for x in read_name_list]
                long_reads = ' '.join(reads_list)
                return long_reads

              
        if self.pre_process_reads=="yes" and self.seq_platforms=="pe":
            pe_left_reads,pe_right_reads= idba_illumina(pe_sample_list,cleaned_pe_read_folder)

        if self.pre_process_reads=="yes" and self.seq_platforms=="pe-mp":
            pe_left_reads,pe_right_reads= idba_illumina(pe_sample_list,cleaned_pe_read_folder)
            mp_left_reads,mp_right_reads= idba_illumina(pe_sample_list,cleaned_mp_read_folder)

        if self.pre_process_reads=="yes" and self.seq_platforms=="pe-ont":
            pe_left_reads,pe_right_reads = idba_illumina(pe_sample_list,cleaned_pe_read_folder)
            ont_reads = idba_longread(ont_sample_list,cleaned_ont_read_folder)   

        if self.pre_process_reads=="yes" and self.seq_platforms=="pe-pac":
            pe_left_reads,pe_right_reads = idba_illumina(pe_sample_list,cleaned_pe_read_folder)
            pac_reads = idba_longread(pac_sample_list,cleaned_pac_read_folder)

        if self.pre_process_reads=="no" and self.seq_platforms=="pe":
            pe_left_reads,pe_right_reads= idba_illumina(pe_sample_list,verified_pe_read_folder)

        if self.pre_process_reads=="no" and self.seq_platforms=="pe-mp":
            pe_left_reads,pe_right_reads= idba_illumina(pe_sample_list,verified_pe_read_folder)
            mp_left_reads,mp_right_reads= idba_illumina(pe_sample_list,verified_mp_read_folder)

        if self.pre_process_reads=="no" and self.seq_platforms=="pe-ont":
            pe_left_reads,pe_right_reads = idba_illumina(pe_sample_list,verified_pe_read_folder)
            ont_reads = idba_longread(ont_sample_list,verified_ont_read_folder)   

        if self.pre_process_reads=="no" and self.seq_platforms=="pe-pac":
            pe_left_reads,pe_right_reads = idba_illumina(pe_sample_list,verified_pe_read_folder)
            pac_reads = idba_longread(pac_sample_list,verified_pac_read_folder)
            
       

        if self.seq_platforms=="pe":
            run_merge_pe_fastq= "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                            "cat {pe_left_reads} > {assembly_name}_pe_R1.fastq; cat {pe_right_reads} > {assembly_name}_pe_R2.fastq;" \
                            .format(idba_assembly_folder=idba_assembly_folder,pe_left_reads=pe_left_reads,pe_right_reads=pe_right_reads, assembly_name=self.assembly_name)


            run_cmd_idba_pe_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa --merge {idba_assembly_folder}{assembly_name}_pe_R1.fastq {idba_assembly_folder}{assembly_name}_pe_R2.fastq {assembly_name}_pe.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)



            run_cmd_idba_pe_ud = "cd {idba_assembly_folder} ; " \
                           "mkdir -p {idba_assembly_log_folder}; " \
                           "/usr/bin/time -v idba_ud " \
                           "-r {assembly_name}_pe.fas  " \
                           "--num_threads {threads} " \
                           "-o {idba_assembly_folder} " \
                           "2>&1 | tee {idba_assembly_log_folder}idba_assembly.log " \
            .format(idba_assembly_folder=idba_assembly_folder,
                    idba_assembly_log_folder=idba_assembly_log_folder,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name)

            cmd_rename_contigs = "cp {idba_assembly_folder}contig.fa " \
                             "{idba_assembly_folder}{assembly_name}_contigs.fasta ".format(idba_assembly_folder=idba_assembly_folder,
                              assembly_name=self.assembly_name)


            print("****** NOW RUNNING COMMAND ******: " + run_merge_pe_fastq)
            run_cmd(run_merge_pe_fastq)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_merge)
            run_cmd(run_cmd_idba_pe_merge)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_ud)
            run_cmd(run_cmd_idba_pe_ud)
            run_cmd(cmd_rename_contigs)

        if self.seq_platforms=="pe-mp":
            run_merge_pe_fastq= "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                            "cat {pe_left_reads} > {assembly_name}_pe_R1.fastq; cat {pe_right_reads} > {assembly_name}_pe_R2.fastq;" \
                            .format(idba_assembly_folder=idba_assembly_folder,pe_left_reads=pe_left_reads,pe_right_reads=pe_right_reads, assembly_name=self.assembly_name)

            run_merge_mp_fastq= "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                "cat {mp_left_reads} > {assembly_name}_mp_R1.fastq; cat {mp_right_reads} > {assembly_name}_mp_R2.fastq;" \
                                .format(idba_assembly_folder=idba_assembly_folder,mp_left_reads=mp_left_reads,mp_right_reads=mp_right_reads, assembly_name=self.assembly_name)

            run_cmd_idba_pe_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa --merge {idba_assembly_folder}{assembly_name}_pe_R1.fastq {idba_assembly_folder}{assembly_name}_pe_R2.fastq {assembly_name}_pe.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)

            run_cmd_idba_mp_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa --merge {idba_assembly_folder}{assembly_name}_mp_R1.fastq {idba_assembly_folder}{assembly_name}_mp_R2.fastq {assembly_name}_mp.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)

            run_cmd_idba_pe_mp_hybrid = "cd {idba_assembly_folder} ; " \
                           "mkdir -p {idba_assembly_log_folder}; " \
                           "/usr/bin/time -v idba_hybrid " \
                           "-r {assembly_name}_pe.fas  " \
                           "--read_level_2 {assembly_name}_mp.fas "\
                           "--num_threads {threads} " \
                           "-o {idba_assembly_folder} " \
                           "2>&1 | tee {idba_assembly_log_folder}idba_assembly.log " \
            .format(idba_assembly_folder=idba_assembly_folder,
                    idba_assembly_log_folder=idba_assembly_log_folder,
                    threads=GlobalParameter().threads,
                    assembly_name=self.assembly_name)


            cmd_rename_contigs = "cp {idba_assembly_folder}contig.fa " \
                             "{idba_assembly_folder}{assembly_name}_contigs.fasta ".format(idba_assembly_folder=idba_assembly_folder,
                              assembly_name=self.assembly_name)

            print("****** NOW RUNNING COMMAND ******: " + run_merge_pe_fastq)
            run_cmd(run_merge_pe_fastq)
            print("****** NOW RUNNING COMMAND ******: " + run_merge_mp_fastq)
            run_cmd(run_merge_mp_fastq)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_merge)
            run_cmd(run_cmd_idba_pe_merge)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_mp_merge)
            run_cmd(run_cmd_idba_mp_merge)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_mp_hybrid)
            run_cmd(run_cmd_idba_pe_mp_hybrid)
            run_cmd(cmd_rename_contigs)



        if self.seq_platforms=="pe-ont":
            run_merge_pe_fastq= "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                            "cat {pe_left_reads} > {assembly_name}_pe_R1.fastq; cat {pe_right_reads} > {assembly_name}_pe_R2.fastq;" \
                            .format(idba_assembly_folder=idba_assembly_folder,pe_left_reads=pe_left_reads,pe_right_reads=pe_right_reads, assembly_name=self.assembly_name)

            run_cmd_idba_pe_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa --merge {idba_assembly_folder}{assembly_name}_pe_R1.fastq {idba_assembly_folder}{assembly_name}_pe_R2.fastq {assembly_name}_pe.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)

            run_merge_ont_fastq="[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                             "cat {ont_reads} > {assembly_name}_ont.fastq; ".format(idba_assembly_folder=idba_assembly_folder,ont_reads=ont_reads,assembly_name=self.assembly_name)

            run_cmd_idba_ont_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa {idba_assembly_folder}{assembly_name}_ont.fastq {assembly_name}_ont.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)

            run_cmd_idba_hybrid_pe_ont= "cd {idba_assembly_folder} ; " \
                                "mkdir -p {idba_assembly_log_folder}; " \
                                "/usr/bin/time -v idba_hybrid " \
                                "-r {assembly_name}_pe.fas  " \
                                "-l {assembly_name}_ont.fas "\
                                "--num_threads {threads} " \
                                "-o {idba_assembly_folder} " \
                                "2>&1 | tee {idba_assembly_log_folder}idba_assembly.log " \
                                .format(idba_assembly_folder=idba_assembly_folder,
                                idba_assembly_log_folder=idba_assembly_log_folder,
                                threads=GlobalParameter().threads,
                                assembly_name=self.assembly_name)


            cmd_rename_contigs = "cp {idba_assembly_folder}contig.fa " \
                             "{idba_assembly_folder}{assembly_name}_contigs.fasta ".format(idba_assembly_folder=idba_assembly_folder,
                              assembly_name=self.assembly_name)


            print("****** NOW RUNNING COMMAND ******: " + run_merge_pe_fastq)
            run_cmd(run_merge_pe_fastq)
            print("****** NOW RUNNING COMMAND ******: " + run_merge_ont_fastq)
            run_cmd(run_merge_ont_fastq)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_merge)
            run_cmd(run_cmd_idba_pe_merge)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_ont_merge)
            run_cmd(run_cmd_idba_ont_merge)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_hybrid_pe_ont)
            run_cmd(run_cmd_idba_hybrid_pe_ont)
            run_cmd(cmd_rename_contigs)


        if self.seq_platforms=="pe-pac":
            run_merge_pe_fastq= "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                            "cat {pe_left_reads} > {assembly_name}_pe_R1.fastq; cat {pe_right_reads} > {assembly_name}_pe_R2.fastq;" \
                            .format(idba_assembly_folder=idba_assembly_folder,pe_left_reads=pe_left_reads,pe_right_reads=pe_right_reads, assembly_name=self.assembly_name)

            run_cmd_idba_pe_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa --merge {idba_assembly_folder}{assembly_name}_pe_R1.fastq {idba_assembly_folder}{assembly_name}_pe_R2.fastq {assembly_name}_pe.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)

            run_merge_pac_fastq="[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                             "cat {pac_reads} > {assembly_name}_pac.fastq; ".format(idba_assembly_folder=idba_assembly_folder,pac_reads=pac_reads,assembly_name=self.assembly_name)

            run_cmd_idba_pac_merge = "[ -d  {idba_assembly_folder} ] || mkdir -p {idba_assembly_folder}; cd {idba_assembly_folder};" \
                                " fq2fa {idba_assembly_folder}{assembly_name}_pac.fastq {assembly_name}_pac.fas ". \
                                format(idba_assembly_folder=idba_assembly_folder,assembly_name=self.assembly_name)


            run_cmd_idba_hybrid_pe_pac= "cd {idba_assembly_folder} ; " \
                                "mkdir -p {idba_assembly_log_folder}; " \
                                "/usr/bin/time -v idba_hybrid " \
                                "-r {assembly_name}_pe.fas  " \
                                "-l {assembly_name}_pac.fas "\
                                "--num_threads {threads} " \
                                "-o {idba_assembly_folder} " \
                                "2>&1 | tee {idba_assembly_log_folder}idba_assembly.log " \
                                .format(idba_assembly_folder=idba_assembly_folder,
                                idba_assembly_log_folder=idba_assembly_log_folder,
                                threads=GlobalParameter().threads,
                                assembly_name=self.assembly_name)

            cmd_rename_contigs = "cp {idba_assembly_folder}contig.fa " \
                              "{idba_assembly_folder}{assembly_name}_contigs.fasta ".format(idba_assembly_folder=idba_assembly_folder,
                              assembly_name=self.assembly_name)

            print("****** NOW RUNNING COMMAND ******: " + run_merge_pe_fastq)
            run_cmd(run_merge_pe_fastq)
            print("****** NOW RUNNING COMMAND ******: " + run_merge_pac_fastq)
            run_cmd(run_merge_pac_fastq)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pe_merge)
            run_cmd(run_cmd_idba_pe_merge)
            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_pac_merge)
            run_cmd(run_cmd_idba_pac_merge)

            print("****** NOW RUNNING COMMAND ******: " + run_cmd_idba_hybrid_pe_pac)
            run_cmd(run_cmd_idba_hybrid_pe_pac)
            run_cmd(cmd_rename_contigs)