#!/usr/bin/python
# -*- coding: utf-8 -*-
# CONDA ENV SOFTWARE: python (3.7), sratoolkit (v2.10.3)

import os, sys, subprocess, argparse, time, shutil
from core import CAPTURE, SBATCH, REVIEW, ezSub

pipe_script = os.path.realpath(__file__)
pipe_dir, pipe_name = os.path.split(pipe_script)
pipe_flag = '-pipe'
buddy_script = f'{pipe_dir}/sra_buddy.py'

parser = argparse.ArgumentParser(description='SRAToolkit Pipe', prog=pipe_name, usage='%(prog)s [options]', epilog='see readme file for further details')
basic = parser.add_mutually_exclusive_group(required=True) # basic modes
basic.add_argument('-setup', action='store_const', const='SETUP', help='setup initial directories')
basic.add_argument('-stage', metavar='<number>', type=int, choices=range(1,3), help='specify processing stage')

checks = parser.add_argument_group(description='check progress:') # checking modes
checks.add_argument('-progress', action='store_true', help='review task status')
checks.add_argument('-accounting', action='store_true', help='review task accounting data')

additional = parser.add_argument_group(description='additional modes:') # additional modes
additional.add_argument(pipe_flag, action='store_true', help='submit pipe as task')
additional.add_argument('-overwrite', action='store_true', help='overwrite existing files')
additional.add_argument('-buddy', action='store_true', help='assistance from sra buddy')

inputs = parser.add_argument_group(description='user inputs:') # user inputs
inputs.add_argument('-u', metavar='<name>', type=str, help='specify user name (default: find current user)')
inputs.add_argument('-p', metavar='</path/to/>', type=str, help='specify path (default: path to home directory)')
inputs.add_argument('-d', metavar='<name>', type=str, default='sra_download', help='specify working directory (default: sra_download)')
inputs.add_argument('-a', metavar='<name>', type=str, default='SraAccList.txt', help='specify accessions list file (default: SraAccList.txt)')
inputs.add_argument('-e', metavar='<name>', type=str, default='sra_env', help='specify conda environment (default: sra_env)')
inputs.add_argument('-l', metavar='<number>', type=int, default=100, help='specify parallel task limit (default: 500)')

hpcc = parser.add_argument_group(description='hpcc settings:') # hpcc settings
hpcc.add_argument('-pa', metavar='<name>', type=str, default='defq', help='specify partition (default: defq)')
hpcc.add_argument('-no', metavar='<number>', type=int, default=1, help='specify nodes (default: 1)')
hpcc.add_argument('-nt', metavar='<number>', type=int, default=1, help='specify ntasks (default: 1)')
hpcc.add_argument('-me', metavar='<number[units]>', type=str.lower, default='4g', help='specify memory [k|m|g|t] (default: stage-specific)')
hpcc.add_argument('-wt', metavar='<HH:MM:SS>', type=str, default='24:00:00', help='specify walltime (default: stage-specific)')

setup, stage, checking_progress, reviewing_accounting, submitting_self, overwriting, submitting_buddy, user, path, directory, accessions_list, environment, limit, *hpcc_settings = vars(parser.parse_args()).values() # define user inputs
partition, nodes, ntasks, memory, walltime = hpcc_settings
system_script, *arguments = sys.argv

error_no_files = '\n*** ERROR: Problem finding {}. ***\n'

if memory and (sum(memory.count(unit) for unit in ['k','m','g','t']) != 1 or not memory.rstrip('kmgt').isdigit()): parser.error(f"argument -me: invalid str format: '{memory}'") # check memory format
if walltime and ( not 8 <= len(walltime) <= 9 or walltime[-6:-2:3].count(':') != 2 or not walltime.replace(':','').isdigit()): parser.error(f"argument -wt: invalid str format: '{walltime}'") # check walltime format
if submitting_buddy:
    if not os.path.exists(buddy_script): print(error_no_files.format(f'\"{buddy_script}\"')) # check sra buddy script exists
    if stage != 1: print('\nSRA Buddy only for use with stage 1.\n') # check sra buddy only being used when downloading

user, path = [ os.getenv(bash) if not argument else argument for argument,bash in [ (user,'USER'),(path,'HOME')] ] # find user, path & directory
path = f'/{path.strip("/")}/{directory.strip("/")}' # ensure correct path format

output_labels = ['accessions','downloads','tmp', 'fastqs'] # specify output directory names
batch_labels = [ f'x.slurm/{label}' for label in ['scripts','out.err','ids'] ] # specify batch directory names
accession_dir, download_dir, tmp_dir, fastq_dir, sh_dir, oe_dir, id_dir = base_dirs = [ f'{path}/{name}' for name in [*output_labels, *batch_labels] ] # specify output & batch directory paths
id_file, accessions_file = [ f'{root}/{name}' for root,name in [(id_dir,f'{stage}-task.ids'), (accession_dir, accessions_list)] ] # specify id & accession files
download_suffix, downloading_suffix, converting_suffix, compressed_suffix = '.sra', ('.tmp', '.lock'), '.fastq', '.fastq.gz' # specify output file extensions

if setup:
    os.makedirs(accession_dir, exist_ok=True) # create initial directories
    print('\nSETUP COMPLETE\n')

else:
    if not os.path.isdir(path): print(error_no_files.format(f'\"{path}\"')) # check working directory exists

    else: # proceed with reviewing / downloading / converting

        if reviewing_accounting: # proceed with reviewing task accounting data
            if not os.path.exists(id_file): print(error_no_files.format(f'\"{id_file}\"')) # check id file exists
            else: REVIEW(id_file) # review task accounting data via job ids

        else: # proceed with finding tasks to process

            [ os.makedirs(path, exist_ok=True) for path in base_dirs ] # create remaining output & batch directories

            accessions = [ accession.strip() for accession in open(accessions_file,'r').readlines() if accession.strip() ] # read accessions from accessions file

            searches = [ # specify file extensions to search & the directories they are located
                (download_dir, downloading_suffix),
                (download_dir, download_suffix),
                (fastq_dir, converting_suffix),
                (fastq_dir, compressed_suffix)]

            started_downloading, started_converting = [ os.listdir(directory) for directory in [download_dir, fastq_dir] ] # list output subdirectories created by sra toolkit
            partially_downloaded, downloaded, partially_converted, converted = [ { os.path.basename(root): [root, *files] for root,dirs,files in os.walk(directory, topdown=False) if any(name.endswith(suffixes) for name in files) } for directory,suffixes in searches ] # list output files found         

            stage_formats = { # specify stage specific input / output directories / files
                1: (accessions, download_dir, started_downloading, partially_downloaded, downloaded),
                2: (started_downloading, fastq_dir, started_converting, partially_converted, converted)}
            stage_inputs, stage_output_dir, stage_output_subdirs, stage_partial_outputs, stage_outputs = stage_formats[stage] # specify stage specific directories & files
            
            subsets = [
                set(stage_outputs).difference(set(stage_partial_outputs)), # partial outputs removed from complete outputs (both exist simultaneously during fastq compression)
                set(stage_outputs).union(set(stage_partial_outputs))] # partial outputs combied with complete outputs
            stage_completed_outputs, stage_any_outputs = [ set(stage_inputs).intersection( subset )  for subset in subsets ] # filter subsets for only intended inputs

            comparisons = [ # specify files lists to compare
                (stage_inputs, stage_output_subdirs), # tasks to process vs tasks started
                (stage_output_subdirs, stage_completed_outputs), # tasks started vs tasks completed
                (started_downloading, downloaded)] # accession subdirs to convert vs available accession files
            pending, progressing, inputs_missing = [ set(stage_inputs).intersection( set(superset).difference(set(subset)) ) for superset, subset in comparisons ] # intersect files lists to determine task status

            if checking_progress: # proceed with reviewing task progress
                print('\nCHECKING TASKS...')
                if len(stage_completed_outputs) == len(stage_inputs): # check tasks are currently in progress
                    print(f'\nAll {len(stage_completed_outputs)} tasks completed.')
                else:
                    print(f'\n{len(stage_inputs)} Tasks found.')
                    summary = { # task statuses
                        'Pending': pending,
                        'In Progress': progressing,
                        'Completed': stage_completed_outputs}
                    for category, samples in summary.items():
                        if samples:
                            print(f'\n{len(samples)} {category}:\n')
                            for name in sorted(samples): print(name)
                print('\nCHECK COMPLETE\n')
                
            else: # proceed with processing tasks

                if stage == 2 and inputs_missing: # check download subdirectories contain downloads
                    print(error_no_files.format('all input files'))
                    print('Downloads were missing for the following accessions:')
                    for name in sorted(inputs_missing): print(name) # accession names with absent downloads

                elif len(stage_any_outputs) == len(stage_inputs) and not overwriting: # check all tasks not already in progress / completed
                    print('\nAll accessions either completed or being processed.\n')

                else: # proceed with processing tasks

                    resubmitting = len(stage_output_subdirs) > 0 and not overwriting

                    description = f'stage{stage}' # specify task description
                    
                    scripts = []

                    for accession in sorted(stage_inputs): # cycle through accessions

                        accession_subdir = f'{stage_output_dir}/{accession}'

                        if overwriting: # check output files need to be overwritten
                            shutil.rmtree(accession_subdir, ignore_errors=True) # remove output subdirectory

                        sh_script = f'{sh_dir}/{accession}-{description}.sh'

                        if accession in stage_output_subdirs and not overwriting: # check task not already in progress / completed
                            print(f'skipping {accession} (already {"completed" if accession in stage_completed_outputs else "in progress"})')

                        else: # proceed with processing tasks

                            scripts.append(sh_script)
                            with open(sh_script,'w') as sh:
                                
                                hpcc_directives = SBATCH(job_id=accession, partition=partition, nodes=nodes, ntasks=ntasks, memory=memory, walltime=walltime, out_err=oe_dir, task=description, conda_env=environment) # specify hpcc directives (slurm)
                                sh.write(hpcc_directives)
                                
                                if stage == 1: 
                                    sh.write(f'prefetch {accession} -O {accession_subdir}\n') # download sra spot files
                                
                                if stage == 2: 
                                    download = os.path.join(*downloaded[accession])
                                    sh.write(f'fasterq-dump {download} -O {accession_subdir} -o {accession} -t {tmp_dir}\n'+ # convert sra spot files to fastqs; default settings equivalent to --split-3 (-3) & --skip-technical 
                                    f'find {accession_subdir} -type f ! -name "*.*" -exec mv {{}} {{}}.fastq \;\n'+ # add fastq extension to any exntensionless files (single read files converted via -3)
                                    f'gzip {accession_subdir}/*.fastq\n') # compress fastqs
                
                    if not scripts: 
                        print('\nNo tasks to process.\n')

                    else: # proceed with submitting tasks

                        if len(scripts) > limit and not submitting_self: # check tasks being submitted not excessive
                            print(f'\nThe number of tasks ({len(scripts)}) exceeds the parallel task limit ({limit});  either increase the limit or re-submit the pipeline as a batch job.')

                        if submitting_self: # proceed with submiting self so tasks can be submitted from hpcc

                            sh_script = f'{sh_dir}/pipe.sh'
                            with open(sh_script, 'w') as sh:  

                                hpcc_directives = SBATCH(job_id='pipe', partition=partition, nodes=1, ntasks=1, memory='2g', walltime='168:00:00', out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
                                sh.write(hpcc_directives)                
                                
                                arguments.remove(pipe_flag) # remove self submission flag
                                print('python', pipe_script, *arguments, sep=' ', file=sh) # submit self as task
                            pipe_id = CAPTURE(f'sbatch {sh_script}')
                            print(f'PIPELINE SUBMITTED ({pipe_id})')

                        else:
                            print('\nSUBMITTING TASKS...')
                            
                            with open(id_file, 'a+' if resubmitting  else 'w') as id_output:
                                for i,script in enumerate(scripts,1): 
                                    ezSub(i=i, user=user, limit=limit) # maintain tasks below parellel task limit
                                    sub_id = CAPTURE(f'sbatch -d singleton {script}') # submit task
                                    print(sub_id, script, sep='\t', file=id_output, flush=True) # record task job id & shell script
                            
                            print('\nALL TASKS SUBMITTED\n')

                            if stage == 1 and submitting_buddy: # proceed with submitting sra buddy
                            
                                buddy_num = int(CAPTURE(f'squeue -u {user} | grep -c "buddy"')) # find any sra buddy scripts already submitted
                            
                                if buddy_num == 0: # # proceed with submitting sra buddy if not alredy submitted
                            
                                    sh_script = f'{sh_dir}/buddy.sh'
                                    with open (sh_script,'w') as sh:

                                        hpcc_directives = SBATCH(job_id='buddy', partition=partition, nodes=1, ntasks=1, memory='2g', walltime='168:00:00', out_err=oe_dir, conda_env=environment) # specify hpcc directives (slurm)
                                        sh.write(hpcc_directives)                
                                        
                                        sh.write(f'echo BUDDY CALLED `date` \n'+
                                        f'python {buddy_script} -i {id_file} -d {download_dir} -a {accessions_file}\n') # specify sra buddy script instructions
                                    
                                    buddy_id = CAPTURE(f'sbatch {sh_script}') # submit sra buddy
                                    print(f'SUB BUDDY SUBMITTED ({buddy_id})')

