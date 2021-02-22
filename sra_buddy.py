import os, sys, subprocess, argparse, time
from core import CAPTURE

parser = argparse.ArgumentParser(description='SRAToolkit Buddy', prog=os.path.basename(__file__), usage='%(prog)s [options]', epilog='to assist with downloading files using the SRAToolkit')

parser.add_argument('-i', metavar='</path/to/id_file>', type=str, help='specify path to slurm id file')
parser.add_argument('-d', metavar='</path/to/download_dir>', type=str, help='specify path to download directory')
parser.add_argument('-a', metavar='</path/to/accessions_list>', type=str, help='specify path to accessions list file')

id_file, download_dir, accessions_file = vars(parser.parse_args()).values() # define user inputs

subprocess.call('echo BUDDY STARTED `date`', shell=True); continue_loop = True # log script start

download_suffix, downloading_suffix, fastq_suffix = '.sra', ('.tmp', '.lock'), '.fastq.gz'
accessions = [ accession.strip() for accession in open(accessions_file,'r').readlines() if accession.strip() ]
#time.sleep(120) # wait 2 mins before checking

id_info = [ (sub_id, script) for sub_id, script, *_ in  [ line.strip().split('\t') for line in open(id_file, 'r').readlines() if line.strip() ] ]

while continue_loop is True:

    to_search = [
        (download_dir, download_suffix),
        (download_dir, downloading_suffix)]

    started_downloading = os.listdir(download_dir) # accession directories created by sra toolkit when downloading
    
    downloaded, partially_downloaded = [ [ os.path.basename(root) for root,dirs,files in os.walk(directory, topdown=False) if any(name.endswith(suffixes) for name in files) ] for directory,suffixes in to_search ]            
    
    to_compare = [
        (accessions, started_downloading), # to process vs started
        (started_downloading, downloaded), # started vs completed
        (started_downloading, set(downloaded).union(set(partially_downloaded)))] # started vs completed & partial

    pending, missing, stalled = [ set(total).difference(set(found)) for total, found in to_compare ] 

    if stalled:
        
        for accession in stalled:
            
            script_file, *_ = [ script for sub_id, script in id_info if accession in os.path.basename(script) ]
            
            with open(f'{id_file}', 'a+') as id_update: # append id file
                resub_id = CAPTURE(f'sbatch {script_file}') # resubmit script
                print(resub_id, script_file, '(RESUBMITTED)', sep='\t', file=id_update) # append resubmitted id        

    if not pending and not missing: break # all tasks completed; exit

    time.sleep(300) # wait 5 minutes to check again

subprocess.call('echo BUDDY COMPLETED `date`', shell=True) # log finish
