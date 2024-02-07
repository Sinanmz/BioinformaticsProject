import subprocess
import pandas as pd
import os

def DESeq2_study_group(input_file_path, output_file_path, social_settings=[]):
    r_script_path = '/Users/sina/Sem-4/bioinformatics/project/task_3/tasks/the-songbird-study/rscript.R'
    if len(social_settings) != 0:
        # filter the input file based on the social settings
        df = pd.read_csv(input_file_path)
        df = df[df['social_settting'].isin(social_settings)]
        df.to_csv('temp', index=False)
        input_file_path = 'temp'
    
    command = ['Rscript', r_script_path, '0', input_file_path, output_file_path]
    result = subprocess.run(command, capture_output=True, text=True)
    # delete the temp file
    if len(social_settings) != 0:
        os.remove('temp')
    summary_file = output_file_path + ".summary"
    
    # Check if the script ran successfully
    if result.returncode == 0:
        print("DESeq2 analysis completed successfully.")
        with open(summary_file, "r") as f:
            summary_output = f.read()
            print(summary_output)
        
    else:
        print("DESeq2 analysis failed. Error message:", result.stderr)
    
    if os.path.exists(summary_file):
        os.remove(summary_file)



def DESeq2_social_setting(input_file_path, output_file_path, base_setting, second_setting, study_groups='both'):
    r_script_path = '/Users/sina/Sem-4/bioinformatics/project/task_3/tasks/the-songbird-study/rscript.R'

    if study_groups != 'both':
        # filter the input file based on the study groups
        df = pd.read_csv(input_file_path)
        df = df[df['study_group'].isin([study_groups])]
        df.to_csv('temp', index=False)
        input_file_path = 'temp'

    command = ['Rscript', r_script_path, '1', input_file_path, output_file_path, base_setting, second_setting]

    result = subprocess.run(command, capture_output=True, text=True)
    # delete the temp file
    if study_groups != 'both':
        os.remove('temp')
    summary_file = output_file_path + ".summary"
    if result.returncode == 0:
        print("DESeq2 analysis completed successfully.")
        with open(summary_file, "r") as f:
            summary_output = f.read()
            print(summary_output)
    else:
        print("DESeq2 analysis failed. Error message:", result.stderr)

    if os.path.exists(summary_file):
        os.remove(summary_file)
