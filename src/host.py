import subprocess
import os.path

def run_shell_command(cmd_lst) :
    """Returns a tuple: (stdout, stderr).
    Note: you must separate the options from their values.
    Note: the output is returned when all calculations are done
    Example: you want 'plink --bfile <file> ...'. Must provide cmd_list=['plink', '--bfile', '<file>', ...]"""
    process = subprocess.Popen(cmd_lst, stdout=subprocess.PIPE, 
                               stderr=subprocess.PIPE, universal_newlines=True)
    return process.communicate()

def run_plink(command, file, out, extension, plink2=True, verbose_err=True, force=False):
    """Returns a tuple: (stdout, stderr).
    Note: command input must have the options/values separated by a single space"""
    # Check is result already exist
    stdout, stderr = '', ''
    if extension[0] == '.': # avoid unexpected behaviour
        stderr = 'ERROR: run_plink function, extension must not have a dot at beginning.'
    if force==False and os.path.exists(out+'.'+extension): 
        print("run_plink: command '{}', the file '{}' already exists (force is set to False)."
             .format(command, out+'.'+extension))
        # ABORT !
        return ('run_plink was aborted (not an error)', '')
    if not stderr :
        # Construct the shell commands
        lst = []
        # Determine plink version to use
        lst.append('plink2' if plink2 else 'plink')
        # Input and output
        lst.append('--bfile'); lst.append(file)
        lst.append('--out');   lst.append(out)
        # Must separate the options from their values (see run_shell_command)
        commands = command.split(' ')
        for c in commands :
            lst.append(c)
        # Provide values to the run_shell_command function
        stdout, stderr = run_shell_command(lst)
    # Display errors
    if stderr and verbose_err : print(stderr)
    # Return the outputs
    return stdout, stderr


def join_position(df):
    pass