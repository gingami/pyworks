import os
import sys
import subprocess
input_file="sample.dimacs"
output_file="output"

#例外処理に変更するかも
if not(os.path.isfile(input_file)):
    print('cannot find "samle.dimacs". ')
    sys.exit()

subprocess.call(["minisat", input_file, output_file])