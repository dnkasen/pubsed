import os, re

fout = open('all_parameters.rst','w')
fout.write("===================================\n")
fout.write("All Runtime Parameters\n")
fout.write("===================================\n\n")

fout.write("|sedona| uses cgs units everywhere\n")
fout.write("\n|\n\n")




for filename in os.listdir('./'):
	if ('.rst' in filename and filename not in ['all_parameters.rst', 'python_tools.rst']):
		fin = open(filename,'r')

		count = 0
		for line in fin:
			if ('list-table' in line):
				count += 1
			if (count):
				fout.write(line)
			if (line == '\n' and count):
				count += 1
			if (count == 3):
				fout.write("|\n")
				count = 0
		fin.close()
	fout.write('\n')

fout.write('\n')
fout.close()
