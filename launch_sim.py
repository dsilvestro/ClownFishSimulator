import os

w_dir = "/Users/danielesilvestro/Software/micro2macroEvolution"

for i in range(10):
	cmd = "cd " + w_dir + "; python ClownSimulator.py -seed " + str(i)
	os.system(cmd)