import os
import glob


path_in = input("Type the path of files: ")
os.chdir(path_in)
ca_result = [i for i in glob.glob('ca*.txt')]  # only case-insensitive in Windows
print(ca_result)

with open('CA merged.txt', 'w') as outfile:
    with open(ca_result[0]) as f1:
        for line in f1:
            outfile.write(line)
    for ca in ca_result[1:]:
        print(ca)
        with open(ca) as f1:
            for line in f1:
                if not line.startswith('L'):
                    outfile.write(line)
