import os
import glob
# still to find a way to remove the headers of the following files

path_in = input("Type the path of files: ")
os.chdir(path_in)
ca_result = [i for i in glob.glob('ca*.txt')]  # only case-insensitive in Windows
print(ca_result)


with open("result.txt", "wb") as outfile:
    for f in ca_result:
        with open(f, "rb") as infile:
            outfile.write(infile.read())
