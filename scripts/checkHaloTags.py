import re

count = 0
for i in range(0,16):
  file = open("file-" + str(i) + ".txt")
  for line in file.xreadlines():
    line = line.strip()
    if not line.startswith("#"):
      data = re.findall(r'\d+',line)
      print "..." + line.rstrip() + " >> " + data[1]
      if(data[1] < 0):
        print line.rstrip()
        print "TAG: " + data[1]
  file.close()
