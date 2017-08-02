t = open("myOut176_Oc.txt","r");
z = open("myOut176_Fr.txt","r");
tempfiles = [t,z]

f = open("myOut176_Wh.txt", "w")
for tempfile in tempfiles:
    f.write(tempfile.read())
