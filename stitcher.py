t = open("myOut174_Oc.txt","r");
z = open("myOut174_Fr.txt","r");
tempfiles = [t,z]

f = open("myOut174_Wh.txt", "w")
for tempfile in tempfiles:
    f.write(tempfile.read())
