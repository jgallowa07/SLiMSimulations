run:
	rm myOut171.txt
	rm myOut172.txt
	slim MyRecipe1_7_1.E
	slim MyRecipe1_7_2.E
plot1:
	Python3 Hist1_7_1.py
plot2:
	Python3 Hist1_7_2.py
clean:
	rm myOut173_Oc.txt
	rm myOut173_Fr.txt
	
