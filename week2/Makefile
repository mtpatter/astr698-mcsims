#
#	ASTR698: MC Simulations
#
#	Maria Patterson		21Mar08
#	last updated:	21Mar08
#
TARGET1=mcs01
TARGET2=mcs02
TARGET3=mcs03
TARGET4=mcs04

all:	code2 

code1:	$(TARGET1).c $(TARGET1).o
	gcc -c $(TARGET1).c -o $(TARGET1).o
	gcc -lm $(TARGET1).o -o $(TARGET1)

code2:	$(TARGET2).c $(TARGET2).o
	gcc -c $(TARGET2).c -o $(TARGET2).o -g
	gcc -lm $(TARGET2).o -o $(TARGET2) 

code3:	$(TARGET3).c $(TARGET3).o
	gcc -c $(TARGET3).c -o $(TARGET3).o -g
	gcc -lm $(TARGET3).o -o $(TARGET3) 

code4:	$(TARGET4).c $(TARGET4).o
	gcc -c $(TARGET4).c -o $(TARGET4).o -g
	gcc -lm $(TARGET4).o -o $(TARGET4) 

clean:
	rm -f *.aux *.log *.dvi *.bak *~