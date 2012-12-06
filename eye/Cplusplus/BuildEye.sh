
rm tessar-tessar.o
#rm tessar.exe
#rm layout.svg
#rm spot.svg
#rm longitudinal_abber.svg


g++ -c tessar.cc
mv tessar.o tessar-tessar.o

make

#./tessar


#cygstart spot.svg
#cygstart longitudinal_abber.svg
#cygstart layout.svg
