rm serial.txt
./base_serial -n 1 -no -s serial.txt
./base_serial -n 10 -no -s serial.txt
./base_serial -n 50 -no -s serial.txt
./base_serial -n 100 -no -s serial.txt
./base_serial -n 250 -no -s serial.txt
./base_serial -n 500 -no -s serial.txt
./base_serial -n 1000 -no -s serial.txt
./base_serial -n 2000 -no -s serial.txt
./base_serial -n 4000 -no -s serial.txt
./base_serial -n 8000 -no -s serial.txt
./base_serial -n 16000 -no -s serial.txt
./base_serial -n 32000 -no -s serial.txt
./autograder -v serial -s serial.txt