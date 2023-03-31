compile： 
    g++ -O2 genData.cpp -o ./bin/genData
    g++ -O2 col_pivot_serial.cpp -o ./bin/col_pivot_serial 
    mpicxx -O2 col_pivot_para.cpp -o ./bin/col_pivot_para
run:
    mkdir data
    ./bin/genData 1024
    /bin/col_pivot_serial ./data/Axb_1024.txt
    mpirun -np 4 /bin/col_pivot_para ./data/Axb_1024.txt    