compile： <br>
    g++ -O2 genData.cpp -o ./bin/genData<br>
    g++ -O2 col_pivot_serial.cpp -o ./bin/col_pivot_serial <br>
    mpicxx -O2 col_pivot_para.cpp -o ./bin/col_pivot_para<br>
run:<br>
    mkdir data<br>
    ./bin/genData 1024<br>
    /bin/col_pivot_serial ./data/Axb_1024.txt<br>
    mpirun -np 4 /bin/col_pivot_para ./data/Axb_1024.txt    <br>
