rm -r report
mkdir report
cp *.csv report
cp *.err report
cp *.dat report
#cp *.out report
cp *.ini report
cp submit.sh report
cp main.cpp report
cp hist-iter* report
tar cf out.tar report
