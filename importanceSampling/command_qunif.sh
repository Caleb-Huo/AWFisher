cd '/data/home/zhuo/research/Tseng/AW/sysdata'

n0=1e3
n1=1e6

for a in {51..200}
do
let "k=a"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling_qunif.R
done 


cd '/data/home/zhuo/research/Tseng/AW/sysdata'

n0=1e3
n1=1e6

for a in 80 120 180 300 500 1000
do
let "k=a"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling_qunif.R
done 




cd '/data/home/zhuo/research/Tseng/AW/sysdata'

n0=1e3
n1=1e7

for a in {2..100}
do
let "k=a"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling_qunif.R
done 

cd '/data/home/zhuo/research/Tseng/AW/sysdata'
k=2
n1=1000
n0=1000
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling_qunif.R
