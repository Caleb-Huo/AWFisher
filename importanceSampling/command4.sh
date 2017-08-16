cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'

n0=1e3
n1=1e7

for a in {1..50}
do
let "k=a*2"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling4.R 
done 


cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'

n0=1e3
n1=1e7

for a in {1..49}
do
let "k=a*2 + 1"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling4.R 
done 

cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'
k=2
n1=1000
n0=1000
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling4.R 
