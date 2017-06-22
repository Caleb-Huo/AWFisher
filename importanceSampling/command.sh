cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'

n0=1e3
n1=1e8

for k in {2..10}
do
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling3.R &
done 


cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'

n0=1e3
n1=1e8

for b in {1..10}
do
for a in {1..9}
do	
let "k=a*10+b"  # Increment inner loop counter.
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling3.R 
done &
done 

