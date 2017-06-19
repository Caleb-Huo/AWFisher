cd '/mnt/glusterfs/zhh18/AW/AWFisher/importanceSampling'

n0=1e3
n1=1e6
k=2

for k in {2..30}
do
echo "Pass $k to args."
nohup R --no-save --quiet --slave --args $k $n1 $n0  < importanceSampling.R &
done 




